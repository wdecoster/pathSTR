# This script is used to serve a Dash app with the information from STR calls of the 1000 Genomes project,
# focusing on the pathogenic STRs as genotyped by STRdust.


import pandas as pd
import dash
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import dash_daq as daq


from argparse import ArgumentParser
import pathSTR.parse_input as parse
from pathSTR.repeats import Repeats
import pathSTR.plot as plot
from pathSTR.version import __version__
from pathSTR.count_kmers import parse_kmers
import os
import sys
import logging
from math import floor, ceil
import dash_bio as dashbio


def main():
    args = get_args()
    # Set up logging, with time stamps, to a file with max size of 10MB and simultaneously to stderr
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler("pathSTR-1000G-dash.log", mode="w"),
            logging.StreamHandler(sys.stderr),
        ],
    )
    logging.info("Starting pathSTR-1000G-dash")
    if args.db:
        df = pd.read_hdf(args.db, key="df").sort_values("gene")
        kmers = {
            (dataset, gene): pd.read_hdf(args.db, key=f"kmer_{dataset}_{gene}")
            for gene in df["gene"].unique()
            for dataset in df["dataset"].unique()
        }
        stat = pd.read_hdf(args.db, key="stat")
        repeats = Repeats(df=pd.read_hdf(args.db, key="repeats"))
        detail_df = pd.read_hdf(args.db, key="details")
        db_version = pd.read_hdf(args.db, key="version").values[0]
        logging.info("Finished reading pathSTR_db file.")
    elif args.vcf and args.sample_info:
        # read in the BED file with the STRs from STRchive
        repeats = Repeats()
        logging.info("Finished parsing repeats.")

        # parse the VCF and sample info file. The VCF argument is a FOFN with columns for vcf, build and genotyper
        df = parse.parse_input(args.vcf, args.sample_info, repeats).sort_values("gene")
        logging.info("Finished parsing VCFs.")

        # Calculate mean and standard deviation per repeat for the comparison with uploaded data
        stat = parse.stats(df)
        logging.info("Finished calculating stats.")

        # parse the kmers for each gene, keeping datasets separate
        # this results in kmers being a dictionary, with the (dataset, gene) tuple as key and the kmer df as value
        # kmer counts cannot trivially be stored in a single DataFrame, as the number and names of columns would be variable
        kmers = {
            (dataset, gene): parse_kmers(df, repeats, gene, dataset)
            for gene in df["gene"].unique()
            for dataset in df["dataset"].unique()
        }
        logging.info("Finished parsing kmers.")

        detail_df = parse.create_details_table(df, repeats)
        logging.info("Finished creating details table.")

        if args.save_db and not args.force:
            if os.path.exists(args.save_db):
                logging.warning(
                    f"pathSTR_db file {args.save_db} already exists, skipping."
                )
            else:
                df.to_hdf(args.save_db, key="df", mode="w")
                stat.to_hdf(args.save_db, key="stat", mode="a")
                for (dataset, gene), kmer_df in kmers.items():
                    kmer_df.to_hdf(args.save_db, key=f"kmer_{dataset}_{gene}", mode="a")
                repeats.df.to_hdf(args.save_db, key="repeats", mode="a")
                detail_df.to_hdf(args.save_db, key="details", mode="a")
                # use the date of today as the database version identifier and save it to the database as a pd.Series
                import datetime

                pd.Series(datetime.date.today().strftime("%Y-%m-%d")).to_hdf(
                    args.save_db, key="version", mode="a"
                )
                logging.info(f"Saved parsed data to {args.save_db}.")
    else:
        logging.error("Provide a database, or a fofn of VCFs and sample info file.")
        raise ValueError("Provide a database, or a fofn of VCFs and sample info file")
    if args.store_only:
        return

    num_samples = len(df["sample"].unique())
    # Create Dash app
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    # Define app layout
    gene_options = [
        {"label": gene, "value": gene} for gene in sorted(df["gene"].unique().tolist())
    ]
    dataset_options = [
        {"label": dataset, "value": dataset}
        for dataset in df["dataset"].unique().tolist()
    ]
    app.layout = html.Div(
        [
            html.Div(
                [
                    html.H1(
                        "Medically-relevant tandem repeats in nanopore sequencing of control cohorts"
                    ),
                ],
                style={
                    "backgroundColor": "#f9f9f9",
                    "border": "1px solid #ddd",
                    "textAlign": "center",
                },
            ),
            dcc.Tabs(
                [
                    dcc.Tab(
                        label="Overview",
                        children=[
                            html.H1("Overview", style={"bottommargin": "0px"}),
                            html.Img(src=dash.get_asset_url("overview-strip.png")),
                            dcc.Loading(
                                id="loading-strip-2",
                                type="cube",
                                children=[
                                    html.Div(
                                        id="strip-plot-log-container",
                                        style={"margin": "0px"},
                                    ),
                                ],
                                style={"margin": "0px", "height": "20vh"},
                            ),
                            html.Br(),
                            html.H1("About"),
                            html.Div(
                                [
                                    html.P(
                                        [
                                            "Welcome to the pathSTR web app, a tool to visualize and compare STR genotypes from nanopore sequencing in control cohorts. ",
                                            "The STR genotypes have been obtained using ",
                                            html.A(
                                                "STRdust",
                                                href="https://github.com/wdecoster/STRdust",
                                                target="_blank",
                                            ),
                                            " and ",
                                            html.A(
                                                "LongTR",
                                                href="https://github.com/gymrek-lab/LongTR",
                                                target="_blank",
                                            ),
                                            " from samples of the 1000 Genomes project, sequenced on the Oxford Nanopore Technologies PromethION. "
                                            "Please let me know if you know a suitable dataset of control individuals to add, for which long read alignments are available online. ",
                                            "If this resource is useful to you, please cite our ",
                                            html.A(
                                                "publication",
                                                href="https://www.medrxiv.org/content/10.1101/2024.03.06.24303700v1/",
                                                target="_blank",
                                            ),
                                            ", as well as the references to the underlying datasets: ",
                                            html.A(
                                                "Noyvert et al. 2023",
                                                href="https://www.medrxiv.org/content/10.1101/2023.12.20.23300308v1",
                                                target="_blank",
                                            ),
                                            ", ",
                                            html.A(
                                                "Schloissnig et al. 2024",
                                                href="https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1",
                                                target="_blank",
                                            ),
                                            " and ",
                                            html.A(
                                                "Gustafson et al. 2024",
                                                href="https://www.medrxiv.org/content/10.1101/2024.03.05.24303792v1",
                                                target="_blank",
                                            ),
                                            ".",
                                        ],
                                        style={"textAlign": "justify"},
                                    ),
                                    html.P(
                                        [
                                            "You can upload your own VCF.gz file(s) (genotyped with STRdust or LongTR) to show alongside the 1000 Genomes data for comparison. Uploading data is currently limited to 100kb files, please let me know if more would be required. For genotyping with STRdust, use the --pathogenic flag. Note that as the app cannot figure out the sex of your individuals, males will show two alleles on haploid chromosomes. This is corrected for samples in the 1000 Genomes cohort. "
                                            "The repeat coordinates and motifs used in this app are obtained from ",
                                            html.A(
                                                "STRchive",
                                                href="https://harrietdashnow.com/STRchive/",
                                                target="_blank",
                                            ),
                                            ". Other dependencies are Python and the Dash/plotly, pandas, hdf5 and cyvcf2 modules for the web app, and snakemake to orchestrate the variant calling.",
                                        ],
                                        style={"textAlign": "justify"},
                                    ),
                                    html.P(
                                        [
                                            f"This web app is developed and maintained by Wouter De Coster. The hosting and deployment is arranged by Svenn D'Hert, as well as some nice layout fixes. The current app version is v{__version__}, and the database of {num_samples} individuals was generated on {db_version}. ",
                                            "The source code is available on ",
                                            html.A(
                                                "GitHub",
                                                href="https://github.com/wdecoster/pathSTR",
                                                target="_blank",
                                            ),
                                            ". Feedback is welcome in the form of an ",
                                            html.A(
                                                "issue on GitHub",
                                                href="https://github.com/wdecoster/pathSTR/issues",
                                                target="_blank",
                                            ),
                                            " or by sending me an email. The repository also contains documentation for the aSTRonaut companion script, for a more flexible kmer visualization in the ",
                                            html.I("sequence"),
                                            " mode.",
                                        ],
                                        style={"textAlign": "justify"},
                                    ),
                                ],
                                style={"width": "80%", "margin": "auto"},
                            ),
                            html.H1("Options", className="my-3"),
                            dcc.Checklist(
                                id="publication-ready",
                                options=[
                                    {
                                        "label": "Publication-ready figures",
                                        "value": "publication-ready",
                                    }
                                ],
                                value=[],
                                inline=True,
                                inputStyle={"margin-left": "15px"},
                            ),
                            # add a checkbox for making the strip plot dynamic
                            dcc.Checklist(
                                id="strip-dynamic",
                                options=[
                                    {
                                        "label": "Dynamic strip plot [rather slow]",
                                        "value": "dynamic",
                                    }
                                ],
                                value=[],
                                inline=True,
                                inputStyle={"margin-left": "15px"},
                            ),
                            # add a dropdown for selecting the dataset, with default STRdust
                            dbc.Row(
                                [
                                    dbc.Col(
                                        html.Label(
                                            "Select dataset for pathSTR:",
                                            style={"margin-left": "15px"},
                                        ),
                                        width=2,
                                    ),
                                    dbc.Col(
                                        dcc.Dropdown(
                                            id="dropdown-dataset",
                                            options=dataset_options,
                                            value="STRdust_hg38",
                                            clearable=False,
                                        ),
                                        width=2,
                                    ),
                                ],
                                align="center",
                            ),
                            html.H1("Downloads", className="my-3"),
                            html.Div(
                                [
                                    html.Div(
                                        children=[
                                            html.Button(
                                                "Download Data as TSV",
                                                id="btn",
                                                className="btn btn-primary btn-lg",
                                            ),
                                            html.Button(
                                                "Download VCFs",
                                                id="download-zip-button",
                                                className="btn btn-success btn-lg mx-3",
                                            ),
                                        ],
                                        style={
                                            "display": "flex",
                                            "justifyContent": "flex-start",
                                        },
                                    ),
                                    dcc.Download(id="download"),
                                    dcc.Download(id="download-zip"),
                                ],
                            ),
                            html.H1("Repeats", className="my-3"),
                            html.Div(
                                [
                                    html.P(
                                        "Repeats used in this app are obtained from STRchive, genotyped for the coordinates below:"
                                    ),
                                    html.Div(
                                        dash_table.DataTable(
                                            repeats.df[
                                                [
                                                    "chrom",
                                                    "start",
                                                    "end",
                                                    "name",
                                                ]
                                            ].to_dict("records"),
                                            [
                                                {"name": i, "id": i}
                                                for i in [
                                                    "chrom",
                                                    "start",
                                                    "end",
                                                    "name",
                                                ]
                                            ],
                                            style_data={"user-select": "auto"},
                                            export_format="csv",
                                        ),
                                        className="table",
                                        style={
                                            "width": "30%",
                                            "margin": "auto",
                                            "display": "flex",
                                            "justifyContent": "flex-start",
                                        },
                                    ),
                                ],
                                style={
                                    "display": "flex",
                                    "justifyContent": "flex-start",
                                },
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Repeat length",
                        children=[
                            html.H1("Repeat length"),
                            html.Div(
                                children=[
                                    dcc.Dropdown(
                                        id="dropdown-gene-length",
                                        options=gene_options,
                                        value=gene_options[0]["value"],
                                        clearable=False,
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                dcc.Checklist(
                                                    id="violin_options",
                                                    options=[
                                                        {
                                                            "label": "Split by population",
                                                            "value": "population",
                                                        },
                                                        {
                                                            "label": "Split by sex",
                                                            "value": "sex",
                                                        },
                                                        {
                                                            "label": "Show repeat length relative to reference genome",
                                                            "value": "ref_diff",
                                                        },
                                                        {
                                                            "label": "Show log-transformed length",
                                                            "value": "log",
                                                        },
                                                        {
                                                            "label": "Show pathogenic length",
                                                            "value": "pathlen",
                                                        },
                                                        {
                                                            "label": "Show density",
                                                            "value": "density",
                                                        },
                                                    ],
                                                    value=["density"],
                                                    inline=True,
                                                    inputStyle={"margin-left": "15px"},
                                                )
                                            ),
                                            dbc.Col(
                                                html.Div(id="warning-pathogenic-length")
                                            ),
                                        ]
                                    ),
                                ]
                            ),
                            html.Div(dcc.Graph(id="violin-plot")),
                            html.Div(
                                dcc.Graph(
                                    id="length-scatter",
                                    style={
                                        "display": "flex",
                                        "justifyContent": "center",
                                        "alignItems": "center",
                                    },
                                ),
                                style={
                                    "justifyContent": "center",
                                    "alignItems": "center",
                                },
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Repeat Composition",
                        children=[
                            html.H1("Repeat Composition"),
                            dbc.Container(
                                [
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label(
                                                    "Select visualization mode:"
                                                ),
                                                width=3,
                                            ),
                                            dbc.Col(
                                                dcc.RadioItems(
                                                    id="kmer_mode",
                                                    options=[
                                                        {
                                                            "label": "Collapsed",
                                                            "value": "collapsed",
                                                        },
                                                        {
                                                            "label": "Raw",
                                                            "value": "raw",
                                                        },
                                                        {
                                                            "label": "Sequence",
                                                            "value": "sequence",
                                                        },
                                                    ],
                                                    value="collapsed",
                                                    inline=True,
                                                    style={"display": "flex"},
                                                ),
                                                width=3,
                                            ),
                                        ],
                                        align="center",
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Select gene:"), width=3
                                            ),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="dropdown-gene-composition",
                                                    options=gene_options,
                                                    value=gene_options[0]["value"],
                                                    clearable=False,
                                                ),
                                                width=2,
                                            ),
                                        ],
                                        align="center",
                                    ),
                                ]
                            ),
                            dbc.Container(
                                [
                                    dbc.Row(
                                        [
                                            html.Div(
                                                id="kmer-options-container",
                                                children=[
                                                    html.Div(
                                                        id="kmer-options-collapsed",
                                                        children=[
                                                            dbc.Row(
                                                                [
                                                                    dbc.Col(
                                                                        html.Label(
                                                                            "Minimal group size:"
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                    dbc.Col(
                                                                        daq.NumericInput(
                                                                            id="kmer-options-group-min-size",
                                                                            value=0,
                                                                            max=1000,
                                                                        ),
                                                                        width="auto",
                                                                    ),
                                                                ],
                                                                align="center",
                                                            )
                                                        ],
                                                        style={"display": "block"},
                                                    ),
                                                    html.Div(
                                                        id="kmer-options-sequence",
                                                        children=[
                                                            dbc.Row(
                                                                [
                                                                    dbc.Col(
                                                                        html.Label(
                                                                            "Direction:"
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                    dbc.Col(
                                                                        dcc.Dropdown(
                                                                            id="kmer-options-sequence-direction",
                                                                            options=[
                                                                                "left-to-right",
                                                                                "right-to-left",
                                                                            ],
                                                                            multi=False,
                                                                            value="left-to-right",
                                                                            placeholder="",
                                                                        ),
                                                                        width=2,
                                                                    ),
                                                                ],
                                                                align="center",
                                                            ),
                                                            dbc.Row(
                                                                [
                                                                    dbc.Col(
                                                                        html.Label(
                                                                            "Show pathogenic length:"
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                    dbc.Col(
                                                                        dcc.Dropdown(
                                                                            id="kmer-options-sequence-pathlen",
                                                                            options=[
                                                                                "on",
                                                                                "off",
                                                                            ],
                                                                            value="off",
                                                                            clearable=False,
                                                                        ),
                                                                        width=2,
                                                                    ),
                                                                ],
                                                                align="center",
                                                            ),
                                                        ],
                                                        style={"display": "none"},
                                                    ),
                                                    html.Div(
                                                        id="kmer-options-raw",
                                                        children=[
                                                            dbc.Row(
                                                                [
                                                                    dbc.Col(
                                                                        html.Label(
                                                                            "Kmer to sort heatmap for:"
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                    dbc.Col(
                                                                        dcc.Dropdown(
                                                                            id="kmer-options-raw-sort",
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                ],
                                                                align="center",
                                                            ),
                                                            dbc.Row(
                                                                [
                                                                    dbc.Col(
                                                                        html.Label(
                                                                            "Number of columns to split the heatmap on:"
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                    dbc.Col(
                                                                        daq.NumericInput(
                                                                            id="kmer-options-raw-num-columns",
                                                                            value="",
                                                                        ),
                                                                        width=3,
                                                                    ),
                                                                ],
                                                                align="center",
                                                            ),
                                                        ],
                                                        style={"display": "none"},
                                                    ),
                                                ],
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        html.Label(
                                                            "Minimal repeat length [in units]:"
                                                        ),
                                                        width=3,
                                                    ),
                                                    dbc.Col(
                                                        dcc.RangeSlider(
                                                            id="repeat-len-slider",
                                                        ),
                                                        width=8,
                                                    ),
                                                ],
                                                align="center",
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                children=[
                                    html.Div(
                                        id="kmer-composition-raw-container",
                                        children=[dcc.Graph(id="kmer-composition-raw")],
                                        style={"display": "none"},
                                    ),
                                    html.Div(
                                        id="kmer-composition-sequence-container",
                                        children=[
                                            dcc.Graph(id="kmer-composition-sequence")
                                        ],
                                        style={"display": "none"},
                                    ),
                                    html.Div(
                                        id="kmer-composition-collapsed-container",
                                        children=[
                                            dcc.Graph(id="kmer-composition-collapsed")
                                        ],
                                        style={"display": "block"},
                                    ),
                                ]
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Details per individual",
                        children=[
                            html.H1("Details per individual and repeat"),
                            dbc.Container(
                                [
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.Label("Select individual:"),
                                                    dcc.Dropdown(
                                                        id="dropdown-details-individual",
                                                        options=[
                                                            {
                                                                "label": individual,
                                                                "value": individual,
                                                            }
                                                            for individual in df[
                                                                "sample"
                                                            ].unique()
                                                        ],
                                                        value=df["sample"].unique()[0],
                                                        clearable=False,
                                                        multi=True,
                                                    ),
                                                ]
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Label("Select gene:"),
                                                    dcc.Dropdown(
                                                        id="dropdown-details-gene",
                                                        options=gene_options,
                                                        value=gene_options[0]["value"],
                                                        clearable=False,
                                                    ),
                                                ]
                                            ),
                                        ]
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    dash_table.DataTable(
                                                        id="details-table",
                                                        style_cell={
                                                            "fontSize": 20,
                                                            "font-family": "arial",
                                                            "textAlign": "left",
                                                            "overflow": "hidden",
                                                            "textOverflow": "ellipsis",
                                                            "minWidth": "180px",
                                                            "width": "180px",
                                                            "maxWidth": "180px",
                                                        },
                                                        style_header={
                                                            "backgroundColor": "white",
                                                            "fontWeight": "bold",
                                                            "font-family": "sans-serif",
                                                            "fontSize": 24,
                                                        },
                                                        tooltip_duration=None,
                                                    ),
                                                ],
                                            ),
                                        ],
                                        className="py-4",
                                    ),
                                ]
                            ),
                            dbc.Button(
                                "Show in IGV",
                                id="igv-button",
                                color="primary",
                                className="mt-2",
                            ),
                            dcc.Loading(
                                id="igv-output",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Upload your data",
                        children=[
                            html.H1("Upload your data"),
                            dcc.Store(id="stored-df"),
                            dcc.Upload(
                                id="upload-data",
                                children=html.Div(
                                    [
                                        "Drag and drop or ",
                                        html.A(
                                            "click to upload VCF.gz files (hg38) with STR genotypes"
                                        ),
                                        " to show your data in the plots",
                                    ]
                                ),
                                style={
                                    "width": "100%",
                                    "height": "60px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "margin": "10px",
                                },
                                multiple=True,
                                max_size=100000,
                            ),
                            html.Div(id="upload-status"),
                            dash_table.DataTable(
                                id="user-data-table",
                                sort_action="native",
                                filter_action="native",
                                row_deletable=True,
                                page_action="native",
                                page_size=50,
                                style_cell={
                                    "fontSize": 14,
                                    "font-family": "sans-serif",
                                },
                                style_header={
                                    "backgroundColor": "white",
                                    "fontWeight": "bold",
                                    "font-family": "sans-serif",
                                    "fontSize": 18,
                                },
                                style_data_conditional=[
                                    {
                                        "if": {"row_index": "odd"},
                                        "backgroundColor": "rgb(248, 248, 248)",
                                    },
                                    {
                                        "if": {
                                            "column_id": "max_z_score",
                                            "filter_query": "{max_z_score} > 3",
                                        },
                                        "fontWeight": "bold",
                                        "color": "rgd(155, 0, 0)",
                                    },
                                ],
                                export_format="xlsx",
                                export_headers="display",
                            ),
                        ],
                    ),
                ]
            ),
        ]
    )

    @app.callback(
        Output("download", "data"), Input("btn", "n_clicks"), prevent_initial_call=True
    )
    def generate_csv(n_clicks):
        """
        Enable users to download the data as a TSV file
        """
        return dcc.send_data_frame(df.to_csv, "pathSTR-1000G.tsv", sep="\t")

    @app.callback(
        [Output("stored-df", "data"), Output("upload-status", "children")],
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
        State("upload-data", "last_modified"),
        State("dropdown-dataset", "value"),
    )
    def store_uploaded_data(
        list_of_contents, list_of_filenames, list_of_dates, dataset
    ):
        """
        Process the uploaded user-VCF files
        """
        if list_of_contents is not None:
            caller, build = dataset.split("_")
            dfs = [
                parse.parse_uploaded_vcf(content, filename, repeats, build, caller)
                for (content, filename, _) in zip(
                    list_of_contents, list_of_filenames, list_of_dates
                )
            ]
            for df, filename in zip(dfs, list_of_filenames):
                if df is None:
                    logging.warning(
                        f"Error parsing {filename}, please make sure it is a valid VCF.gz file generated by the specified genotyper"
                    )
                    # Create a dummy dataframe to communicate error
                    return (
                        (pd.DataFrame(()).to_dict("records")),
                        f"Error processing file {filename}, please make sure it is a valid VCF.gz file, generated by the specified genotyper for the STRchive STRs. Please report unsolved issues on GitHub at https://github.com/wdecoster/pathSTR/issues/.",
                    )
            uploaded_df = pd.concat(dfs, ignore_index=True)
            plural = "s" if len(dfs) > 1 else ""
            return uploaded_df.to_dict("records"), f"Uploaded {len(dfs)} file{plural}"
        else:
            return pd.DataFrame(()).to_dict("records"), ""

    @app.callback(
        [
            Output("violin-plot", "figure"),
            Output("warning-pathogenic-length", "children"),
        ],
        [
            Input("dropdown-gene-length", "value"),
            Input("stored-df", "data"),
            Input("violin_options", "value"),
            Input("publication-ready", "value"),
            Input("dropdown-dataset", "value"),
        ],
    )
    def update_violin(
        selected_gene, stored_df, violin_options, publication_ready, dataset
    ):
        """
        Create repeat length violin plot
        """
        if stored_df is None:
            filtered_df = df[(df["gene"] == selected_gene) & (df["dataset"] == dataset)]
        else:
            stored_df = pd.DataFrame(stored_df)
            combined_df = pd.concat([df, stored_df], ignore_index=True)
            filtered_df = combined_df[
                (combined_df["gene"] == selected_gene)
                & (combined_df["dataset"] == dataset)
            ]
        warning = (
            html.P(
                "Interpret the pathogenic length with caution and always take repeat composition into account.",
                style={
                    "color": "red",
                    "fontSize": 14,
                    "fontStyle": "italic",
                    "margin-left": 15,
                },
            )
            if "pathlen" in violin_options
            else ""
        )
        return (
            plot.violin_plot(
                filtered_df,
                repeats=repeats,  # show optionally the pathogenic length
                selected_gene=selected_gene,
                violin_options=violin_options,
                publication_ready="publication-ready" in publication_ready,
            ),
            warning,
        )

    ## the text in the upload data field should be adapted based on the selected dataset in About
    @app.callback(
        Output("upload-data", "children"),
        Input("dropdown-dataset", "value"),
    )
    def update_upload_data_text(dataset):
        return [
            "Drag and drop or ",
            html.A(
                "click",
                href="#",
            ),
            f" to upload {dataset} VCF.gz files to show your data in the plots",
        ]

    @app.callback(
        Output("length-scatter", "figure"),
        [
            Input("dropdown-gene-length", "value"),
            Input("stored-df", "data"),
            Input("violin_options", "value"),
            Input("publication-ready", "value"),
            Input("dropdown-dataset", "value"),
        ],
    )
    def update_length_scatter(
        selected_gene, stored_df, violin_options, publication_ready, dataset
    ):
        """
        Create repeat length scatter plot
        """
        if stored_df is None:
            filtered_df = df[(df["gene"] == selected_gene) & (df["dataset"] == dataset)]
        else:
            stored_df = pd.DataFrame(stored_df)
            combined_df = pd.concat([df, stored_df], ignore_index=True)
            filtered_df = combined_df[
                (combined_df["gene"] == selected_gene)
                & (combined_df["dataset"] == dataset)
            ]
        return plot.length_scatter(
            filtered_df,
            selected_gene,
            path_length=repeats.pathogenic_min_length(selected_gene),  # optional
            violin_options=violin_options,
            publication_ready="publication-ready" in publication_ready,
        )

    # Create a callback that changes the visibility of the kmer options and kmer-plots based on the mode
    @app.callback(
        [
            Output("kmer-options-raw", "style"),
            Output("kmer-options-sequence", "style"),
            Output("kmer-options-collapsed", "style"),
            Output("kmer-composition-raw-container", "style"),
            Output("kmer-composition-sequence-container", "style"),
            Output("kmer-composition-collapsed-container", "style"),
        ],
        Input("kmer_mode", "value"),
    )
    def show_kmer_options(mode):
        """
        Show the corret kmer options and correct plot based on the selected mode
        I get it, this is rather awful
        """
        if mode == "raw":
            return (
                {"display": "block"},
                {"display": "none"},
                {"display": "none"},
                {
                    "display": "flex",
                    "justifyContent": "center",
                    "alignItems": "center",
                },
                {"display": "none"},
                {"display": "none"},
            )
        elif mode == "sequence":
            return (
                {"display": "none"},
                {"display": "block"},
                {"display": "none"},
                {"display": "none"},
                {"display": "flex", "justifyContent": "center", "alignItems": "center"},
                {"display": "none"},
            )
        else:
            return (
                {"display": "none"},
                {"display": "none"},
                {"display": "block"},
                {"display": "none"},
                {"display": "none"},
                {"display": "flex", "justifyContent": "center", "alignItems": "center"},
            )

    ## changing the dropdown-gene-composition should clear the kmer-options-raw-sort values
    @app.callback(
        Output("kmer-options-raw-sort", "value"),
        Input("dropdown-gene-composition", "value"),
    )
    def clear_kmer_sort(selected_gene):
        return None

    @app.callback(
        [
            Output("kmer-options-raw-sort", "options"),
        ],
        [
            Input("dropdown-gene-composition", "value"),
            Input("dropdown-dataset", "value"),
        ],
    )
    def specify_kmers_to_sort_on(selected_gene, dataset):
        """
        Return the options for the kmer dropdown in the raw mode
        """
        return ([i for i in kmers[(dataset, selected_gene)].columns if i != "length"],)

    @app.callback(
        Output("kmer-composition-raw", "figure"),
        [
            Input("kmer_mode", "value"),
            Input("dropdown-gene-composition", "value"),
            Input("repeat-len-slider", "value"),
            Input("stored-df", "data"),
            Input("kmer-options-raw-sort", "value"),
            Input("kmer-options-raw-num-columns", "value"),
            State("publication-ready", "value"),
            Input("dropdown-dataset", "value"),
        ],
        suppress_callback_exceptions=True,
    )
    def update_kmer_composition_raw(
        kmer_mode,
        selected_gene,
        length_range,
        stored_df,
        raw_sort,
        num_columns,
        publication_ready,
        dataset,
    ):
        """
        Create a kmer composition plot
        :param kmer_mode: mode to show the kmer composition plot in (raw, collapsed or sequence)
        :param selected_gene: gene to show the kmer composition for
        :param length_range: minimal and maximal length of the repeats to show
        :param stored_df: uploaded data
        :param raw_sort: column to sort the raw kmer composition on
        :param publication_ready: whether to show a publication-ready plot
        :param dataset: dataset to show the kmer composition for
        """
        if kmer_mode == "raw":
            if len(stored_df) == 0:
                kmer_df = kmers[(dataset, selected_gene)]
            else:
                stored_df = pd.DataFrame(stored_df)
                kmer_df = pd.concat(
                    [
                        kmers[(dataset, selected_gene)],
                        parse_kmers(stored_df, repeats, selected_gene, dataset),
                    ]
                ).fillna(0.0)
            return plot.kmer_plot_raw(
                kmer_df,
                selected_gene=selected_gene,
                length_range=length_range,
                sort=raw_sort,
                num_columns=num_columns,
                publication_ready="publication-ready" in publication_ready,
            )
        else:
            return dash.no_update

    @app.callback(
        Output("kmer-composition-collapsed", "figure"),
        [
            Input("kmer_mode", "value"),
            Input("dropdown-gene-composition", "value"),
            Input("repeat-len-slider", "value"),
            Input("stored-df", "data"),
            Input("kmer-options-group-min-size", "value"),
            State("publication-ready", "value"),
            Input("dropdown-dataset", "value"),
        ],
        suppress_callback_exceptions=True,
    )
    def update_kmer_composition_collapsed(
        kmer_mode,
        selected_gene,
        length_range,
        stored_df,
        collapsed_group_min_size,
        publication_ready,
        dataset,
    ):
        """
        Create a kmer composition plot
        :param kmer_mode: mode to show the kmer composition plot in (raw, collapsed or sequence)
        :param selected_gene: gene to show the kmer composition for
        :param length_range: minimal and maximal length of the repeats to show
        :param stored_df: uploaded data
        :param collapsed_group_min_size: minimal group size to show in the collapsed mode
        :param raw_sort: column to sort the raw kmer composition on
        :param sequence_direction: direction to show the sequence kmer composition in
        :param publication_ready: whether to show a publication-ready plot
        :param dataset: dataset to show the kmer composition for
        """
        if kmer_mode == "collapsed":
            if len(stored_df) == 0:
                kmer_df = kmers[(dataset, selected_gene)]
            else:
                stored_df = pd.DataFrame(stored_df)
                kmer_df = pd.concat(
                    [
                        kmers[(dataset, selected_gene)],
                        parse_kmers(stored_df, repeats, selected_gene, dataset),
                    ]
                ).fillna(0.0)
            return plot.kmer_plot_collapsed(
                kmer_df,
                selected_gene=selected_gene,
                length_range=length_range,
                min_group_size=collapsed_group_min_size,
                publication_ready="publication-ready" in publication_ready,
            )
        else:
            return dash.no_update

    @app.callback(
        Output("kmer-composition-sequence", "figure"),
        [
            Input("kmer_mode", "value"),
            Input("dropdown-gene-composition", "value"),
            Input("repeat-len-slider", "value"),
            Input("stored-df", "data"),
            Input("kmer-options-sequence-direction", "value"),
            Input("kmer-options-sequence-pathlen", "value"),
            State("publication-ready", "value"),
            Input("dropdown-dataset", "value"),
        ],
        suppress_callback_exceptions=True,
    )
    def update_kmer_composition_sequence(
        kmer_mode,
        selected_gene,
        length_range,
        stored_df,
        sequence_direction,
        show_pathogenic_length,
        publication_ready,
        dataset,
    ):
        """
        Create a kmer composition plot
        :param kmer_mode: mode to show the kmer composition plot in (raw, collapsed or sequence)
        :param selected_gene: gene to show the kmer composition for
        :param length_range: minimal and maximal length of the repeats to show
        :param stored_df: uploaded data
        :param sequence_direction: direction to show the sequence kmer composition in
        :param publication_ready: whether to show a publication-ready plot
        :param dataset: dataset to show the kmer composition for
        """
        if kmer_mode == "sequence":
            if len(stored_df) == 0:
                kmer_df = kmers[(dataset, selected_gene)]
            else:
                stored_df = pd.DataFrame(stored_df)
                kmer_df = pd.concat(
                    [
                        kmers[(dataset, selected_gene)],
                        parse_kmers(stored_df, repeats, selected_gene, dataset),
                    ]
                ).fillna(0.0)
            filtered_df = df[(df["gene"] == selected_gene) & (df["dataset"] == dataset)]
            # pathogenic lenght is in units of the motif length, so has to be converted to basepairs for this plot
            pathogenic_length = (
                repeats.pathogenic_min_length(selected_gene)
                * repeats.motif_length(selected_gene)
                if show_pathogenic_length == "on"
                else None
            )
            return plot.kmer_plot_sequence(
                kmer_df,
                repeat_df=filtered_df,
                selected_gene=selected_gene,
                pathogenic_length=pathogenic_length,
                length_range=length_range,
                direction=sequence_direction,
                publication_ready="publication-ready" in publication_ready,
            )
        else:
            return dash.no_update

    @app.callback(
        Output("repeat-len-slider", "min"),
        Output("repeat-len-slider", "max"),
        Input("dropdown-gene-composition", "value"),
        Input("dropdown-dataset", "value"),
    )
    def update_slider_range(selected_gene, dataset):
        """Update the slider based on the minimal and maximal length of the repeat"""
        data = kmers[(dataset, selected_gene)].dropna(subset=["length"])
        return floor(data["length"].min()), ceil(data["length"].max())

    @app.callback(
        Output("strip-plot-log-container", "children"),
        [
            Input("stored-df", "data"),
            Input("strip-dynamic", "value"),
            Input("dropdown-dataset", "value"),
        ],
    )
    def update_stripplot(stored_df, dynamic, dataset):
        """
        Create a strip plot of the repeat lengths
        By default, only a static image is shown
        If the dynamic checkbox is checked, the plot is updated with the uploaded data
        This makes the app quite a bit faster
        """
        if dynamic:
            if stored_df is None:
                return dcc.Graph(
                    id="strip-plot-log",
                    figure=plot.create_strip_plot(
                        df[df["dataset"] == dataset], log=True
                    ),
                )
            else:
                strip_df = pd.concat([df, pd.DataFrame(stored_df)], ignore_index=True)
                return dcc.Graph(
                    id="strip-plot-log",
                    figure=plot.create_strip_plot(
                        strip_df[strip_df["dataset"] == dataset], log=True
                    ),
                )
        else:
            dash.no_update

    @app.callback(
        [
            Output("details-table", "data"),
            Output("details-table", "columns"),
            Output("details-table", "tooltip_data"),
        ],
        [
            Input("dropdown-details-individual", "value"),
            Input("dropdown-details-gene", "value"),
            Input("dropdown-dataset", "value"),
        ],
    )
    def update_details_table(individuals, gene, dataset):
        """
        Update the details table based on the selected individual and gene
        """
        if not individuals:
            return (
                dash.no_update,
                dash.no_update,
                dash.no_update,
            )
        # if only one individual is selected, it is not a list
        if isinstance(individuals, str):
            individuals = [individuals]

        details_table = (
            detail_df[
                (detail_df["dataset"] == dataset)
                & (detail_df["gene"] == gene)
                & (detail_df.index.isin(individuals))
            ]
            .drop(columns=["dataset", "gene"], level=0)
            .transpose()
        )
        # convert the multi-index to a single column (with empty name)
        details_table.insert(
            0, "", [f"{j}:{i}" if j else i for i, j in details_table.index]
        )
        columns = [({"name": c, "id": c}) for c in details_table.columns]
        details = details_table.to_dict("records")
        tooltip_data = [
            {
                column: {"value": str(value), "type": "markdown"}
                for column, value in row.items()
            }
            for row in details
        ]
        return details, columns, tooltip_data

    @app.callback(
        Output("igv-output", "children"),
        Input("igv-button", "n_clicks"),
        [
            State("dropdown-details-individual", "value"),
            State("dropdown-details-gene", "value"),
            Input("dropdown-dataset", "value"),
        ],
    )
    def return_igv(n_clicks, individuals, gene, dataset):
        """
        When the button is clicked, show the IGV plot for the selected individual and gene
        """
        if n_clicks is not None:
            if not individuals:
                return html.Div()
            if isinstance(individuals, str):
                individuals = [individuals]
            chrom, start, end = repeats.coords(gene)
            build = dataset.split("_")[1]
            return html.Div(
                [
                    dashbio.Igv(
                        id="igv",
                        genome="hg38_1kg" if build == "hg38" else "hs1",
                        locus=f"{chrom}:{start-25}-{end+25}",
                        tracks=(
                            [
                                make_igv_alignment_track(individual, build)
                                for individual in individuals
                            ]
                        ),
                    )
                ]
            )

    def make_igv_alignment_track(individual, build):
        # use the df to get the hg38 url from the individual (sample)
        if not build == "hg38":
            sys.stderr.write("Only hg38 alignments are supported - FIXME\n")
        hg38_path = df.loc[df["sample"] == individual, "hg38_path"].values[0]
        alignment_type = hg38_path.split(".")[-1]
        index_extension = "crai" if alignment_type == "cram" else "bai"
        return {
            "name": individual,
            "type": "alignment",
            "url": hg38_path,
            "indexURL": f"{hg38_path}.{index_extension}",
            "height": 300,
            "format": alignment_type,
            "showSoftClips": True,
            "showInsertionText": True,
            "showDeletionText": True,
            "showSoftClips": True,
            "colorBy": "strand",
        }

    @app.callback(
        Output("download-zip", "data"),
        [Input("download-zip-button", "n_clicks"), State("dropdown-dataset", "value")],
        prevent_initial_call=True,
    )
    def download_zip(n_clicks, dataset):
        caller = dataset.split("_")[0]
        return dcc.send_file(f"pathSTR_{caller}_good_samples.zip")

    @app.callback(
        Output("user-data-table", "data"),
        Output("user-data-table", "columns"),
        Input("stored-df", "data"),
        State("dropdown-dataset", "value"),
    )
    def update_table(data, dataset):
        if data is None or len(data) == 0:
            return (
                dash.no_update,
                dash.no_update,
            )  # Don't update the table if no data is uploaded
        user_df = (
            pd.DataFrame(data)
            .drop(columns=["chrom", "sequence", "Superpopulation", "Group", "Sex"])
            .set_index("gene", drop=False)
            .join(stat.loc[dataset])
        )
        # Calculate z-scores for the uploaded data
        user_df["z_score"] = (
            (user_df["length"] - user_df["mean"]) / user_df["std"]
        ).round(1)
        # pivot the data to show a row per gene and sample - not per allele
        user_df = user_df.pivot(
            index=["gene", "sample"],
            columns=["allele"],
            values=["length", "z_score", "mean", "std"],
        )
        # Flatten the multi-index columns
        user_df.columns = [
            "_".join([str(v) for v in c]) for c in user_df.columns.to_flat_index()
        ]

        # Calculate a max z-score per row for ease of sorting
        user_df["max_z_score"] = user_df[["z_score_Allele1", "z_score_Allele2"]].max(
            axis=1
        )
        # drop duplicated mean and std columns and rename remaining
        user_df = user_df.drop(columns=["mean_Allele2", "std_Allele2"]).rename(
            columns={
                "mean_Allele1": "1000G mean",
                "std_Allele1": "1000G std",
            }
        )
        # reorder the columns
        user_df = user_df.reset_index().rename(
            columns={
                "sample": "individual",
            }
        )[
            [
                "individual",
                "gene",
                "length_Allele1",
                "z_score_Allele1",
                "length_Allele2",
                "z_score_Allele2",
                "max_z_score",
                "1000G mean",
                "1000G std",
            ]
        ]
        # make numeric columns numeric in the datatable
        columns = [
            (
                {"name": c, "id": c, "deletable": True}
                if c in ["individual", "gene"]
                else {"name": c, "id": c, "deletable": True, "type": "numeric"}
            )
            for c in user_df.columns
        ]

        return user_df.to_dict("records"), columns

    app.title = "pathSTR"
    # Run the app
    app.run(host=args.host, port=args.port, debug=args.debug)


def get_args():
    parser = ArgumentParser(description="Get STRs from VCF")
    parser.add_argument(
        "--vcf",
        help="Fofn of Input VCF(s) with a column for build and genotyper, without header",
    )
    parser.add_argument(
        "--sample_info",
        help="Sample info file created by the workflow (pathSTR_samples.tsv)",
    )
    parser.add_argument("--db", help="Input is in one pathSTR_db file")
    parser.add_argument("--save_db", help="Save the parsed data to a pathSTR_db file")
    parser.add_argument(
        "--store_only",
        help="Only store the parsed data and do not start the Dash app",
        action="store_true",
    )
    parser.add_argument(
        "--force", help="Overwrite the existing pathSTR_db file", action="store_true"
    )
    parser.add_argument(
        "--debug",
        help="Run the app in debug mode",
        action="store_true",
    )
    parser.add_argument(
        "--host", help="Host IP used to serve the application", default="127.0.0.1"
    )
    parser.add_argument(
        "--port", help="Port used to serve the application", default=8050, type=int
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
