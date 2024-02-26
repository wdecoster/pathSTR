# This script is used to serve a Dash app with the information from STR calls of the 1000 Genomes project,
# focusing on the pathogenic STRs as genotyped by STRdust.


import pandas as pd
import dash
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from argparse import ArgumentParser
import pathSTR.parse_input as parse
from pathSTR.repeats import Repeats
import pathSTR.plot as plot
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
            gene: pd.read_hdf(args.db, key=f"kmer_{gene}")
            for gene in df["gene"].unique()
        }
        stat = pd.read_hdf(args.db, key="stat")
        repeats = Repeats(df=pd.read_hdf(args.db, key="repeats"))
        logging.info("Finished reading pathSTR_db file.")
    elif args.vcf and args.sample_info:
        # read in the BED file with the STRs from STRchive
        repeats = Repeats()
        logging.info("Finished parsing repeats.")
        # parse the VCF and sample info file
        df = parse.parse_input(args.vcf, args.sample_info, repeats).sort_values("gene")
        logging.info("Finished parsing --vcf.")
        # Calculate mean and standard deviation per repeat for the comparison with uploaded data
        stat = parse.stats(df)
        logging.info("Finished calculating stats.")
        kmers = {gene: parse_kmers(df, repeats, gene) for gene in df["gene"].unique()}

        logging.info("Finished parsing kmers.")
        if args.save_db:
            if os.path.exists(args.save_db):
                logging.warning("pathSTR_db file already exists, skipping.")
            else:
                df.to_hdf(args.save_db, key="df", mode="w")
                stat.to_hdf(args.save_db, key="stat", mode="a")
                for gene, kmer_df in kmers.items():
                    kmer_df.to_hdf(args.save_db, key=f"kmer_{gene}", mode="a")
                repeats.df.to_hdf(args.save_db, key="repeats", mode="a")
                logging.info(f"Saved parsed data to {args.save_db}.")
    else:
        logging.error(
            "Provide either a pathSTR_db file, or VCF, sample info and BED file."
        )
        raise ValueError(
            "Provide either a pathSTR_db file, or VCF, sample info and BED file"
        )
    if args.store_only:
        return
    # Create Dash app
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    # Define app layout
    gene_options = [
        {"label": gene, "value": gene} for gene in sorted(df["gene"].unique().tolist())
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
                            dcc.Store(id="strip-plot-log-store"),
                            dcc.Loading(
                                id="loading-strip-2",
                                type="cube",
                                children=[
                                    html.Div(
                                        dcc.Graph(id="strip-plot-log"),
                                        style={"margin": "0px"},
                                    ),
                                ],
                                style={"margin": "0px", "height": "20vh"},
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
                                    dcc.Checklist(
                                        id="violin_options",
                                        options=[
                                            {
                                                "label": "Split by population",
                                                "value": "population",
                                            },
                                            {"label": "Split by sex", "value": "sex"},
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
                                                dcc.Dropdown(
                                                    id="dropdown-gene-composition",
                                                    options=gene_options,
                                                    value=gene_options[0]["value"],
                                                    clearable=False,
                                                ),
                                            ),
                                            dbc.Col(
                                                dcc.Dropdown(id="dropdown-kmer-options")
                                            ),
                                        ]
                                    )
                                ]
                            ),
                            dbc.Container(
                                [
                                    dbc.Row(
                                        [
                                            dcc.RadioItems(
                                                id="kmer_mode",
                                                options=[
                                                    {
                                                        "label": "Collapsed",
                                                        "value": "collapsed",
                                                    },
                                                    {"label": "Raw", "value": "raw"},
                                                    {
                                                        "label": "Sequence",
                                                        "value": "sequence",
                                                    },
                                                ],
                                                value="collapsed",
                                                inline=True,
                                            ),
                                            html.Label(
                                                "Minimal repeat length [in units]:"
                                            ),
                                            dcc.RangeSlider(
                                                id="repeat-len-slider",
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                dcc.Graph(id="kmer-composition"),
                                style={
                                    "display": "flex",
                                    "justifyContent": "center",
                                    "alignItems": "center",
                                },
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
                                    )
                                ]
                            ),
                            dash_table.DataTable(
                                id="details-table",
                                style_cell={
                                    "fontSize": 14,
                                    "font-family": "sans-serif",
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
                                    "fontSize": 18,
                                },
                                # style_table={"overflowX": "auto"},
                                tooltip_duration=None,
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
                        label="Your data",
                        children=[
                            dcc.Store(id="stored-df"),
                            dcc.Upload(
                                id="upload-data",
                                children=html.Div(
                                    [
                                        "Drag and drop or ",
                                        html.A("click to upload STRdust VCF.gz files"),
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
                    dcc.Tab(
                        label="About",
                        children=[
                            html.H1("About"),
                            html.Div(
                                [
                                    html.P(
                                        [
                                            "This web app was developed by Wouter De Coster. ",
                                            "The source code is available on ",
                                            html.A(
                                                "GitHub",
                                                href="https://github.com/wdecoster/pathSTR-1000G",
                                                target="_blank",
                                            ),
                                            ". Feedback is welcome in the form of an ",
                                            html.A(
                                                "issue on GitHub",
                                                href="https://github.com/wdecoster/pathSTR-1000G/issues",
                                                target="_blank",
                                            ),
                                            " or by sending me an email.",
                                        ]
                                    ),
                                    html.P(
                                        [
                                            "The STR genotypes have been obtained using ",
                                            html.A(
                                                "STRdust",
                                                href="https://github.com/wdecoster/STRdust",
                                                target="_blank",
                                            ),
                                            " from samples of the 1000 Genomes project, sequenced on the Oxford Nanopore Technologies PromethION. "
                                            "Please let me know if you know a suitable dataset of control individuals to add, for which long read alignments are available online. ",
                                            "If this resource is useful to you, please cite our ",
                                            html.A(
                                                "publication",
                                                href="https://github.com/wdecoster/pathSTR-1000G",
                                                target="_blank",
                                            ),
                                            ", as well as the references to the underlying datasets: ",
                                            html.A(
                                                "Noyvert et al. 2023",
                                                href="https://www.medrxiv.org/content/10.1101/2023.12.20.23300308v1",
                                                target="_blank",
                                            ),
                                            " and ",
                                            html.A(
                                                "Gustafson et al. 2024",
                                                href="https://github.com/wdecoster/STRdust",
                                                target="_blank",
                                            ),
                                            ".",
                                        ]
                                    ),
                                    html.P(
                                        [
                                            "You can upload your own STRdust VCF.gz file(s) to show alongside the 1000 Genomes data for comparison, but this is currently limited to 100kb files, please let me know if more would be required. Note that as the app cannot figure out the sex of your individuals, males will show two alleles on haploid chromosomes. This is corrected for samples in the 1000 Genomes cohort. "
                                            "The repeat coordinates and motifs used in this app are obtained from ",
                                            html.A(
                                                "STRchive",
                                                href="https://strhive.com/",
                                                target="_blank",
                                            ),
                                            ". Other dependencies are Python and the Dash/plotly, pandas, hdf5 and cyvcf2 modules for the web app, and snakemake to orchestrate the variant calling.",
                                        ],
                                        style={"textAlign": "justify"},
                                    ),
                                ],
                                style={"width": "80%", "margin": "auto"},
                            ),
                            html.H1("Downloads"),
                            html.Div(
                                [
                                    html.Div(
                                        children=[
                                            html.Button(
                                                "Download Data as TSV",
                                                id="btn",
                                                style={
                                                    "backgroundColor": "#86888a",  # background color
                                                    "border": "none",  # remove border
                                                    "color": "white",  # text color
                                                    "padding": "15px 32px",  # padding
                                                    "textAlign": "center",  # center text
                                                    "textDecoration": "none",  # remove underline
                                                    "fontSize": "16px",  # font size
                                                    "margin": "4px 2px",  # margin
                                                    "cursor": "pointer",  # cursor style
                                                    "borderRadius": "12px",  # rounded corners
                                                    "marginRight": "10px",
                                                },
                                            ),
                                            html.Button(
                                                "Download STRdust VCFs",
                                                id="download-zip-button",
                                                style={
                                                    "backgroundColor": "#86888a",  # background color
                                                    "border": "none",  # remove border
                                                    "color": "white",  # text color
                                                    "padding": "15px 32px",  # padding
                                                    "textAlign": "center",  # center text
                                                    "textDecoration": "none",  # remove underline
                                                    "fontSize": "16px",  # font size
                                                    "margin": "4px 2px",  # margin
                                                    "cursor": "pointer",  # cursor style
                                                    "borderRadius": "12px",  # rounded corners
                                                },
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
        return dcc.send_data_frame(df.to_csv, "pathSTR-1000G.tsv")

    @app.callback(
        [Output("stored-df", "data"), Output("upload-status", "children")],
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
        State("upload-data", "last_modified"),
    )
    def store_uploaded_data(list_of_contents, list_of_filenames, list_of_dates):
        if list_of_contents is not None:
            dfs = [
                parse.parse_uploaded_vcf(content, filename, repeats)
                for (content, filename, _) in zip(
                    list_of_contents, list_of_filenames, list_of_dates
                )
            ]
            for df, filename in zip(dfs, list_of_filenames):
                if df is None:
                    print(
                        f"Error parsing {filename}, please make sure it is a valid STRdust VCF.gz file"
                    )
                    # Create a dummy dataframe to communicate error
                    return (
                        pd.DataFrame(()).to_dict("records")
                    ), f"Error processing file {filename}"
            uploaded_df = pd.concat(dfs, ignore_index=True)
            plural = "s" if len(dfs) > 1 else ""
            return uploaded_df.to_dict("records"), f"Uploaded {len(dfs)} file{plural}"
        else:
            return pd.DataFrame(()).to_dict("records"), ""

    # @app.callback(
    #     Output("dropdown-gene-composition", "value"),
    #     Input("dropdown-gene-length", "value"),
    #     State("dropdown-gene-composition", "value"),
    # )
    # def update_gene_composition_dropdown(value_other, value_current):
    #     if value_current == value_other:
    #         return dash.no_update
    #     else:
    #         return value_other

    # @app.callback(
    #     Output("dropdown-gene-length", "value"),
    #     Input("dropdown-gene-composition", "value"),
    #     State("dropdown-gene-length", "value"),
    # )
    # def update_gene_length_dropdown(value_other, value_current):
    #     if value_current == value_other:
    #         return dash.no_update
    #     else:
    #         return value_other

    @app.callback(
        Output("violin-plot", "figure"),
        [
            Input("dropdown-gene-length", "value"),
            Input("stored-df", "data"),
            Input("violin_options", "value"),
        ],
    )
    def update_violin(selected_gene, stored_df, violin_options):
        if stored_df is None:
            filtered_df = df[df["gene"] == selected_gene]
        else:
            stored_df = pd.DataFrame(stored_df)
            combined_df = pd.concat([df, stored_df], ignore_index=True)
            filtered_df = combined_df[combined_df["gene"] == selected_gene]
        return plot.violin_plot(
            filtered_df,
            repeats=repeats,  # show optionally the pathogenic length
            selected_gene=selected_gene,
            violin_options=violin_options,
        )

    @app.callback(
        Output("length-scatter", "figure"),
        [
            Input("dropdown-gene-length", "value"),
            Input("stored-df", "data"),
            Input("violin_options", "value"),
        ],
    )
    def update_length_scatter(selected_gene, stored_df, violin_options):
        if stored_df is None:
            filtered_df = df[df["gene"] == selected_gene]
        else:
            stored_df = pd.DataFrame(stored_df)
            combined_df = pd.concat([df, stored_df], ignore_index=True)
            filtered_df = combined_df[combined_df["gene"] == selected_gene]
        return plot.length_scatter(
            filtered_df,
            path_length=repeats.pathogenic_min_length(selected_gene),  # optional
            violin_options=violin_options,
        )

    # change the settings for the kmer dropdown based on the mode
    @app.callback(
        [
            Output("dropdown-kmer-options", "multi"),
            Output("dropdown-kmer-options", "options"),
            Output("dropdown-kmer-options", "value"),
            Output("dropdown-kmer-options", "placeholder"),
        ],
        [Input("kmer_mode", "value"), Input("dropdown-gene-composition", "value")],
    )
    def disable_kmer_motifs_dropdown(mode, selected_gene):
        if mode == "raw":
            multi = True
            options = [i for i in kmers[selected_gene].columns if i != "length"]
            placeholder = "Select kmers to sort on"
            value = ""
        elif mode == "sequence":
            multi = False
            options = ["left-to-right", "right-to-left"]
            placeholder = ""
            value = "left-to-right"
        else:
            multi = False
            options = ["hide singletons", "show everything"]
            placeholder = ""
            value = "show everything"
        return multi, options, value, placeholder

    @app.callback(
        Output("kmer-composition", "figure"),
        [
            Input("dropdown-gene-composition", "value"),
            Input("dropdown-kmer-options", "value"),
            Input("repeat-len-slider", "value"),
            Input("kmer_mode", "value"),
            Input("stored-df", "data"),
        ],
    )
    def update_kmer_composition(
        selected_gene, kmer_options, length_range, kmer_mode, stored_df
    ):
        if len(stored_df) == 0:
            kmer_df = kmers[selected_gene]
        else:
            # I am not too sure if this part with user data is going to work - especially if the axis to concatenate on is correct
            stored_df = pd.DataFrame(stored_df)
            kmer_df = pd.concat(
                [kmers[selected_gene], parse_kmers(stored_df, repeats, selected_gene)]
            ).fillna(0.0)
        filtered_df = df[df["gene"] == selected_gene]
        return plot.kmer_plot(
            kmer_df,
            repeat_df=filtered_df,
            mode=kmer_mode,
            length_range=length_range,
            kmer_options=kmer_options,
        )

    @app.callback(
        Output("repeat-len-slider", "min"),
        Output("repeat-len-slider", "max"),
        Input("dropdown-gene-composition", "value"),
    )
    def update_slider_range(selected_gene):
        """Update the slider based on the minimal and maximal length of the repeat"""
        return floor(kmers[selected_gene]["length"].min()), ceil(
            kmers[selected_gene]["length"].max()
        )

    @app.callback(
        Output("strip-plot-log-store", "data"),
        Input("stored-df", "data"),
    )
    def update_stripplot(stored_df):
        if stored_df is None:
            return plot.create_strip_plot(df, log=True)
        else:
            strip_df = pd.concat([df, pd.DataFrame(stored_df)], ignore_index=True)
            return plot.create_strip_plot(strip_df, log=True)

    @app.callback(
        Output("strip-plot-log", "figure"),
        Input("strip-plot-log-store", "data"),
    )
    def update_strip_plot_from_store(strip_log_data):
        return strip_log_data

    @app.callback(
        [
            Output("details-table", "data"),
            Output("details-table", "columns"),
            Output("details-table", "tooltip_data"),
        ],
        [
            Input("dropdown-details-individual", "value"),
            Input("dropdown-details-gene", "value"),
        ],
    )
    def update_details_table(individuals, gene):
        if not individuals:
            return (
                dash.no_update,
                dash.no_update,
                dash.no_update,
            )
        if isinstance(individuals, str):
            individuals = [individuals]

        detail_df = (
            df[(df["gene"] == gene) & (df["sample"].isin(individuals))]
            .reset_index(names="sample.1")
            .drop(columns=["Group", "allele", "gene", "chrom"])
            .round(1)
            .groupby("sample")
            .transform(
                lambda x: (
                    ",".join([str(i) for i in set(x)])
                    if x.name in ["sample.1", "Sex", "Superpopulation"]
                    else ",".join([str(i) for i in list(x)])
                )
            )
            .drop_duplicates()
            .transpose()
        )
        detail_df = (
            detail_df.rename(columns=detail_df.loc["sample.1"])
            .drop("sample.1")
            .reset_index(names="")
        )
        columns = [({"name": c, "id": c}) for c in detail_df.columns]
        tooltip_data = [
            {
                column: {"value": str(value), "type": "markdown"}
                for column, value in row.items()
            }
            for row in detail_df.to_dict("records")
        ]
        return detail_df.to_dict("records"), columns, tooltip_data

    @app.callback(
        Output("igv-output", "children"),
        Input("igv-button", "n_clicks"),
        [
            State("dropdown-details-individual", "value"),
            State("dropdown-details-gene", "value"),
        ],
    )
    def return_igv(n_clicks, individuals, gene):
        if n_clicks is not None:
            if not individuals:
                return html.Div()
            if isinstance(individuals, str):
                individuals = [individuals]
            locus = repeats.gene_to_coords(gene)
            # Would be slightly more elegant to be able to get those values from the repeats object
            # but this will do for now
            # need the values separately to slightly extend the locus for better visualization
            chrom = locus.split(":")[0]
            start, end = [int(i) for i in locus.split(":")[1].split("-")]
            return html.Div(
                [
                    dashbio.Igv(
                        id="igv",
                        genome="hg38_1kg",
                        # reference={
                        #     "id": "hg38",
                        #     "name": "Human (GRCh38/hg38)",
                        #     "fastaURL": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/1KG_ONT_VIENNA_hg38.fa.gz",
                        #     "tracks": [
                        #         {
                        #             "name": "Refseq Genes",
                        #             "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz",
                        #             "indexed": False,
                        #         }
                        #     ],
                        # },
                        locus=f"{chrom}:{start-25}-{end+25}",
                        tracks=(
                            [
                                make_igv_alignment_track(individual)
                                for individual in individuals
                            ]
                        ),
                    )
                ]
            )

    def make_igv_alignment_track(individual):
        return {
            "name": individual,
            "type": "alignment",
            "url": f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/{individual}.hg38.cram",
            "indexURL": f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/{individual}.hg38.cram.crai",
            "height": 300,
            "format": "cram",
            "showSoftClips": True,
            "showInsertionText": True,
            "showDeletionText": True,
            "showSoftClips": True,
            "colorBy": "strand",
        }

    @app.callback(
        Output("download-zip", "data"),
        [Input("download-zip-button", "n_clicks")],
        prevent_initial_call=True,
    )
    def download_zip(n_clicks):
        return dcc.send_file("pathSTR_STRdust_good_samples.zip")

    @app.callback(
        Output("user-data-table", "data"),
        Output("user-data-table", "columns"),
        Input("stored-df", "data"),
    )
    def update_table(data):
        if data is None or len(data) == 0:
            return (
                dash.no_update,
                dash.no_update,
            )  # Don't update the table if no data is uploaded
        user_df = (
            pd.DataFrame(data)
            .drop(columns=["chrom", "sequence", "Superpopulation", "Group", "Sex"])
            .set_index("gene", drop=False)
            .join(stat)
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
        # user_df.to_csv("pathSTR-1000G-user-data.tsv", sep="\t")
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
    app.run(debug=args.debug)


def get_args():
    parser = ArgumentParser(description="Get STRs from VCF")
    parser.add_argument(
        "--vcf",
        nargs="+",
        help="Input VCFs",
    )
    parser.add_argument("--sample_info", help="Sample info file")
    parser.add_argument("--db", help="Input is in one pathSTR_db file")
    parser.add_argument("--save_db", help="Save the parsed data to a pathSTR_db file")
    parser.add_argument(
        "--store_only",
        help="Only store the parsed data and do not start the Dash app",
        action="store_true",
    )
    parser.add_argument(
        "--debug",
        help="Run the app in debug mode",
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
