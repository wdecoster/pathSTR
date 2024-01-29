# This script is used to serve a Dash app with the information from STR calls of the 1000 Genomes project,
# focusing on the pathogenic STRs as genotyped by STRdust.
# The Dash app has one panel focusing on the STR lengths, and one panel focusing on the STR sequence composition.

import base64
import gzip
import pandas as pd
import plotly.express as px
import dash
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
from cyvcf2 import VCF
from argparse import ArgumentParser
import os
from itertools import chain
import zipfile
from flask import send_file


class Repeats(object):
    def __init__(self, bed):
        self.df = self.get_repeat_info(bed)

    def get_repeat_info(self, bed):
        """
        This function parses a bed file as obtained from STRchive, which is intended for TRGT,
        but also works for STRdust.
        """
        bed = pd.read_csv(
            bed, sep="\t", header=None, names=["chrom", "start", "end", "info"]
        )
        # the TRGT bed file has as ID the <disease>_<gene> and we only care about the gene
        bed["name"] = bed["info"].apply(
            lambda x: [i.split("_")[1] for i in x.split(";") if i.startswith("ID=")][
                0
            ].replace("ID=", "")
        )
        # in case there are duplicates in the name column, use the original ID from the info field
        dups = bed.duplicated(subset="name", keep=False)
        bed.loc[dups, "name"] = bed.loc[dups, "info"].apply(
            lambda x: [
                Repeats.fix_name(i) for i in x.split(";") if i.startswith("ID=")
            ][0].replace("ID=", "")
        )
        bed["motifs"] = bed["info"].apply(
            lambda x: [
                i.replace("MOTIF=", "").split(",")
                for i in x.split(";")
                if i.startswith("MOTIFS=")
            ][0]
        )
        bed["motif_length"] = bed["motifs"].apply(lambda x: len(x[0]))
        bed["id"] = (
            bed["chrom"] + ":" + bed["start"].astype(str) + "-" + bed["end"].astype(str)
        )
        return bed.drop(columns=["info", "chrom", "start", "end"]).set_index("id")

    @staticmethod
    def fix_name(name):
        return f"{name.split('_')[1]}_{name.split('_')[0]}"

    def motif_length(self, gene):
        return self.df.loc[self.df["name"] == gene, "motif_length"].values[0]

    def gene(self, id):
        return self.df.loc[id, "name"]


def main():
    args = get_args()
    if args.vcf and args.bed:
        # read in the BED file with the STRs
        repeats = Repeats(args.bed)
        # read in the VCFs
        lengths = [get_lengths_from_vcf(vcf, repeats) for vcf in args.vcf]
        # make a dataframe
        df = pd.DataFrame(
            flatten(lengths), columns=["gene", "sample", "length", "ref_diff"]
        )
        # for every repeat in the dataframe, divide the length by the motif length
        df["length"] = df.apply(
            lambda x: x["length"] / repeats.motif_length(x["gene"]), axis=1
        )
        df["ref_diff"] = df.apply(
            lambda x: x["ref_diff"] / repeats.motif_length(x["gene"]), axis=1
        )
        df["Group"] = "1000 Genomes"
        df.to_feather("pathSTR-1000G.feather")
    elif args.arrow:
        df = pd.read_feather(args.arrow)
    else:
        raise ValueError("Please provide either --vcf and --bed or --arrow")

    # Create Dash app
    app = Dash(__name__)

    # Define app layout
    gene_options = [
        {"label": gene, "value": gene} for gene in sorted(df["gene"].unique().tolist())
    ]
    app.layout = html.Div(
        [
            html.Div(
                [
                    html.H1(
                        "Pathogenic Repeats from 1000 Genomes Project nanopore resequencing"
                    ),
                ],
                style={
                    "backgroundColor": "#f9f9f9",  # Change as needed
                    "border": "1px solid #ddd",  # Change as needed
                    "textAlign": "center",  # Center the text
                },
            ),
            dcc.Tabs(
                [
                    dcc.Tab(
                        label="Overview",
                        children=[
                            html.H1("Overview", style={"bottommargin": "0px"}),
                            dcc.Store(id="strip-plot-store"),
                            dcc.Store(id="strip-plot-log-store"),
                            dcc.Loading(
                                id="loading-strip-1",
                                type="cube",
                                children=[
                                    html.Div(
                                        dcc.Graph(id="strip-plot"),
                                        style={"margin": "0px"},
                                    )
                                ],
                                style={"margin": "0px", "height": "20vh"},
                            ),
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
                                dcc.Dropdown(
                                    id="dropdown-gene",
                                    options=gene_options,
                                    value=gene_options[0]["value"],
                                )
                            ),
                            html.Div(dcc.Graph(id="violin-plot")),
                            html.Div(dcc.Graph(id="violin-plot-log")),
                        ],
                    ),
                    dcc.Tab(
                        label="Repeat Composition",
                        children=[
                            html.H1("Repeat Composition"),
                            # Add your code for the repeat composition tab here
                            # ...
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
                                        html.A("click to upload STRdust VCF files"),
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
                            ),
                            dash_table.DataTable(id="user-data-table"),
                        ],
                    ),
                    dcc.Tab(
                        label="Downloads",
                        children=[
                            html.H1("Download"),
                            html.Div(
                                [
                                    # This button triggers the download
                                    html.Button("Download Data as TSV", id="btn"),
                                    # This component handles the download
                                    dcc.Download(id="download"),
                                ]
                            ),
                            html.Div(
                                [
                                    # Another button to download all VCFs as a ZIP
                                    html.Button(
                                        "Create ZIP file for download",
                                        id="download-zip-button",
                                    ),
                                    html.A(
                                        "Download ZIP",
                                        id="download-zip-link",
                                        download="pathSTR-1000G-vcfs.zip",
                                        href="",
                                        target="_blank",
                                    ),
                                ]
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
                                            "This web app was developed by Wouter De Coster. The STR genotypes have been obtained using ",
                                            html.A(
                                                "STRdust",
                                                href="https://github.com/wdecoster/STRdust",
                                                target="_blank",
                                            ),
                                            " from samples of the 1000 Genomes project, sequenced on the Oxford Nanopore Technologies PromethION. In the length plot, each dot represents a repeat length in total repeat motifs, including those in the reference genome. "
                                            "You can upload your own STRdust VCF file(s) to show alongside the 1000 Genomes data for comparison. The source code is available on ",
                                            html.A(
                                                "GitHub",
                                                href="https://github.com/wdecoster/pathSTR-1000G",
                                                target="_blank",
                                            ),
                                            ". If this resource is useful to you, please cite our ",
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
                                            ". Feedback is welcome in the form of an ",
                                            html.A(
                                                "issue on GitHub",
                                                href="https://github.com/wdecoster/pathSTR-1000G/issues",
                                                target="_blank",
                                            ),
                                            ". The repeat coordinates and motifs used in this app are obtained from ",
                                            html.A(
                                                "STRchive",
                                                href="https://strhive.com/",
                                                target="_blank",
                                            ),
                                            ". Other dependencies are Python and the Dash and cyvcf2 modules for the web app, and snakemake to orchestrate the variant calling.",
                                        ],
                                        style={"textAlign": "justify"},
                                    ),
                                ],
                                style={"width": "80%", "margin": "auto"},
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
        Output("stored-df", "data"),
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
        State("upload-data", "last_modified"),
    )
    def store_uploaded_data(list_of_contents, list_of_filenames, list_of_dates):
        if list_of_contents is not None:
            dfs = [
                get_lengths_from_uploaded_vcf(content, filename, repeats)
                for (content, filename, _) in zip(
                    list_of_contents, list_of_filenames, list_of_dates
                )
            ]
            uploaded_df = pd.concat(dfs, ignore_index=True)
            return uploaded_df.to_dict("records")

    @app.callback(
        [Output("violin-plot", "figure"), Output("violin-plot-log", "figure")],
        [Input("dropdown-gene", "value"), Input("stored-df", "data")],
    )
    def update_violin(selected_gene, stored_df):
        if stored_df is None:
            filtered_df = df[df["gene"] == selected_gene]
        else:
            stored_df = pd.DataFrame(stored_df)
            combined_df = pd.concat([df, stored_df], ignore_index=True)
            filtered_df = combined_df[combined_df["gene"] == selected_gene]
        return violin_plot(filtered_df), violin_plot(filtered_df, log=True)

    @app.callback(
        [Output("strip-plot-store", "data"), Output("strip-plot-log-store", "data")],
        Input("stored-df", "data"),
    )
    def update_stripplot(stored_df):
        if stored_df is None:
            strip_df = df.sort_values("gene")
        else:
            stored_df = pd.DataFrame(stored_df)
            strip_df = pd.concat([df, stored_df], ignore_index=True).sort_values("gene")
        strip = create_strip_plot(strip_df)
        strip_log = create_strip_plot(strip_df, log=True)
        return strip, strip_log

    @app.callback(
        [Output("strip-plot", "figure"), Output("strip-plot-log", "figure")],
        [Input("strip-plot-store", "data"), Input("strip-plot-log-store", "data")],
    )
    def update_strip_plot_from_store(strip_data, strip_log_data):
        return strip_data, strip_log_data

    @app.callback(
        Output("download-zip-link", "href"),
        [Input("download-zip-button", "n_clicks")],
    )
    def generate_zip(n_clicks):
        if n_clicks is not None:
            with zipfile.ZipFile("pathSTR-1000G-vcfs.zip", "w") as zipf:
                for file in args.vcf:
                    zipf.write(file)
            return "/download_zip/"

    @app.server.route("/download_zip/")
    def download_zip():
        return send_file(
            "pathSTR-1000G-vcfs.zip",
            mimetype="zip",
            as_attachment=True,
            attachment_filename="pathSTR-1000G-vcfs.zip",
        )

    @app.callback(
        Output("user-data-table", "data"),
        Output("user-data-table", "columns"),
        Input("stored-df", "data"),
    )
    def update_table(data):
        if data is None:
            return (
                dash.no_update,
                dash.no_update,
            )  # Don't update the table if no data is uploaded

        user_df = pd.DataFrame(data)
        columns = [{"name": i, "id": i} for i in user_df.columns]
        return user_df.to_dict("records"), columns

    # Run the app
    app.run_server(debug=True)


def violin_plot(filtered_df, log=False):
    fig = px.violin(
        filtered_df,
        x="gene",
        y="length",
        color="Group",
        points="all",
        hover_data=["sample"],
    )
    if log:
        fig.update_layout(
            yaxis_type="log",
            xaxis_title="",
            yaxis_title="Repeat length [log(units)]",
        )
        # fig.update_layout(xaxis_range=[1, filtered_df["length"].max()])
    else:
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [units]")
    fig.update_traces(marker=dict(size=3))
    if filtered_df["Group"].nunique() > 1:
        fig.update_layout(legend_title_text="Group")
    else:
        fig.update_layout(showlegend=False)
    return fig


def create_strip_plot(strip_df, log=False):
    fig = px.strip(
        strip_df,
        x="gene",
        y="length",
        color="Group",
        stripmode="overlay",
        hover_data=["sample"],
    )
    if log:
        fig.update_layout(yaxis_type="log")
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [log(units)]")
    else:
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [units]")
    fig.update_traces(marker=dict(size=2))
    if strip_df["Group"].nunique() > 1:
        fig.update_layout(legend_title_text="Group")
    else:
        fig.update_layout(showlegend=False)
    return fig


def flatten(it):
    return chain.from_iterable(it)


def get_lengths_from_vcf(vcf, repeats):
    calls = []
    name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        gene = repeats.gene(f"{v.CHROM}:{str(v.POS)}-{str(v.end)}")
        full_lengths = v.INFO.get("FRB")
        ref_diff = v.INFO.get("RB")
        calls.append((gene, name, full_lengths[0], ref_diff[0]))
        calls.append((gene, name, full_lengths[1], ref_diff[1]))
    return calls


def get_lengths_from_uploaded_vcf(contents, filename, repeats):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    # write the decoded file to a temporary file
    # make sure the file ends with .gz and is gzipped
    # and read it back in using cyvcf2
    if filename.endswith(".gz"):
        tempfile = os.path.join("/tmp", os.path.basename(filename))
        with open(tempfile, "wb") as f:
            f.write(decoded)
    else:
        tempfile = os.path.join("/tmp", os.path.basename(filename) + ".gz")
        with gzip.open(tempfile, "wb") as f:
            f.write(decoded)
    try:
        calls = get_lengths_from_vcf(tempfile, repeats)
    except OSError:
        # this happens when the file is not a VCF, or malformed
        os.remove(tempfile)
        return None
    df = pd.DataFrame(calls, columns=["gene", "sample", "length"])
    df["length"] = df.apply(
        lambda x: x["length"] / repeats.motif_length(x["gene"]),
        axis=1,
    )
    df["Group"] = "Uploaded"
    os.remove(tempfile)
    return df


def get_args():
    parser = ArgumentParser(description="Get STRs from VCF")
    parser.add_argument(
        "--vcf",
        nargs="+",
        help="Input VCFs",
    )
    parser.add_argument("--arrow", help="Input is in one Arrow file")
    parser.add_argument("--bed", help="STRchive BED file with STRs")
    return parser.parse_args()


if __name__ == "__main__":
    main()
