# This script is used to serve a Dash app with the information from STR calls of the 1000 Genomes project,
# focusing on the pathogenic STRs as genotyped by STRdust.
# The Dash app has one panel focusing on the STR lengths, and one panel focusing on the STR sequence composition.

import base64
import gzip
import pandas as pd
import plotly.express as px
import dash
from dash import Dash, html, dcc, ctx
from dash.dependencies import Input, Output, State
from cyvcf2 import VCF
from argparse import ArgumentParser
import os
from itertools import chain
import zipfile
from flask import send_file


def main():
    args = get_args()
    # read in the BED file with the STRs
    coords_to_names, names_to_motif_length = get_repeat_info(args.bed)
    # read in the VCFs
    lengths = [get_lengths_from_vcf(vcf, coords_to_names) for vcf in args.vcf]
    # make a dataframe
    df = pd.DataFrame(flatten(lengths), columns=["gene", "sample", "length"])
    # for every repeat in the dataframe, divide the length by the motif length
    df["length"] = df.apply(
        lambda x: x["length"] / names_to_motif_length[x["gene"]], axis=1
    )
    df["Group"] = "1000 Genomes"
    # df.to_csv("pathSTR-dash.tsv", index=False, sep="\t")

    # Create Dash app
    app = Dash(__name__)

    # Define app layout
    gene_options = [
        {"label": gene, "value": gene} for gene in sorted(df["gene"].unique().tolist())
    ]
    app.layout = html.Div(
        [
            html.H1("Pathogenic Repeats from 1000 Genomes Project ONT resequencing"),
            dcc.Tabs(
                [
                    dcc.Tab(
                        label="Overview",
                        children=[
                            html.H1("Overview"),
                            dcc.Loading(
                                id="loading-strip",
                                type="circle",
                                children=[
                                    html.Div(dcc.Graph(id="strip-plot")),
                                    html.Div(dcc.Graph(id="strip-plot-log")),
                                ],
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
                                            "This web app was developed by Wouter De Coster. "
                                            "The STR genotypes have been obtained using ",
                                            html.A(
                                                "STRdust",
                                                href="https://github.com/wdecoster/STRdust",
                                                target="_blank",
                                            ),
                                            " from samples of the 1000 Genomes project, sequenced on ONT. "
                                            "In the length plot, each dot represents a repeat length in total repeat motifs, including those in the reference genome. "
                                            "You can upload your own STRdust VCF file(s) to show alongside the 1000 Genomes data for comparison. "
                                            "The source code is available on ",
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
                get_lengths_from_uploaded_vcf(
                    content, filename, coords_to_names, names_to_motif_length
                )
                for (content, filename, _) in zip(
                    list_of_contents, list_of_filenames, list_of_dates
                )
            ]
            uploaded_df = pd.concat(dfs, ignore_index=True)
            return uploaded_df.to_dict("records")

    @app.callback(
        Output("violin-plot", "figure"),
        [Input("dropdown-gene", "value"), Input("stored-df", "data")],
    )
    def update_violin(selected_gene, stored_df):
        if stored_df is None:
            filtered_df = df[df["gene"] == selected_gene]
        else:
            stored_df = pd.DataFrame(stored_df)
            combined_df = pd.concat([df, stored_df], ignore_index=True)
            filtered_df = combined_df[combined_df["gene"] == selected_gene]
        return create_violin_plot(filtered_df)

    @app.callback(
        [Output("strip-plot", "figure"), Output("strip-plot-log", "figure")],
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

    # Run the app
    app.run_server(debug=True)


def create_violin_plot(filtered_df):
    fig = px.violin(
        filtered_df,
        x="gene",
        y="length",
        color="Group",
        points="all",
        hover_data=["sample"],
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(xaxis_title="", yaxis_title="Repeat length [units]")
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
        fig.update_layout(xaxis_title="", yaxis_title="log-Repeat length [units]")
    else:
        fig.update_layout(xaxis_title="", yaxis_title="Repeat length [units]")
    fig.update_traces(marker=dict(size=2))
    if strip_df["Group"].nunique() > 1:
        fig.update_layout(legend_title_text="Group")
    else:
        fig.update_layout(showlegend=False)
    return fig


def get_repeat_info(bed):
    """
    This function parses a bed file as obtained from STRchive, which is intended for TRGT
    """
    coords_to_gene = {}
    gene_to_motif_length = {}
    bed = pd.read_csv(bed, sep="\t", header=None)
    bed.columns = ["chrom", "start", "end", "info"]
    bed["name"] = bed["info"].apply(
        lambda x: [i.split("_")[1] for i in x.split(";") if i.startswith("ID=")][
            0
        ].replace("ID=", "")
    )
    bed["motif"] = bed["info"].apply(
        lambda x: [i.split(",")[0] for i in x.split(";") if i.startswith("MOTIFS=")][
            0
        ].replace("MOTIF=", "")
    )
    # in case there are duplicates in the name column, use the original ID from the info field
    dups = bed.duplicated(subset="name", keep=False)
    bed.loc[dups, "name"] = bed.loc[dups, "info"].apply(
        lambda x: [change_repeat_name(i) for i in x.split(";") if i.startswith("ID=")][
            0
        ].replace("ID=", "")
    )
    for chrom, start, end, name, motif in bed[
        ["chrom", "start", "end", "name", "motif"]
    ].values:
        coords_to_gene[(chrom, str(start), str(end))] = name
        gene_to_motif_length[name] = len(motif)
    return coords_to_gene, gene_to_motif_length


def change_repeat_name(name):
    return f"{name.split('_')[1]}_{name.split('_')[0]}"


def flatten(it):
    return chain.from_iterable(it)


def get_lengths_from_vcf(vcf, coords_to_names):
    calls = []
    name = os.path.basename(vcf).replace(".vcf.gz", "")
    for v in VCF(vcf):
        gene = coords_to_names[(v.CHROM, str(v.POS), str(v.end))]
        lengths = v.INFO.get("FRB")
        calls.append((gene, name, lengths[0]))
        calls.append((gene, name, lengths[1]))
    return calls


def get_lengths_from_uploaded_vcf(
    contents, filename, coords_to_names, names_to_motif_length
):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    # write the decoded file to a temporary file
    # make sure the file ends with .gz
    # and read it back in using cyvcf2
    # now perhaps this could be done more efficiently, but it works
    if filename.endswith(".gz"):
        base = os.path.basename(filename)
        tempfile = os.path.join("/tmp", base)
        with open(tempfile, "wb") as f:
            f.write(decoded)
    else:
        base = os.path.basename(filename) + ".gz"
        tempfile = os.path.join("/tmp", base)
        with gzip.open(tempfile, "wb") as f:
            f.write(decoded)
    calls = get_lengths_from_vcf(tempfile, coords_to_names)
    df = pd.DataFrame(calls, columns=["gene", "sample", "length"])
    df["length"] = df.apply(
        lambda x: x["length"] / names_to_motif_length[x["gene"]], axis=1
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
    parser.add_argument("--bed", help="STRchive BED file with STRs")
    return parser.parse_args()


if __name__ == "__main__":
    main()
