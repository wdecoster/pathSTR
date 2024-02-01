# This script is used to serve a Dash app with the information from STR calls of the 1000 Genomes project,
# focusing on the pathogenic STRs as genotyped by STRdust.


import pandas as pd
import dash
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
from argparse import ArgumentParser
import zipfile
from flask import send_file
import parse_input as parse
from repeats import Repeats
import plot


def main():
    args = get_args()
    # read in the BED file with the STRs
    repeats = Repeats(args.bed)
    df = parse.parse_input(args.vcf, args.sample_info, args.feather, repeats)

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
                                children=[
                                    dcc.Dropdown(
                                        id="dropdown-gene",
                                        options=gene_options,
                                        value=gene_options[0]["value"],
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
                                        ],
                                        value=[],
                                        inline=True,
                                    ),
                                ]
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
                                            "You can upload your own STRdust VCF.gz file(s) to show alongside the 1000 Genomes data for comparison, but this is currently limited to 100kb files, please let me know if more would be required. The source code is available on ",
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
                                            ". Please let me know if you know a suitable dataset to add. Feedback is welcome in the form of an ",
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
        [Output("stored-df", "data"), Output("upload-status", "children")],
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
        State("upload-data", "last_modified"),
    )
    def store_uploaded_data(list_of_contents, list_of_filenames, list_of_dates):
        if list_of_contents is not None:
            dfs = [
                parse.get_lengths_from_uploaded_vcf(content, filename, repeats)
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
            return uploaded_df.to_dict("records"), f"Uploaded {len(dfs)} files"
        else:
            return dash.no_update, dash.no_update

    @app.callback(
        [Output("violin-plot", "figure"), Output("violin-plot-log", "figure")],
        [
            Input("dropdown-gene", "value"),
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
            log=False,
            violin_options=violin_options,
        ), plot.violin_plot(
            filtered_df,
            log=True,
            violin_options=violin_options,
        )

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
        strip = plot.create_strip_plot(strip_df)
        strip_log = plot.create_strip_plot(strip_df, log=True)
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


def get_args():
    parser = ArgumentParser(description="Get STRs from VCF")
    parser.add_argument(
        "--vcf",
        nargs="+",
        help="Input VCFs",
    )
    parser.add_argument("--bed", help="STRchive BED file with STRs")
    parser.add_argument("--sample_info", help="Sample info file")
    parser.add_argument("--feather", help="Input is in one feather file")
    return parser.parse_args()


if __name__ == "__main__":
    main()
