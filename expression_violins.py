import plotly.graph_objects as go
from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t") \
        .drop(columns=["GeneId", "Chromosome", "Start", "End", "GeneBiotype"]) \
        .set_index("GeneName")

    df = df.drop(columns=['Length']).div((df["Length"] / 1000), axis=0)
    df = df.div(df.sum() / 1e6, axis=1)

    plots = [make_violins(df, gene) for gene in args.genes]
    html_head = """<!DOCTYPE html>
    <html>
    <head>
    <meta charset="UTF-8">
    <style>
    </style>
    <title>Violinplots expression</title>
    </head>"""
    with open(args.output, 'w') as output:
        output.write(html_head)
        for p in plots:
            if p:
                output.write(p)


def make_violins(df, gene):
    try:
        expression = df.loc[gene].to_frame()
    except KeyError:
        print("Gene {} not found, skipping.".format(gene))
        return None
    expression.columns = [gene]
    expression["group"] = expression.index.str.split('-').str[0]
    expression.sort_values(by=["group"], inplace=True)
    fig = go.Figure(data=go.Violin(y=expression[gene],
                                   x=expression["group"],
                                   pointpos=0,
                                   jitter=0.5,
                                   points="all"),
                    layout=dict(title=gene,
                                yaxis=dict(title='TPM normalized expression',
                                           rangemode='nonnegative')))
    return fig.to_html(full_html=False, include_plotlyjs='cdn')


def get_args():
    parser = ArgumentParser(description="make PCA")
    parser.add_argument("-i", "--input", required=True, help="tsv to plot violins from")
    parser.add_argument("-g", "--genes", required=True, nargs='+',
                        help="gene name(s) to make violins for")
    parser.add_argument("-o", "--output", default="violins.html", help="html file to write to")
    return parser.parse_args()


if __name__ == '__main__':
    main()
