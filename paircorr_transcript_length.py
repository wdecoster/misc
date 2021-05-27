from argparse import ArgumentParser
from pysam import AlignmentFile
from scipy.stats import pearsonr
import plotly.express as px
import pandas as pd
from random import random


def main():
    args = get_args()
    df_plots = []
    names = args.names if args.names else [bamf.replace('.bam', '') for bamf in args.bam]
    for bamf, name in zip(args.bam, names):
        df = get_reads(bamf, sampling=args.sample)
        df_plots.append(df.sample(n=min(args.maxdots, len(df))).assign(sample=name))

        pearson_correlations(df, name)

    scatter(pd.concat(df_plots))


def get_reads(bamf, sampling=0.05):
    bam = AlignmentFile(bamf)
    read_lengths = []
    transcript_lengths = []
    for read in bam.fetch():
        if read.mapq > 0 and random() < sampling:
            read_lengths.append(read.reference_length)
            transcript_lengths.append(bam.get_reference_length(read.reference_name))
    return pd.DataFrame(list(zip(read_lengths, transcript_lengths)),
                        columns=['reads', 'transcripts'])


def pearson_correlations(df, name):
    corr, pval = pearsonr(x=df['reads'], y=df['transcripts'])
    print(f"Pearson correlation coefficient overall for {name}: {corr}, p-value {pval}")

    df_shorter = df.loc[df["transcripts"] < 2000]
    corr, pval = pearsonr(x=df_shorter['reads'], y=df_shorter['transcripts'])
    print("Pearson correlation coefficient for transcripts <2kb"
          f" for {name}: {corr}, p-value {pval}")

    df_shorter = df_shorter.loc[df_shorter["transcripts"] < 1000]
    corr, pval = pearsonr(x=df_shorter['reads'], y=df_shorter['transcripts'])
    print("Pearson correlation coefficient for transcripts <1kb"
          f" for {name}: {corr}, p-value {pval}")


def scatter(all_df):
    fig = px.scatter(all_df, x="reads", y="transcripts", color="sample", trendline="ols",
                     marginal_x='violin', marginal_y='violin')
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(title="Read length vs MANE",
                      xaxis_title="Read lengths",
                      yaxis_title="MANE transcript lengths")
    fig.update_yaxes(range=[0, 1e4])

    fig_logged = px.scatter(all_df, x="reads", y="transcripts", color="sample",
                            trendline="ols", log_x=True, log_y=True,
                            marginal_x='violin', marginal_y='violin')
    fig_logged.update_traces(marker=dict(size=3))
    fig_logged.update_layout(title="Read length vs MANE",
                             xaxis_title="log-transformed read lengths",
                             yaxis_title="log-transformed MANE transcript lengths")
    with open("read-length-plot.html", 'w') as output:
        output.write(fig.to_html(include_plotlyjs='cdn'))
        output.write(fig_logged.to_html(include_plotlyjs='cdn'))


def get_args():
    parser = ArgumentParser(description='pairwise correlation of length of transcript and reads')
    parser.add_argument('bam', nargs='+', help='bam file of reads aligned to transcriptome')
    parser.add_argument('-n', '--names', help='names for bam files', nargs='*')
    parser.add_argument('-s', '--sample',
                        help='fraction to subsample reads on during collection',
                        default=0.05)
    parser.add_argument('--maxdots',
                        help='maximum number of dots to plot for each dataset',
                        default=10000)
    return parser.parse_args()


if __name__ == '__main__':
    main()
