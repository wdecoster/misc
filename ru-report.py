from argparse import ArgumentParser
import tempfile
import pysam
import numpy as np


def main():
    args = get_args()
    off_target_bed = complement_bed(args.bam, args.bed, args.workdir)
    coverage_table = coverage(args.bam, args.bed, off_target_bed, args.workdir)
    length_plots = overlay_histogram(args.bam, args.bed, off_target_bed)
    make_report(content=[coverage_table, length_plots], outputf=args.output)


def get_args():
    parser = ArgumentParser(description="Make a report from a read until experiment")
    parser.add_argument("--bam", help="Sorted and indexed bam file.", required=True)
    parser.add_argument("--bed", help="Targets in bed-4 format.", required=True)
    parser.add_argument("--workdir", help="Directory for intermediate file",
                        default=tempfile.mkdtemp())
    parser.add_argument("--output", help="Output filename (.html) with existing path.",
                        default="ru-report.html")
    return parser.parse_args()


def coverage(bam, bed, off_target_bed, workdir):
    coverage = mosdepth(bam, bed, workdir, prefix="ontarget")
    return coverage.append(off_target_coverage(bam, off_target_bed, workdir),
                           ignore_index=True) \
        .to_html(index=False,
                 float_format=lambda x: f'{x:.2f}',
                 columns=['name', 'chrom', 'begin', 'end', 'coverage'],
                 na_rep='-')


def mosdepth(bam, bed, dir, prefix):
    """
    Call mosdepth with a certain bed file (which can have 3 or 4 columns)
    And return the .regions.bed.gz file as a pandas DataFrame
    """
    import pandas as pd
    from subprocess import call
    import shlex
    if len(open(bed).readline().split('\t')) == 4:
        names = ['chrom', 'begin', 'end', 'name', 'coverage']
    else:
        names = ['chrom', 'begin', 'end', 'coverage']
    call(shlex.split(f"mosdepth --fast-mode --no-per-base --by {bed} -t4 {dir}/{prefix} {bam}"))
    return pd.read_csv(f"{dir}/{prefix}.regions.bed.gz", sep="\t", header=None,
                       names=names)


def off_target_coverage(bam, bed, dir):
    """
    Calculate (using mosdepth and the off_target_bed) the mean coverage in the off target regions,
    with coverage weighted by length of the intervals.
    Return as dictionary to be inserted in the coverage DataFrame
    """
    df = mosdepth(bam, bed, dir, prefix="offtarget")
    df["size"] = df["end"] - df["begin"]
    return {'name': 'offtarget',
            'coverage': (df["coverage"] * df["size"]).sum() / df["size"].sum()}


def read_lengths(bamf, bed):
    """
    Extract the aligned read lengths in the bam file based on a bed file
    """
    bam = pysam.AlignmentFile(bamf)
    lengths = []
    for line in open(bed):
        if line:
            chrom, start, end = line.strip().split('\t')[0:3]
            for read in bam.fetch(reference=chrom, start=int(start), end=int(end)):
                lengths.append(read.query_alignment_length)
    return np.array(lengths)


def complement_bed(bamf, bed, dir):
    """
    Create a complement of the bed file with targets
    To be used for off target reads
    """
    from pybedtools import BedTool
    bam = pysam.AlignmentFile(bamf)
    chromsizesf = f"{dir}/chromsizes.txt"
    complement = f"{dir}/off_target_bed"
    with open(chromsizesf, 'w') as chromsizes:
        for contig, length in zip(bam.references, bam.lengths):
            chromsizes.write(f"{contig}\t{length}\n")
    BedTool(bed).sort(g=chromsizesf).complement(g=chromsizesf).saveas(complement)
    return complement


def downsample(l, n):
    """
    Downsample a numpy array l to a new array with length n
    Or return l entirely if it's smaller than n
    """
    if l.size <= n:
        return l
    else:
        return np.random.choice(l, size=n)


def overlay_histogram(bam, bed, off_target_bed):
    """
    Create a two-column subplot with normalised density plots of read lengths
    with left the normal and right the log-transformed lengths
    """
    import plotly
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Read lengths", "Log-transformed read lengths"))
    on = downsample(read_lengths(bam, bed), 1000)
    off = downsample(read_lengths(bam, off_target_bed), 1000)

    fig.add_histogram(x=on, opacity=0.4, name="On Target",
                      histnorm="density", row=1, col=1, marker=dict(color="green"))
    fig.add_histogram(x=off, opacity=0.4, name="Off Target",
                      histnorm="density", row=1, col=1, marker=dict(color="blue"))
    fig.add_histogram(x=np.log10(on), opacity=0.4, name="On Target",
                      histnorm="density", row=1, col=2, marker=dict(color="green"))
    fig.add_histogram(x=np.log10(off), opacity=0.4, name="Off Target",
                      histnorm="density", row=1, col=2, marker=dict(color="blue"))
    fig.update_layout(barmode='overlay')
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * max(np.amax(on), np.amax(off))]
    fig["layout"]["xaxis2"].update(tickvals=np.log10(xtickvals), ticktext=xtickvals)
    fig.update_xaxes(title_text="read lengths", row=1, col=1, rangemode='nonnegative')

    return plotly.offline.plot(fig, output_type="div", show_link=False, include_plotlyjs='cdn')


def make_report(content, outputf):
    with open(outputf, 'w') as output:
        for item in content:
            output.write(item)


if __name__ == '__main__':
    main()
