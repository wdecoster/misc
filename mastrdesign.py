from cyvcf2 import VCF
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import pyranges as pr
import sys
from math import ceil


def main():
    args = get_args()
    targets = get_input(args.input, extension=args.extension)
    make_targets_file(targets.df, flanks=args.flanking, prefix=args.prefix)
    get_fasta(targets, args.reference, flanks=args.flanking, prefix=args.prefix)
    if args.masked_reference:
        get_fasta(targets, args.masked_reference, flanks=args.flanking,
                  prefix=args.prefix, suffix=".masked")


def get_input(input_files, extension=15):

    def is_vcf(filename):
        return True if filename.endswith(('.vcf', '.vcf.gz')) else False

    def vcf2df(vcf):
        return pd.DataFrame.from_records([(v.CHROM, v.start, v.end, get_name(v)) for v in VCF(vcf)],
                                         columns=['Chromosome', 'Start', 'End', 'Name'])

    def get_name(variant):
        """tailored to snpeff and vep, to be adjusted for others"""
        if variant.INFO.get('ANN'):  # SnpEff
            gene = variant.INFO.get('ANN').split('|')[3]
        elif variant.INFO.get('CSQ'):  # VEP
            gene = variant.INFO.get('CSQ').split('|')[3]
        else:
            sys.exit("ERROR: Did not recognize annotation!")
        return gene if gene else 'no_gene'

    def is_gtf(filename):
        return True if filename.endswith(('.gtf', '.gtf.gz')) else False

    def gtf2df(gtf):
        df = pr.read_gtf(gtf, as_df=True)
        df = df.loc[df["Feature"] == 'exon', ['Chromosome', 'Start', 'End', 'gene_name']]
        df.columns = ['Chromosome', 'Start', 'End', 'Name']
        return df

    def split_too_long_targets(gr, maxlen=200):
        '''
        Split targets which are too long in targets of maximally maxlen size
        '''
        long_targets = gr.lengths() > maxlen
        gr_long = gr[long_targets]
        new_intervals = []
        for chrom, start, end in gr_long.df.itertuples(index=False, name=None):
            num_targets = ceil((end - start) / maxlen)
            new_target_size = ceil((end - start) / num_targets)
            for i in range(num_targets):
                new_intervals.append((chrom, start+i*new_target_size, start+(i+1)*new_target_size))
        new_intervals = pd.DataFrame(data=new_intervals, columns=['Chromosome', 'Start', 'End'])
        return pr.PyRanges(pd.concat([gr[~long_targets].df, new_intervals]))

    def merge_near_targets(gr, distance=100):
        clustered = gr.cluster(slack=distance)
        clustered_targets = clustered.df["Cluster"].duplicated(keep=False)
        new_intervals = []
        for cluster in clustered[clustered_targets].df["Cluster"].unique():
            cl = clustered[clustered.df["Cluster"] == cluster]
            new_intervals.append(
                (cl.Chromosome.values[0], cl.Start.values[0], cl.End.values[len(cl)-1]))
        new_intervals = pd.DataFrame(data=new_intervals, columns=['Chromosome', 'Start', 'End'])
        return pr.PyRanges(pd.concat([gr[~clustered_targets].df, new_intervals]))

    vcf_targets = [vcf2df(f) for f in input_files if is_vcf(f)]
    gtf_targets = [gtf2df(f) for f in input_files if is_gtf(f)]
    targets = pd.concat(vcf_targets + gtf_targets, ignore_index=True)
    gr = pr.PyRanges(targets).sort().merge()
    gr = merge_near_targets(gr)
    gr = gr.slack(extension)
    gr = split_too_long_targets(gr)
    gr = gr.nearest(pr.PyRanges(targets), suffix="_ori")
    gr.Identifier = pd.Series(np.arange(1, len(gr) + 1)).astype(str) + '_' + gr.Name.values
    return gr


def make_targets_file(df, flanks=400, prefix="design"):
    '''
    Create .fas.target file for mpcr.
    '''
    df["size"] = df["End"] - df["Start"]
    with open(f"{prefix}.fas.target", "w") as targets:
        for identifier, size in df[['Identifier', 'size']].itertuples(index=False, name=None):
            targets.write(f">{identifier}\n{flanks},{size}\n")
    return df


def get_fasta(gr, fas, flanks=400, prefix="design", suffix=""):
    gr_temp = gr.copy()
    gr_temp.Start -= flanks
    gr_temp.End += flanks
    gr_temp.seq = pr.get_fasta(gr_temp, fas)
    with open(f"{prefix}.fas{suffix}", "w") as fas:
        for identifier, seq in gr_temp.df[["Identifier", "seq"]].itertuples(index=False, name=None):
            fas.write(f">{identifier}\n{seq}\n")
    gr_temp.to_bed(path=f"{prefix}.bed")


def get_args():
    parser = ArgumentParser(description="Make targets for MASTR design.")
    parser.add_argument("input", nargs='+', help="vcf or gtf file to use for targets")
    parser.add_argument("-f", "--flanking",
                        help="length of flanking sequence for primers",
                        default=400,
                        type=int)
    parser.add_argument("-e", "--extension",
                        help="length of target extension",
                        default=15,
                        type=int)
    parser.add_argument("-p", "--prefix",
                        help="prefix of output files",
                        default="design",
                        type=str)
    parser.add_argument("-r", "--reference",
                        help="path to reference fasta file",
                        required=True,
                        type=str)
    parser.add_argument("-m", "--masked_reference",
                        help="path to masked reference fasta file",
                        type=str)
    return parser.parse_args()


if __name__ == '__main__':
    main()
