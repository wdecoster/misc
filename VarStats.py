#!/complgen/bin/anaconda/bin/python2.7
# Script to list carriers of a certain variant in a new field
# First argument is tsv file in which the field has to be added, alternatively, reads from stdin
# Optional second argument is a sample info file
# Script makes certain assumptions about the structure of variant and sample info file and sample notation compatible with mastr data
# wdecoster

from __future__ import print_function
import sys
import os
import argparse


def getInput():
    '''Checking for input on stdin or file'''
    if args.target in ['-', 'stdin']:
        return sys.stdin
    else:
        if os.path.isfile(args.target):
            return open(args.target, 'r')
        else:
            sys.exit("File not found, is the path correct?")


def extractInfoDict(input, countmode, varcaller):
    '''
    Creating a dictionary with useful information about the input file
    Will break if specific naming pattern is not followed
    '''
    infodict = {'comments': [], countmode: {}}
    for line in input:
        if line.startswith('#'):
            infodict['comments'].append(line)
            continue
        else:
            # The first line after the comments is expected to be the header
            infodict['header'] = makeHeader(line)
            break
    else:
        sys.exit("ERROR: could not find the header line in the file!")
    # Samples are collected by looping over the header based on the countmode and the varcaller, if they're not blanco or ceph
    for index, field in enumerate(infodict['header']):
        if field.startswith(countmode + '-' + varcaller) and not field.split('-')[-1].startswith('blanco') and not field.split('-')[-1].startswith('ceph'):
            # mapping the sample index : sample name
            infodict[countmode][index] = field.split('-')[-1]
    if len(infodict[countmode].keys()) == 0:
        sys.exit("Error: could not find a single '{}-{}' field in the input file.".format(countmode, varcaller))
    return infodict


def extractSampleDict(samples):
    '''
    Extract sample information based on tab separated sample info file
    Performs sanity checks on samples in use
    If ran without a sample info file, returns a bogus dictionary mapping all samples to 'notspecified'
    '''
    # Remove the ngs id from the samples which were found in the variantfile
    simplifiedsamples = ['_'.join(sample.split('_')[:-1]) for sample in samples if not sample.lower(
    ).startswith('blanco') and not sample.lower().startswith('ceph')]
    if args.samples:
        sampleDict = {}
        with args.samples as sampleinfoFile:
            sampleinfo = [line.strip() for line in sampleinfoFile.readlines()
                          if not line == "" and not line.startswith('#')]
            if not sampleinfo[0].lower() == "sample\tcondition":
                sys.exit(
                    "Could not find mandatory tab separated fields 'sample' & 'condition' in sample info file")
            for line in sampleinfo[1:]:
                linelist = line.split('\t')
                if not linelist[1].lower() in ['con', 'pat']:
                    sys.exit("Error for {}:\nThe condition of a sample should be either con or pat (case insensitive)".format(
                        linelist[0]))
                sampleDict[linelist[0]] = linelist[1].lower()  # mapping sample : phenotype
        for sample in simplifiedsamples:  # Check if all samples found in variantfile are present in sample info file
            if not sample in sampleDict.keys():
                sys.exit(
                    "Error: could not find sample {} (present in variantfile) in sample info file".format(sample))
    else:
        sampleDict = {sample: 'notspecified' for sample in simplifiedsamples}
    return sampleDict


def writeOutput(infodict, input, sampleInfo):
    '''
    Writing the result to stdout. Looping over inputfile and getting info out of dict
    In the case a sample info file is presented, a more elaborate output is given.
    '''
    for line in infodict['comments']:
        print(line, end='')
    print('\t'.join(infodict['header']), end='\n')
    groups_spec = ['pat', 'con']
    groups_all = ['pat', 'con', 'notspecified']
    zygosities_v = ['t', 'm', 'c', 'o']  # zygosities considered 'variant'
    zygosities_all = ['r', 't', 'm', 'c', 'o', 'u']  # all zygosities
    for line in input:
        dataToAdd = []
        variantstatus = {  # nested dictionary for both patients and controls and all possible zygosities, and not specified for when ran without sample info
            'pat': {'t': [], 'u': [], 'r': [], 'm': [], 'o': [], 'c': []},
            'con': {'t': [], 'u': [], 'r': [], 'm': [], 'o': [], 'c': []},
            'notspecified': {'t': [], 'u': [], 'r': [], 'm': [], 'o': [], 'c': []}
        }
        linelist = line.split('\t')
        # loop over the indices which map to 'mode' field of samples, fill the nested dictionary
        for key in sorted(infodict['zyg'].keys()):
            sample = infodict['zyg'][key]  # this is the full name of the sample including ngs id
            # look up phenotype of sample using individual id
            phenotype = sampleInfo['_'.join(sample.split('_')[:-1])]
            variantstatus[phenotype][linelist[key]].append(sample)  # fill in dictionary
        # Create a sorted list of all sample IDs which have a zygosity from zygosities_v and as such considered to be variant. List comprehension creates nested list, which is flattened
        carriers = sorted([item for sublist in [variantstatus[group][zyg] for group, zyg in zip(
            sorted(groups_all * len(zygosities_v)), zygosities_v * len(groups_all))] for item in sublist])
        if len(carriers) > 0:
            dataToAdd.extend([str(len(carriers)), ', '.join(carriers)])
        else:
            dataToAdd.extend([str(0), '-'])
        if args.samples:
            # Get the counts of all categories in a list 'counts', in the following order (conforming the header)
            # [0]con_r [1]con_t [2]con_m [3]con_c [4]con_o [5]con_u [6]pat_r [7]pat_t [8]pat_m [9]pat_c [10]pat_o [11]pat_u
            counts = [len(variantstatus[group][zyg]) for group, zyg in zip(
                sorted(groups_spec * len(zygosities_all)), zygosities_all * len(groups_spec))]
            names = [', '.join(variantstatus[group][zyg]) for group, zyg in zip(
                sorted(groups_spec * len(zygosities_all)), zygosities_all * len(groups_spec))]
        else:
            counts = [len(variantstatus['notspecified'][zyg]) for zyg in zygosities_all]
            names = [', '.join(variantstatus['notspecified'][zyg]) for zyg in zygosities_all]
        dataToAdd.extend([str(item) for item in counts])
        if args.names:
            dataToAdd.extend(names)
        if args.stats:
            dataToAdd.extend(getStats(counts))
        print(line.replace('\n', '\t') + '\t'.join(dataToAdd))


def getStats(counts):
    '''
    Executed if var.stats is set using --stats flag.
    The obtained counts are first reorganized.
    '''
    statdata = [counts[2], counts[1], counts[0], counts[8], counts[7], counts[6]]
    if not sum(statdata) == 0:
        fisher_exact_p = str(fisheR(statdata[0], statdata[1],
                                    statdata[2], statdata[3], statdata[4], statdata[5])[0])
        hwe_p = str(haRdy(statdata[0] + statdata[3], statdata[1] +
                          statdata[4], statdata[2] + statdata[5])[0])
    else:
        fisher_exact_p = "NA"
        hwe_p = "NA"
    return [fisher_exact_p, hwe_p]


def makeHeader(headerline):
    '''
    Based on the options set the header will be differentially modified to contain the calculated values
    Header is a list, to which extensions are made depending on flags set.
    '''
    initialheader = [word for word in headerline.replace('\n', '').split('\t')]
    initialheader.extend(['carriernumber', 'carriers'])
    if args.samples:
        initialheader.extend(['count_con_r', 'count_con_t', 'count_con_m', 'count_con_c', 'count_con_o', 'count_con_u',
                              'count_pat_r', 'count_pat_t', 'count_pat_m', 'count_pat_c', 'count_pat_o', 'count_pat_u'])
        if args.names:
            initialheader.extend(['samples_con_r', 'samples_con_t', 'samples_con_m', 'samples_con_c', 'samples_con_o', 'samples_con_u',
                                  'samples_pat_r', 'samples_pat_t', 'samples_pat_m', 'samples_pat_c', 'samples_pat_o', 'samples_pat_u'])
    else:
        initialheader.extend(['count_r', 'count_t', 'count_m', 'count_c', 'count_o', 'count_u'])
        if args.names:
            initialheader.extend(['samples_r', 'samples_t', 'samples_m',
                                  'samples_c', 'samples_o', 'samples_u'])
    if args.stats:
        initialheader.extend(['fisher_exact_p', 'hwe'])
    return initialheader


def getArgs():
    parser = argparse.ArgumentParser(
        description="Add counts and statistics to a .tsv variant file.", prog="Toolbox varstats")
    parser.add_argument("target",
                        help="Variant in tsv file to which information should be added. Use <stdin> for input on stdin.",
                        default='stdin')
    parser.add_argument("--samples",
                        help="A sample info file containing columns sample and condition (which should be pat or con, case insensitive).",
                        type=argparse.FileType('r'))
    parser.add_argument("--stats",
                        help="Statistics (fisher exact test and hwe) to be added to the file, requires --samples to be set. Requires functional rpy2 installation which might be an issue.",
                        action='store_true')
    parser.add_argument('--caller',
                        help="Counting should be performed based on either <gatk>, <sam> or <both>",
                        default='gatk',
                        choices=['gatk', 'sam', 'both'])
    parser.add_argument('--names',
                        help="Add all names to the result in different subgroups",
                        action="store_true")
    args = parser.parse_args()
    if args.stats and not args.samples:
        sys.exit("Using --stats requires --samples and a sample info file.")
    return args


def makeRFunctions():
    '''
    Initiates R functions
    Using rpy2 make a call to R functions for fisher exact test and hardy weinberg equilibrium test, returning the pvalues
    Functions are defined as variables in global scope.
    rpy2 dependency is problematic, although can be solved reasonably using conda
    See also: https://github.com/wdecoster/problemsolved/blob/master/Installingrpy2usingconda.md
    '''
    try:
        import readline
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        importr("HardyWeinberg")
        fisheR = robjects.r(
            "function(conhom, conhet, conref, pathom, pathet, patref){return(fisher.test(matrix(c(conhom, conhet, conref, pathom, pathet, patref), nrow = 2, byrow=TRUE))['p.value']$p.value)}")
        haRdy = robjects.r(
            "function(hom, het, ref){return(HWExact(c(hom, het, ref), verbose=F)$pval)}")
    except:
        sys.exit("Fatal error when attempting to load required modules for statistics.\nYou either don't have an rpy2 installation or it's corrupt.\nTry again without --stats")
    return fisheR, haRdy


def main():
    input = getInput()
    infodict = extractInfoDict(input, 'zyg', args.caller)
    sampleInfo = extractSampleDict(infodict['zyg'].values())
    writeOutput(infodict, input, sampleInfo)


if __name__ == "__main__":
    args = getArgs()
    if args.stats:
        fisheR, haRdy = makeRFunctions()
    main()
