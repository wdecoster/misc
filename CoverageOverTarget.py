#!/complgen/bin/anaconda/bin/python2.7
#Script to obtain coverages and averages for genomic intervals for one or multiple samples
#Assumes bam files are indexed and sorted
#First argument is interval file, following argument(s) is one or more bam files, assumed to be sorted and indexed
#wdecoster

import sys, os, multiprocessing, collections, copy, glob
##### BAD WORKAROUND FOR AVOIDING IMPORT OF WRONG PYTHON MODULE
for dir in ['/storage/breva/complgen/bin/anaconda/lib/python2.7/site-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg', '/complgen/bin/anaconda/lib/python2.7/site-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg']:
	if dir in sys.path:
		sys.path.remove(dir)
import pysam

def checkbams(filelist):
	'''For now, only checks file extension and existance, limited testing, actually ignoring all errors'''
	if filelist[0].lower() == "seqcap":
		bams = [file for file in glob.glob('/complgen/bin/MastrDesign/db/SeqCapBams/*.bam')]
		if not len(bams) == 0:
			print("Using a standard set of {} bam files from SeqCap exomes.".format(len(bams)))
			return bams
		else:
			sys.exit("Error when attempting to grab my standard set of SeqCap bams, could not find a single bam.")
	bams = [file for file in filelist if os.path.isfile(file) and file[-4:] == ".bam"]
	indexfiles = [file for file in filelist if os.path.isfile(file + '.bai') and file[-4:] == ".bam"]
	if not len(indexfiles) == len(bams):
		sys.exit("Could not find .bai index file for {} bam files and this is required. Use samtools index.".format(len(bams)-len(indexfiles)))
	if len(bams) > 0:
		print("Found {} valid bamfiles for analysis.".format(len(bams)))
		return bams
	else:
		sys.exit("Argument Error: Could not find a valid bam file. I naively checked for existance and file extension and did not perform any more elaborate checks.")

def checkintervals(intervalfile):
	'''Check validity of supplied intervalfile and intervals, returning those as a list of tuples'''
	if os.path.isfile(intervalfile):
		with open(intervalfile) as intervalf:
			intervallines = [line.strip() for line in intervalf.readlines() if not line == '\n']
			if not len(intervallines) > 1:
				sys.exit("Maybe unproperly formatted intervalfile, could not find a single interval or header is missing.")
			intervallist = []
			targetdict = {}
			if not intervallines[0] == "chromosome\tbegin\tend\tname":
				print(intervallines[0])
				sys.exit("Could not find the required header of the intervalfile 'chromosome\tbegin\tend\tname'")
			else:
				for line in intervallines[1:]:
					interval = line.split('\t')
					if not len(interval) == 4:
						sys.exit("Insufficient information on line with name {}, could not find 4 fields.".format(interval[3]))
					try:
						if not int(interval[2]) - int(interval[1]) > 0:
							sys.exit("On line with name {} the end is smaller or equal to the begin, invalid interval.".format(interval[3]))
					except TypeError:
						sys.exit("Could not find an integer where expecting begin or end on line with name {}".format(interval[3]))
					intervallist.append((interval[3], interval[0], int(interval[1]), int(interval[2]))) #When the interval survived all checks it is appended to the intervallist as a tuple with elements name, chromosome, int(begin), int(end)
					targetdict[interval[3]] = (interval[0], int(interval[1]), int(interval[2]))
				if not len(targetdict.keys()) == len(set(targetdict.keys())):
					sys.exit("Make sure to use unique identifiers as names in the targetfile.")
				return(intervallist, targetdict) #Returning both in list and in dict form for convenience, perhaps to be reformated later
	else:
		sys.exit("Path specified to intervalfile is incorrect.")

def getAverages(coveragelist, targets):
	averages = {'sample' : 'Averages'}
	for target in targets.keys():
		totcov = 0
		for sample in coveragelist:
			totcov += sample[target]
		averages[target] = totcov / len(coveragelist)
	return(averages)

def makeCoverageDicts(bam):
	'''For each sample, a dict is created which contains key sample with value samplename and a dictionary key per target with value: list of coverages in the interval'''
	workfile = pysam.AlignmentFile(bam, "rb")
	coverages = {'sample': os.path.basename(bam)[:-4].replace('map-rdsbwa-', '')} #Stripping off the extension to get the samplename
	for targetline in intervals: #Using intervals from the global scope
		coverage = []
		for pileupcolumn in workfile.pileup(targetline[1], targetline[2], targetline[3]): #Total coverage over targeted region of each base
			if pileupcolumn.pos >= targetline[2] and pileupcolumn.pos < targetline[3]:
				coverage.append(pileupcolumn.nsegments)
		coverages[targetline[0]] = coverage
	return coverages

def makeRegAverages(coverageslist, targets):
	localcoverageslist = copy.deepcopy(coverageslist) #Need a deep copy in order not to influence the objects for the rest of the processing
	for sample in localcoverageslist:
		for target in targets.keys(): #Since 'sample' is not in targets.keys() this dict-item is skipped
			sample[target] = sum(sample[target])/ float(targets[target][2] - targets[target][1]) #Calculating average
	return localcoverageslist

def perNucleotide(coveragelist):
	for sample in coveragelist:
		for target in targets.keys():
			assert len(sample[target]) <= float(targets[target][2] - targets[target][1]) + 1
			sample[target] = sum([cov < 30 for cov in sample[target]]) / float(targets[target][2] - targets[target][1]) #cutoff from global scope, calculate fraction
	return coveragelist
	
def writeOut(coveragelist, targets, name, round=True):
	with open(name, 'w') as output:
		output.write('Interval\t{}\n'.format('\t'.join([entry['sample'] for entry in coveragelist]))) #List expression generates sample list in order as in coveragelist
		for target in targets:
			if round:
				output.write('{}\t{}\n'.format(target, '\t'.join([str(int(entry[target])) for entry in coveragelist]))) #List expression generates coverageslist in same order as samples
			else:
				output.write('{}\t{}\n'.format(target, '\t'.join([str(entry[target]) for entry in coveragelist]))) #List expression generates coverageslist in same order as samples

def histtrack(coveragelist, targets):
	posaveragedict = {}
	with open("PerNucleotideHistogram.txt", 'w') as outfile:
		outfile.write("browser position chr12:40,618,813-40,763,086\n track name=CoverageHist alwaysZero=on type=bedGraph description=CoverageHist visibility=2 color=255,0,0\n")
		for target in targets:
			perNuclList = zip([sample[target] for sample in coveragelist])
			posaveragedict[target] = [sum(position)/float(len(pos)) for position in perNuclList]
	
if len(sys.argv) >= 3:
	bamlist = checkbams(sys.argv[2:])
	intervals, targets = checkintervals(sys.argv[1])
	print("Preparing per amplicon statistics.")
	pool = multiprocessing.Pool(processes=16)
	res = pool.map_async(makeCoverageDicts, bamlist)
	coveragelist = (res.get())
	#histtrack(coveragelist, targets)
	assert len(coveragelist) == len(bamlist), "I somewhere lost a sample. Shame."
	regAverages = makeRegAverages(coveragelist, targets)
	if len(bamlist) > 1:
		regAverages.append(getAverages(regAverages, targets)) #adding 'averages' as being another sample
	writeOut(regAverages, targets, "CoveragesOverTargets.tsv")
	cutoff = 30
	print("Preparing per nucleotide statistics, arbitrary cut-off is {}x coverage.".format(cutoff))
	nuclCov = perNucleotide(coveragelist)
	assert len(nuclCov) == len(bamlist), "I somewhere lost a sample. Shame."
	if len(bamlist) > 1:
		nuclCov.append(getAverages(nuclCov, targets)) #adding 'averages' as being another sample
	writeOut(nuclCov, targets, "NucleotidesBelow" + str(cutoff) + "x.tsv", round=False)
else:
	sys.exit("Argument Error: Required arguments are a intervalfile with header chromosome begin end name and one or multiple bam files.")