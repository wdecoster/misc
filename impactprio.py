#!/complgen/bin/anaconda/bin/python2.7
#Script to parse and prioritize what's in the refGene_impact field, writes to stdout

import sys, os, time

impactlist = ["GENEDEL", "GENESUB", "CDSSTART", "CDSDEL", "CDSSUB", "CDSSTOP", "CDSNONSENSE", "CDSFRAME",
	"CDSSUBSTART", "CDSDELSTART", "CDSINS", "ESPLICE", "CDSSUBSPLICE", "CDSDELSPLICE", "UTR5SPLICE", "UTR3SPLICE",
	"CDSMIS", "CDSsilent", "splice", "RNAEND", "RNASPLICE", "RNASTART", "RNA", "UTR5", "UTR3", "intron", "upstream", "downstream", "annot", ""]

def getImpactField(inputhandle):
	for line in inputhandle:
		if line.startswith('#'):
			continue
		elif line.startswith('chr'):
			headerline = line.strip().split('\t')
			try:
				impactfield = headerline.index('refGene_impact')
			except ValueError:
				sys.exit("INPUT ERROR: could not find the field 'refGene_impact' in the compar file.")
			break
		else:
			sys.exit("Header not recognized, expected header starting with 'chr' but found line[:15]...")
	return(headerline, impactfield)

def main(comparinfile):
	root, ext = os.path.splitext(comparinfile)
	with open(comparinfile) as inputcompar:
		headerline, impactfield = getImpactField(inputcompar)
		print("{}\t{}".format('\t'.join(headerline), 'refGene_priority_impact'))
		for line in inputcompar:
			varimpacts = list(set(line.strip().split('\t')[impactfield].split(';')))
			try:
				impactsdamage = min([impactlist.index(impact) for impact in varimpacts])
			except ValueError:
				for impact in varimpacts:
					if impact not in impactlist:
						sys.stderr.write("WARNING: unrecognized impact value: '{}', contact Wouter for your exotic impact to get it implemented.\n".format(impact))
				print("{}\t{}".format(line.strip(), "unrecognized"))
			print("{}\t{}".format(line.strip(), impactlist[impactsdamage]))
			
def givehelp():
	print("Prioritization of impact fields, which is done in the following order from damaging to pretty much okay:")
	for i, elem in enumerate([im for im in impactlist if not im in ["", "annot"]]):
		print("{}\t{}".format(i, elem))
	print("Inform me if you disagree.")

if os.path.isfile(sys.argv[1]):
	main(sys.argv[1])
elif sys.argv[1] == "help":
	givehelp()
else:
	sys.exit("USER ERROR: unrecognized input or path to inputfile incorrect.\nUse 'help' as argument to get priority list.")