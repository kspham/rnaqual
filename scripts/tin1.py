'''
Calculate transcript integrity number (TIN) for each transcript in BED file
Only consider exon-exon junctions
'''

import sys,os
import math
from optparse import OptionParser
from qcmodule import SAM
from qcmodule import BED
from qcmodule import mystat
from qcmodule import getBamFiles
from numpy import mean,median,std

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012. All rights reserved."
__credits__ = []
__license__ = "GPL"
__version__="2.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('class.log','a')
	print >>sys.stderr,mesg
	print >>LOG,mesg


def main():

	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input file(s) in BAM format. "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files. 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam file (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools.  [required]')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. [required]")
 
	(options,args)=parser.parse_args()
                
	if not (options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)
	(v_major, v_minor) = sys.version_info[0:2]
	if v_major != 2 or v_minor <7:
		print >>sys.stderr, "This program needs python 2.7.*"

	bamfiles = getBamFiles.get_bam_files(options.input_file)
	for ff in (bamfiles):
		if not os.path.exists(ff):
			print >>sys.stderr, '\n\n' + ff + " does NOT exists" + '\n'
			sys.exit(0)

	print >>sys.stderr, "Total %d BAM files:" % len(bamfiles)
	for f in sorted(bamfiles):
		print >>sys.stderr, "\t" + f
    
	
    
	print "\t".join(['Bam_file','TIN(mean)', 'TIN(median)','TIN(stdev)'])
	for f in sorted(bamfiles):
		OUT = open(os.path.basename(f).replace('bam','') + 'tin.xls','w')
		print >>OUT, "\t".join(["geneID","chrom", "tx_start", "tx_end","entropy","TIN"])
		tin_values=[]
		sam_obj = SAM.ParseBAM(f)
		bed_obj = BED.ParseBED(options.ref_gene_model)
		for gname,i_chr, i_tx_start, i_tx_end, junctions  in bed_obj.getSpliceJunctions():
			if len(junctions) == 0:	#singel exon gene
				entropy = "NA"
				entropy_bar = "NA"
				tin = "NA"
			freq = sam_obj.junction_freq(i_chr, i_tx_start, i_tx_end, junctions)
			entropy = mystat.shannon_entropy(freq.values())
			if entropy == "NA":
				tin = "NA"
			else:
				tin = 10*(math.exp(entropy) / len(junctions))
				tin_values.append(tin)
			print >>OUT, '\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, entropy, tin)])
		OUT.close()
			
		print "\t".join( [str(i) for i in (os.path.basename(f), mean(tin_values), median(tin_values), std(tin_values))])


if __name__ == '__main__':
	main()
