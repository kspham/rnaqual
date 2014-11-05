#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Calculate transcript integrity number (TIN) for each transcript in BED file.
-------------------------------------------------------------------------------------------------'''


import sys,os
import math
from optparse import OptionParser
from qcmodule import SAM
from qcmodule import BED
from qcmodule import mystat
from qcmodule import getBamFiles
from numpy import mean,median,std
from time import strftime
import pysam

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012. All rights reserved."
__credits__ = []
__license__ = "GPL"
__version__="2.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


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
	print >>sys.stderr,mesg

def genebody_percentile(refbed):
	'''
	return percentile points of gene body
	'''
	if refbed is None:
		print >>sys.stderr,"You must specify a bed file representing gene model\n"
		exit(0)
	
	transcript_count = 0
	for line in open(refbed,'r'):
		try:
			if line.startswith(('#','track','browser')):continue  
			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
			geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
				
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
			transcript_count += 1
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue
		gene_all_base=[]
		mRNA_len =0
		flag=0
		for st,end in zip(exon_starts,exon_ends):
			#print chrom,st,end
			gene_all_base.extend(range(st+1,end+1))		#1-based coordinates on genome
		
		yield (geneName, chrom, tx_start, tx_end, mystat.percentile_list(gene_all_base))	#get 100 points from each gene's coordinates


def genebody_coverage(samfile, chrom, positions):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''

	cvg = []
	'''
	for i in positions:
		st = i -1
		end =i
		try:
			for pileupcolumn in samfile.pileup(chrom,st, end, truncate=True):
				if pileupcolumn.n == 0:
					cvg.append(0)
					continue
				cover_read = 0
				for pileupread in pileupcolumn.pileups:
					if pileupread.is_del: continue
					if pileupread.alignment.is_qcfail:continue 
					if pileupread.alignment.is_secondary:continue 
					if pileupread.alignment.is_unmapped:continue
					if pileupread.alignment.is_duplicate:continue
					cover_read +=1
				cvg.append(cover_read)
		except:
			cvg.append(0)	
	return cvg		
	'''	
	try:
		for pileupcolumn in samfile.pileup(chrom, positions[0]-1, positions[-1], truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions: continue
			if pileupcolumn.n == 0:
				cvg.append(0)
				continue				
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read +=1
			cvg.append(cover_read)
	except:
		cvg = []
	return cvg

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input file(s) in BAM format. "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files. 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam file (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools. [required]')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. [required]")
	(options,args)=parser.parse_args()
	
	if not (options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print >>sys.stderr, '\n\n' + options.ref_gene_model + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)
	
	printlog("Read BED file (reference gene model) ...")
	gene_percentiles = genebody_percentile(refbed = options.ref_gene_model)
		
	printlog("Get BAM file(s) ...")
	bamfiles = sorted(getBamFiles.get_bam_files(options.input_files))
	
	print >>sys.stderr, "Total %d BAM files:" % len(bamfiles)
	for f in bamfiles:
		print >>sys.stderr, "\t" + f
	
	print "\t".join(['Bam_file','TIN(mean)', 'TIN(median)','TIN(stdev)'])
	for f in bamfiles:
		OUT = open(os.path.basename(f).replace('bam','') + 'tin.xls','w')
		print >>OUT, "\t".join(["geneID","chrom", "tx_start", "tx_end","entropy","TIN"])
		tin_values=[]
		samfile = pysam.Samfile(f, "rb")
		
		for gname, i_chr, i_tx_start, i_tx_end, positions in genebody_percentile(options.ref_gene_model):
			if len(positions) == 0:
				entropy = "NA"
				entropy_bar = "NA"
				tin = "NA"
			
			coverage = genebody_coverage(samfile, i_chr,positions)
			entropy = mystat.shannon_entropy(coverage)
			if entropy == "NA":
				tin = "NA"
			else:
				tin = 10*(math.exp(entropy) / len(coverage))
				tin_values.append(tin)
			print >>OUT, '\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, entropy, tin)])
		OUT.close()
		print "\t".join( [str(i) for i in (os.path.basename(f), mean(tin_values), median(tin_values), std(tin_values))])

	
	
if __name__ == '__main__':
	
	#samfile = pysam.Samfile('/data2/bsi/staff_analysis/m102324/eRIN/PLOS_one/SRR873822_RIN10.bam', "rb")
	#p = [66942667, 66942668, 66942669,66942670]
	#a=genebody_coverage(samfile, 'chrX',p )
	#for i,j in zip(p,a):
	#	print i,j
	
	#genebody_percentile("/data2/bsi/staff_analysis/m102324/eRIN/PLOS_one/tmp")
	
	main()
