#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
This program is used estimate clipping profile of RNA-seq reads from BAM or SAM file
Note that to use this funciton, CIGAR strings within SAM/BAM file should have 'S' operation 
(This means your reads mapper should support clipped mapping)
-------------------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
from time import strftime
import subprocess
 

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s).")
	(options,args)=parser.parse_args()

	if not (options.input_file):
		parser.print_help()
		sys.exit(0)
	for input_file in ([options.input_file]):
		if not os.path.exists(input_file):
			print >>sys.stderr, '\n\n' + input_file + " does NOT exists" + '\n'
			#parser.print_help()
			sys.exit(0)

	obj = SAM.ParseBAM(options.input_file)
	obj.clipping_profile(outfile=options.output_prefix)
	try:
		subprocess.call("Rscript " + options.output_prefix + '.clipping_profile.r',shell=True)
	except:
		print >>sys.stderr, "Cannot generate pdf file form " + options.output_prefix + '.clipping_profile.r'
		pass

if __name__ == '__main__':
	main()
