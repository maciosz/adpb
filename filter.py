#! /usr/bin/python

import EpickiKonwerter
import argparse

def main():

	parser = argparse.ArgumentParser( description='Filter sequence file according to a given expression.' )
	parser.add_argument( '-i', metavar='<inputfile>', type=str, \
			     required=True, nargs=1, help='input file path' )
	parser.add_argument( '-o', metavar='<outputfile>', type=str, \
			     required=True, nargs=1, help='output file path' )
	parser.add_argument( '-f', metavar='<format>', type=str, \
			     required=True, nargs=1, \
			     choices=['fasta', 'fastq', 'gff2', 'gff3', 'sam'], \
			     help='input and output file format' )
	parser.add_argument( '-e', metavar='<expression>', type=str, \
			     required=True, nargs=1, \
			     help='expression to be evaluated' )

	args = parser.parse_args()

	filterFile( args.i[0], args.o[0], args.f[0], args.e[0] )

	
def filterFile( inFile, outFile, IOFormat, expression ):

	e = EpickiKonwerter.EpickiKonwerter()

	filters = dict( zip( ['fasta', 'fastq', 'gff2', 'gff3', 'sam'], \
			    [e.filterFASTA, e.filterFASTQ, e.filterGFF, e.filterGFF, e.filterSAM] ) )

	filters[ IOFormat ]( inFile, outFile, expression )



if __name__ == "__main__":
	main()
