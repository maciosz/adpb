#! /usr/bin/python

import EpickiKonwerter
import argparse

def main():

	parser = argparse.ArgumentParser( description='Convert sequence file formats with optional filtering.' )
	parser.add_argument( '-i', metavar='<inputfile>', type=str, \
			     required=True, nargs=1, help='input file path' )
	parser.add_argument( '-o', metavar='<outputfile>', type=str, \
			     required=True, nargs=1, help='output file path' )
	parser.add_argument( '-f', metavar='<format>', type=str, \
			     required=True, nargs=2, \
			     choices=['fasta', 'fastq', 'gff2', 'gff3', 'sam'], \
			     help='input and output file formats' )
	parser.add_argument( '-e', metavar='<expression>', type=str, \
			     default=[False], nargs=1, \
                             help='expression to be evaluated' )

	args = parser.parse_args()

	convertFile( args.i[0], args.o[0], args.f[0], args.f[1], args.e[0] )

	
def convertFile( inFile, outFile, inFormat, outFormat, expression=False ):

	e = EpickiKonwerter.EpickiKonwerter()

	reader = dict( zip( ['fasta', 'fastq', 'gff2', 'gff3', 'sam'], \
			    [e.readFASTA, e.readFASTQ, e.readGFF2_GFF3, e.readGFF2_GFF3, e.readSAM] ) )
	writer = dict( zip( ['fasta', 'fastq', 'gff2', 'gff3', 'sam'], \
			    [e.writeFASTA, e.writeFASTQ, e.writeGFF2, e.writeGFF3, e.writeSAM] ) )

	reader[ inFormat ]( inFile )
	if expression:	e.filter( expression )
	writer[ outFormat ]( outFile )



if __name__ == "__main__":
	main()
