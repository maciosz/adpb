import re
import csv

class EpickiKonwerter:

	def __init__( self ):
		self.sequences = []
		self.seqIDs = []
		self.references = []
		self.refURIs = {}	#dictionary {refID: URI}
		self.refLens = {}	#dictionary {refID: length}
		self.sources = []
		self.methods = []
		self.starts = []
		self.stops = []
		self.scores = []
		self.strands = []
		self.phases = []
		self.groups = []	# lista slownikow
		self.descriptions = []
		self.MAPQ = []
		self.CIGAR = []
		self.quality = []

	def readFASTA( self, filename ):
		lines = open( filename ).readlines()
		for line in lines:
			if line.startswith('>'):
				self.descriptions.append( line[1:].strip() )
				self.sequences.append('')
			else:
				self.sequences[-1] += line.strip()



	def writeFASTA( self, filename ):
		return 0


	def readSAM( self, filename ):
	"""Method for reading sequence data from .sam files."""
		for line in open( filename ).readlines():
			line = line.strip().split('\t')

			#this is a header section
			if line[0].startswith('@'):
				tmp = dict(pair.split(':', 1) for pair in line[1:])
				if line[0] == '@HD':
					#nothing to remember really
						
				elif line[0] == '@SQ':
					#metadata on reference sequences
					uri = ''
					if 'UR' in tmp.keys(): uri = tmp['UR']
					references[tmp['SN']] = uri
					refLens[tmp['SN']] = tmp['LN']
						
				elif line[0] == '@RG':
					#metadata on read groups

				elif line[0] == '@PG':
					#metadata on aligning program

				elif line[0] == '@CO':
					#additional comments

			else:
				self.seqIDs.append(line[0])
				self.references.append(line[2])
				self.starts.append(line[3])
				self.ends.append(max([sum([int(i) for i in re.findall('\d+', c)]+[0]), len(line[9])]))
				self.MAPQ.append(line[4])
				self.CIGAR.append(line[5])
				self.sequences.append(line[9])
				self.quality.append(line[10])
				self.descriptions.append(line[11])

	def writeSAM( self, filename ):
	"""Method for writing sequence data in SAM format."""
		with open( filename, 'w' ) as out:
			writer = csv.writer(out, delimiter='\t')

			writer.writerow(['@HD\tVN:1.0\tSO:unknown'])

			refIDs = set(self.refLens.keys()).union( \
				 set(self.refURIs.keys())
			for refID in refIDs:
				row = ['@SQ']
				if refID in refLens:
					row.append('LN:'+self.refLens[refID])
				else:
					row.append('LN:*')
				if refID in refURIs:
					row.append('UR:'+self.refURIs[refID])
					writer.writerow(row)

			writer.writerow(['@PG', 'ID:EpickiKonwerter', 'VN:0.0.0.0.0.0.0.0.1'])

			ref = '*'
			start = '*'
			mapq = '*'
			cigar = '*'
			quality = '*'
			tags = ''
			for i in xrange(len(self.sequences)):
				if self.seqIDs: seqID = self.seqIDs[i]
				else: seqID = 'seq%d' % i
				if self.references: ref = self.references[i]
				if self.starts: start = self.starts[i]
				if self.MAPQ: mapq = self.MAPQ[i]
				if self.CIGAR: cigar = self.CIGAR[i]
				if self.quality: quality = self.quality[i]
				if self.descriptions: tags = self.descriptions[i]
				writer.writerow([seqID, '*', ref, start, 	mapq, cigar, '*', '*', '*', self.sequences[i], quality, tags])
