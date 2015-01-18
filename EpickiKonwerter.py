import re

class EpickiKonwerter:

	def __init__( self ):
		self.sequences = []
		self.SAMsequences = {}	#list of pairs (seqID, seq)
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
				self.sequences.append(line[9])
				self.references.append(line[2])
				self.starts.append(line[3])
				self.ends.append(max([sum([int(i) for i in re.findall('\d+', c)]+[0]), len(line[9])]))
				self.MAPQ.append(line[4])
				self.CIGAR.append(line[5])
				self.quality.append(line[10])
				self.descriptions.append(line[11])

	def writeSAM( self, filename ):
	"""Method for writing sequence data in SAM format."""

