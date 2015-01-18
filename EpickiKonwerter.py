class EpickiKonwerter:

	def __init__( self ):
		self.sequences = []
		self.references = []
		self.sources = []
		self.methods = []
		self.starts = []
		self.stops = []
		self.scores = []
		self.strands = []
		self.phases = []
		self.groups = []	# lista slownikow
		self.descriptions = []

	def readFASTA( self, filename ):
		lines = open( filename ).readlines()
		for line in lines:
			if line.startswith('>'):
				self.descriptions.append( line[1:].strip() )
				self.sequences.append('')
			else:
				self.sequences[-1] += line.strip()


	def splitSequence( self, sequence, length=80 ):
		list_of_sequences = []
		for i in xrange( len(sequence) / length ):
			list_of_sequences.append( sequence[ i*length : (i+1)*length ] )
		if len(sequence) % length != 0:
			list_of_sequences.append( sequence[ (i+1)*length : ] )
		return list_of_sequences


	def writeFASTA( self, filename ):
		output = open( filename, 'w' )
		for i, sequence in enumerate( self.sequences ):
			output.write( '>' + self.descriptions[i] + '\n' )
			sequences_to_write = '\n'.join( self.splitSequence(sequence) )
			output.write( sequences_to_write )
			output.write('\n')


