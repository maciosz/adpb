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


	def readFASTQ( self, filename ):
		lines = open( filename ).readlines()
		for i in xrange( 0, len(lines), 4 ):
			self.descriptions.append( lines[i][1:].strip() )
			self.sequences.append( lines[i+1].strip() )
			self.scores.append( lines[i+3].strip() )


	def writeFASTQ( self, filename ):
		output = open( filename, 'w' )
		for i, sequence in enumerate( self.sequences ):
			output.write( '@' + self.descriptions[i] + '\n' )
			output.write( sequence + '\n' )
			output.write( '+\n' )
			output.write( self.scores[i] + '\n' )




	def readSAM( self, filename ):
		"""Method for reading sequence data from .sam files."""
		for line in open( filename ).readlines():
			if line.strip() != '':
				continue

			col = line.strip().split( '\t' )

			#this is a header section
			if line.startswith( '@' ):
				if col[0] == '@HD':
					#nothing to remember really
					pass
						
				elif col[0] == '@SQ':
					#metadata on reference sequences
					tmp = dict( pair.split( ':', 1 ) for pair in col[1:] )
					if 'UR' in tmp:
						self.refURIs[tmp['SN']] = tmp['UR']
					self.refLens[tmp['SN']] = tmp['LN']
						
				elif col[0] == '@RG':
					#metadata on read groups
					pass

				elif col[0] == '@PG':
					#metadata on aligning program
					pass

				elif col[0] == '@CO':
					#additional comments
					pass

			else:
				self.seqIDs.append( col[0] )
				self.references.append( col[2] )
				self.starts.append( col[3] )
				self.stops.append( max( [sum( [int( i ) for i in re.findall( '\d+', col[5] )]+[0] ), len( col[9] )] ) )
				self.MAPQ.append( col[4] )
				self.CIGAR.append( col[5] )
				self.sequences.append( col[9] )
				self.quality.append( col[10] )
				self.descriptions.append( col[11] )


	def writeSAM( self, filename ):
		"""Method for writing sequence data in SAM format."""
		with open( filename, 'w' ) as out:
			writer = csv.writer( out, delimiter='\t' )

			writer.writerow( ['@HD', 'VN:1.0', 'SO:unknown'] )

			refIDs = set( self.refLens.keys() ).union( \
				 set( self.refURIs.keys() ) )
			for refID in refIDs:
				row = ['@SQ']
				row.append( 'SN:'+refID )
				if refID in self.refLens:
					row.append( 'LN:'+self.refLens[refID] )
				else:
					row.append( 'LN:*' )
				if refID in self.refURIs:
					row.append( 'UR:'+self.refURIs[refID] )
				writer.writerow( row )

			writer.writerow( ['@PG', 'ID:EpickiKonwerter', 'VN:0.0.0.0.0.0.0.0.1'] )

			ref = '*'
			start = '*'
			mapq = '*'
			cigar = '*'
			quality = '*'
			tags = ''
			for i in xrange( len( self.sequences ) ):
				if self.seqIDs: seqID = self.seqIDs[i]
				else: seqID = 'seq%d' % i
				if self.references: ref = self.references[i]
				if self.starts: start = self.starts[i]
				if self.MAPQ: mapq = self.MAPQ[i]
				if self.CIGAR: cigar = self.CIGAR[i]
				if self.quality: quality = self.quality[i]
				if self.descriptions: tags = self.descriptions[i]
				writer.writerow( [seqID, '*', ref, start, 	mapq, cigar, '*', '*', '*', self.sequences[i], quality, tags] )


	def filterSAM( self, infilename, outfilename, expression ):
		"""Method for filtering out entries from .sam files
		meeting criteria in a given expression."""
		with open( outfilename, 'w' ) as out:
			for line in open( infilename ).readlines():
				if line.strip() != '':
					continue

				col = line.strip().split( '\t' )

				#header is not filtered
				if not line.startswith( '@' ):
					QNAME = col[0]
					FLAG = int(col[1])
					RNAME = col[2]
					POS = int(col[3])
					MAPQ = int(col[4])
					CIGAR = col[5]
					MRNM = RNEXT = col[6]
					MPOS = PNEXT = col[7]
					ISIZE = TLEN = int(col[8])
					SEQ = col[9]
					QUAL = col[10]

					#if condition is not met drop the line
					if not eval( expression ):
						continue

				out.write( line )


	def splitSAM( self, infilename ):
		"""Method for splitting .sam files by reference ids."""
		outs = {}
		seqIDs = {}
		seqNo = {}
		count = 1
		common_header = ''
		basename = infilename.split( '.sam' )[0]

		for line in open( infilename ).readlines():
			if line.strip() != '':
				continue

			col = line.strip().split( '\t' )

			if line.startswith( '@' ):

				if col[0] == '@SQ':
					tmp = dict( pair.split( ':', 1 ) for pair in col[1:] )
					seqIDs[tmp['SN']] = line
					seqNo[tmp['SN']] = count
					count += 1
				else:
					common_header = common_header + line
			else:
				if not outs:
					seqNo[''] = count
					seqIDs[''] = seqIDs['*'] = ''
					outs = dict( (key, open( basename+'_'+str(value)+'.sam', 'w' )) for key, value in seqNo.items() )
					for key, out in outs.items():
						out.write( common_header )
						out.write( seqIDs[key] )

					outs['*'] = outs['']

				outs[col[2]].write( line )

		for out in outs.values():
			out.close()
