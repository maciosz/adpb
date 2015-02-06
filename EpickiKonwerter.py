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
		self.group_gff = []
		self.descriptions = []
		self.MAPQ = []
		self.CIGAR = []
		self.quality = []

	"""
	def openFile( self, filename ):
	"""
	#	Opens file, return list of lines.
	"""
		try:
			return open( filename ).readlines()
		except IOError:
			print "Can\'t open file " + filename + " for reading."
	"""

	def readFASTA( self, filename ):
		"""
		Reads file in FASTA format.
		"""
		
		lines = self.open_file( filename )
		for line in lines:
			if line.startswith('>'):
				self.descriptions.append( line[1:].strip() )
				self.sequences.append('')
			else:
				self.sequences[-1] += line.strip()


	def splitSequence( self, sequence, length=80 ):
		"""
		Splits sequence into list of sequences of given length.
		"""
		list_of_sequences = []
		i = -1
		for i in xrange( len(sequence) / length ):
			list_of_sequences.append( sequence[ i*length : (i+1)*length ] )
		if len(sequence) % length != 0:
			list_of_sequences.append( sequence[ (i+1)*length : ] )
		return list_of_sequences


	def getIthAttribute( self, i, list_of_attributes, default = ''):
		"""
		Returns i-th element of list_of_attributes, default if i-th element doesn't exist.
		"""
		try:
			return list_of_attributes[i]
		except IndexError:
			return default

	def writeFASTA( self, filename ):
		"""
		Writes file in FASTA format.
		"""
		output = open( filename, 'w' )
		for i, sequence in enumerate( self.sequences ):
			description = self.getIthAttribute(i, self.descriptions)
			output.write( '>' + description + '\n' )
			sequences_to_write = '\n'.join( self.splitSequence(sequence) )
			output.write( sequences_to_write )
			output.write('\n')


	def readFASTQ( self, filename ):
		"""
		Reads file in FASTQ format.
		"""
		lines = self.open_file( filename )
		for i in xrange( 0, len(lines), 4 ):
			if lines[i][0] == "@" and lines[i+2][0] == '+':
				self.descriptions.append( lines[i][1:].strip() )
				self.sequences.append( lines[i+1].strip() )
				self.scores.append( lines[i+3].strip() )
			else:
				print "File " + filename + " seems not to be in FASTQ format."


	def writeFASTQ( self, filename ):
		"""
		Writes file in FASTQ format.
		"""
		output = open( filename, 'w' )
		for i, sequence in enumerate( self.sequences ):
			description = self.getIthAttribute(i, self.descriptions)
			output.write( '@' + description + '\n' )
			output.write( sequence + '\n' )
			output.write( '+\n' )
			score = self.getIthAttribute( i, self.scores)
			output.write( score + '\n' )

	def filterFASTA( self, infilename, outfilename, length_range ):
		"""
		Reads file in FASTA format and filters out sequences outside the given range
		(given in format "10:100", 0 for unlimited ("0:100" means no longer than 100, "30:0" means no shorter than 30)),
		then writes to FASTA file
		"""
		minimum, maximum = interpretLengthRange( length_range )
		filterFASTAOrFASTQ( self, infilename, outfilename, minimum, maximum, "FASTA")
	
	def filterFASTQ( self, infilename, outfilename, length_range ):
		"""
		Reads file in FASTQ format and filters out sequences outside the given range
		(given in format "10:100", 0 for unlimited ("0:100" means no longer than 100, "30:0" means no shorter than 30)),
		then writes to FASTQ file
		"""
		minimum, maximum = interpretLengthRange( length_range )
		filterFASTAOrFASTQ( self, infilename, outfilename, minimum, maximum, "FASTQ")

	def interpretLengthRange( self, length_range ):
		minimum, maximum = length_range.split(':')
		if maximum == 0:
			maximum = 'a'
		return minimum, maximum



	def filterFASTAOrFASTQ( self, infilename, outfilename, minimum, maximum, fileformat ):
		def checkLength( i ):
			length = len( self.sequences[i] )
			return length >= minimum and length <= maximum
		if fileformat == "FASTA":
			self.readFASTA( infilename )
		if fileformat == "FASTQ":
			self.readFASTQ( infilename )
		for list_of_attributes in vars(self):
			if list_of_attributes == 'sequences':
				continue
			vars(self)[list_of_attributes] = filter( lambda x: checkLength(vars(self)[list_of_attributes].index(x)), vars(self)[list_of_attributes] )
		self.sequences = filter( lambda sequence: len(sequence) >=minimum and len(sequence) <= maximum, self.sequences )
		if fileformat == "FASTA":
			self.writeFASTA( outfilename )
		if fileformat == "FASTQ":
			self.writeFASTQ( outfilename )



	

	def readSAM( self, filename ):
		"""Method for reading sequence data from .sam files."""
		for line in open( filename ).readlines():
			if line.strip() == '':
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

				source = "*"
				if len( col ) == 10:
					desc = ''
				elif len( col ) == 11:
					desc = col[11]
				else:
					desc = ' '.join( col[11:] )

				source = "*"
				method = "*"
				score = "."
				strand = "."
				phase = "."
				group_gff = "*"
				if re.search( "CT", desc ):
					tmp = re.findall( "[.+/-];[A-Za-z0-9|%]+;", desc ).split( ";" )
					if tmp:
						strand = tmp[0]
						method = tmp[1]

					tmp = re.findall( "FSource=[A-Za-z0-9|%_\-]+;", desc )
					if tmp: source = tmp[0].split("=")[1][:-1]

					tmp = re.findall( "FScore=[0-9.]+;", desc )
					if tmp: score = tmp[0].split("=")[1][:-1]

					tmp = re.findall( "FPhase=[A-Za-z0-9|%_\-]+;", desc )
					if tmp: phase = tmp[0].split("=")[1][:-1]

					desc = "".join( re.split( "CT[A-Za-z0-9%|_\-.+:?$^&*@]+", desc ) )

				self.sources.append( source )
				self.methods.append( method )
				self.scores.append( score )
				self.strands.append( strand )
				self.phases.append( phase )
				self.group_gff.append( group_gff )
				self.descriptions.append( desc )


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
			for i in xrange( max( len( self.seqIDs ), len( self.sequences ) ) ):
				if self.seqIDs: seqID = self.seqIDs[i]
				else: seqID = 'seq%d' % i
				if self.references: ref = self.references[i]
				if self.starts: start = self.starts[i]
				if self.MAPQ: mapq = self.MAPQ[i]
				if self.CIGAR: cigar = self.CIGAR[i]
				if self.quality: quality = self.quality[i]
				if self.descriptions:
					tags = self.descriptions[i]

				if self.strands and self.methods:
					tmp = "CT:"
					tmp = tmp + self.strands[i] + ";" + self.methods[i]
					if self.scores:
						tmp = tmp + ";FScore=" + self.scores[i]

					if self.sources:
						tmp = tmp + ";FSource=" + self.sources[i]

					if self.phases:
						tmp = tmp + ";FPhase=" + self.phases[i]

					tags = tmp + " " + tags

				print "sekwencja nr ", i, " wyglada tak: "
				print self.sequences[i]

				writer.writerow( [seqID, '*', ref, start, mapq, cigar, \
						  '*', '*', '*', self.sequences[i], \
						  quality, tags] )


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

					#if condition is met write the line
					if eval( expression ):
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


	def open_file(self, filename):
		"""
		Opens a file.
		"""
		file1 = open(filename, "r")
		try: 
			data = file1.readlines()
		except IOError:
			print "Error occured during reading a file.\n"
		finally:
			file1.close()
    
		return data
      
      
	def save_file(self, filename, content):
		"""
		Writes a given content to a file.
		"""
		file1 = open(filename, "w")
		try: 
			file1.writelines(content)
		except IOError:
			print "Error occured during writing to file.\n"
		finally:
			file1.close()    
  

	def readGFF2_GFF3(self, filename):
		"""
		Reads file in gff2 or gff3 format.
		"""
		data = self.open_file(filename)    
        
		for i in data:
			if i.startswith("#"):
				self.headline.append(i)
			else:
				line = i.strip().split("\t")
	    
				self.seqIDs.append(line[0])
				self.sources.append(line[1])
				self.methods.append(line[2])
				self.starts.append(line[3])
				self.stops.append(line[4])
				if (line[5] == ".") or (float(line[5]) <= 1000) and (float(line[5]) > 0):
					self.scores.append(line[5])
				if (line[6] == "+") or (line[6] == "-") or (line[6] == "."):
					 self.strands.append(line[6])
				if (line[7] == "0") or (line[7] == "1") or (line[7] == "2") or (line[7] == "."):
					 self.phases.append(line[7])
            
				temp = []
				for j in range(8, len(line)):
					if self.headline[0] == "##gff-version 3\n":
						line[j] = line[j].split("=")
					else:
						line[j] = line[j].split(" ")
					for k in range(len(line[j])):
						if re.search(";", line[j][k]) != None:
							temp = line[j][k].split(";") 
							line[j].pop(k)
							line[j].insert(k, temp[0] + ";")
							line[j].insert(k+1, temp[1])
							temp = []
					for l in range(len(line[j])):
						if line[j][l].startswith(r"\""):    
							line[j][l] = line[j][l].translate(None, "\"")
	            
					not_empty = filter(None, line[j])
					self.group_gff.append(not_empty)
			line = []
            
		self.headline.pop(0)
        

	def writeGFF2(self, filename):
		"""
		Writes file in gff2 file format.
		"""
		out = []
		verse = [self.seqIDs, self.sources, self.methods, self.starts, self.stops, self.scores, self.strands, self.phases, self.group_gff]
		SOFA_ids = [ "SO:0000704", "SO:0000234", "SO:0000147", "SO:0000316", "SO:0000188", "SO:0000610", "SO:0000553", "SO:0000204", "SO:0000205"]
		types_gff3 = ["gene", "mRNA", "exon", "cds", "intron", "polyA", "polyA", "five", "three"]
		line = ""
		out.append("##gff-version 2\n")
		if self.headline != []:
			for k in self.headline:
				out.append(k)
		for i in range(len(self.seqIDs)):
			for j in range(len(verse)):
				if j == 2 and (re.search("SO:", verse[j][i]) != None): #self.methods
					for l in range(len(SOFA_ids)):
						if verse[j][i] == SOFA_ids[l]:
							line += types_gff3[l] + "\t"
				elif j == 8: # self.group_gff
					for k in range(len(self.group_gff[i])):
						if self.group_gff[i][k-1].endswith(";") and k != 0:
							if (re.search(" ", self.group_gff[i][k+1]) != None) or (re.search("_", self.group_gff[i][k+1]) != None):
								self.group_gff[i][k+1] = "\"" + self.group_gff[i][k+1] + "\""
						line += self.group_gff[i][k] + " "
				else:
					line += verse[j][i] + "\t"

			out.append(line+"\n")
			line = ""
        
		self.save_file(filename, out)
        
	    
	def writeGFF3(self, filename):
		"""
		 Writes file in gff3 file format.
		"""
		out = []
		verse = [self.seqIDs, self.sources, self.methods, self.starts, self.stops, self.scores, self.strands, self.phases, self.group_gff]
        
		out.append("##gff-version 3\n")
		if self.headline != []:
			for k in self.headline:
				out.append(k)
        
		line = ""
		for i in range(len(self.seqIDs)):
			for j in range(len(verse)):
				if j == 8:
					for k in range(len(self.group_gff[i])):
						if k == 0 or self.group_gff[i][k-1].endswith(";"):
							if self.group_gff[i][k] == "ID" or \
							   self.group_gff[i][k] == "Name" or \
							   self.group_gff[i][k] == "Alias" or \
							   self.group_gff[i][k] == "Parent" or \
							   self.group_gff[i][k] == "Target" or \
							   self.group_gff[i][k] == "Gap" or \
							   self.group_gff[i][k] == "Derives_from" or \
							   self.group_gff[i][k] == "Note" or \
							   self.group_gff[i][k] == "Dbxref" or \
							   self.group_gff[i][k] == "Ontology_term" or \
							   self.group_gff[i][k] == "Is_circular" :
								line += self.group_gff[i][k] + "="
							else:
								line += "Note="
						else:
							line += self.group_gff[i][k]
				else:
					line += verse[j][i] + "\t"

			out.append(line+"\n")
			line = ""
        
		self.save_file(filename, out)

    
	def filterGFF(self, filename_in, filename_out, phrase):
		"""
		Filters gffs files using eval.
		"""
		data = self.open_file(filename_in)
		out = [] 
		line = ""
		headline = 0
              
		for i in data:
			if i.startswith("#"):
				out.append(i)
				headline += 1
			else:
				line = i.strip().split("\t") 

				SEQID = line[0]
				SOURCE = line[1]
				METHOD = line[2]
                
				try:
					START = int(line[3])
				except:
					START = line[3]
				try:
					END = int(line[4])
				except:
					END = line[4]
				try:
					SCORE = int(line[5])
				except:
					SCORE = line[5]

				STRAND = line[6]
                
				try:
					PHASE = int(line[7])
				except:
					PHASE = line[7]
		    
				ATTRIBUTES = line[8]
 
				try:
					if eval(phrase):
						out.append(i)
				except:
					pass


		print str(len(out) - headline) + " lines written to output file.\n"
        
		self.save_file(filename_out, out)


	def makeDictionaryOfAttributes( self, i ):
		variables = {}
		for name_of_attribute in vars(self):
			if len( vars(self)[ name_of_attribute ] ) != 0:
				if name_of_attribute == "refURI" or name_of_attribute == "refLen":
					variables[ name_of_attribute ] = vars(self)[ name_of_attribute ][ self.seqID[i] ]
				else:
					variables[ name_of_attribute ] = vars(self)[ name_of_attribute ][i]
		return variables

	def checkWhichObjectsMatch( self, expression ):
		which_objects_dont_match = []
		length = max( [ len( list_of_attributes ) for list_of_attributes in vars(self).items() ] )
		for i in xrange( length ):
			variables = self.makeDictionaryOfAttributes( i )
			does_object_matches = eval(expression, variables) 
			if not does_object_matches:
				which_objects_dont_match.append(i)
		return which_objects_dont_match


	def filter(self, expression):
		which_objects_dont_match = self.checkWhichObjectsMatch( expression )
		which_objects_dont_match.sort( reverse = True )
		for list_of_attributes in vars( self ).values():
			if len( list_of_attributes ) > 0:
				for which_object in which_objects_dont_match:
					if type(list_of_attributes) == dict:
						which_object = self.seqID[ which_object ]
					list_of_attributes.pop( which_object )

			




