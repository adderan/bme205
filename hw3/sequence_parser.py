from __future__ import print_function
import string, re, sys, itertools


sep = re.compile(r"\s") #regular expression to match any whitespace character

def split_info_line(line):
	
	split = sep.search(line)
	if split: #check if whitespace was found
		name = line[1:split.start()] #set the name to the characters after "@" and before the whitespace
		description = line[split.end():]
	else:
		name = line[1:] #if no whitespace found, set the name to the whole line without the "@" character
		description = '' 
	return (name, description)

def read_fasta(infile, reading_quality_values = False):
	name = None
	description = None
	seq = None

	sep = re.compile(r"\s") #a regular expression that will match any whitespace.
	for line in infile:
		line = line.rstrip()
		if line.startswith(">"): #check if the line is a name/description line
			if name is not None and description is not None and seq is not None:
				yield (name, description, seq) #a '>' character indicates a new sequence, so yield the previous sequence and then clear the variables.
				name = None
				description = None
				seq = None
			name, description = split_info_line(line)


		else: #if this is a sequence line, append the line to the current sequence, first initializing the current sequence to '' if it is None.
			if seq is None: seq = ''
			seq += line
			if reading_quality_values:
				seq += " "
	yield (name, description, seq) #yield the last sequence in the file.

def read_fastq(infile):
	name = None
	description = None
	seq = None
	quality = None

	mode = "nucleotide"
	#The mode will either be "nucleotide" or "quality". In nucleotide mode, encountering a "+" character at the
	#beginning of a line will cause a switch to quality mode. This works because the "+" character cannot appear in nucleotide
	#sequences. In quality mode, the (name, description, seq, quality)
	#tuple will be yielded when the quality string reaches exactly the same length as the nucleotide string, and the mode will
	#switch back to nucleotide. 
	
	for line in infile:
		line = line.rstrip()

		if mode == "nucleotide":
			if line.startswith("@"):
				name, description = split_info_line(line)
			elif line.startswith("+"):
				qualname, qualdescription = split_info_line(line)
				if qualname is not None and (qualname != name or qualdescription != description):
					print("Warning: Names or descriptions do not match.", sys.stderr)
				mode = "quality"
			else:
				if seq is None: seq = ''
				seq += line
		elif mode == "quality":
			if quality is None: quality = ''
			quality += line
			if len(quality) == len(seq):
				mode = "nucleotide"
				yield (name, description, seq, quality)
				seq = None
				description = None
				name = None
				quality = None
def read_fasta_with_quality(sequencefile, qualityfile):
	for (name, description, sequence), (qualname, qualdescription, quality) in itertools.izip(read_fasta(sequencefile), read_fasta(qualityfile, reading_quality_values = True)):
		if qualname != name or qualdescription != description:
			print("Warning: Names or descriptions do not match.", sys.stderr)
		quality = quality.split()
		yield (name, description, sequence, quality)

	

	
