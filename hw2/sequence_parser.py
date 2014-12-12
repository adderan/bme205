from __future__ import print_function
import string, re, sys, itertools
#Alden Deran (adderan) bme 205 homework 2

sep = re.compile(r"\s") #regular expression to match any whitespace character

def split_info_line(line):
	"""Splits a FASTA or FASTQ line containing the name and description
	of a sequence and returns the name and description as a tuple."""
	split = sep.search(line)
	if split: #check if whitespace was found
		name = line[1:split.start()] #set the name to the characters after "@" and before the whitespace
		description = line[split.end():]
	else:
		name = line[1:] #if no whitespace found, set the name to the whole line without the "@" character
		description = '' 
	return (name, description)

def read_fasta(infile, reading_quality_values = False):
	"""Reads either sequences or quality values from a FASTA file depending
	on whether reading_quality_values is True. Yields tuples of (name, description, sequence)
	for each sequence in the file."""
	name = None
	description = None
	seq = None

	sep = re.compile(r"\s") #a regular expression that will match any whitespace.
	for line in infile:
		line = line.rstrip()
		if line.startswith(">"): #check if the line is a name/description line
			if name is not None and description is not None and seq is not None:

				#a '>' character indicates a new sequence, so yield the previous sequence and then clear the variables.
				yield (name, description, seq) 				
				name = None
				description = None
				seq = None
			name, description = split_info_line(line)


		else:
			if seq is None: seq = ''
			seq += line
			if reading_quality_values: #add spaces between the quality values on different lines.
				seq += " "
	yield (name, description, seq) #yield the last sequence in the file.

def read_fastq(infile, phred_value):
	"""Reads a FASTQ file with the provided quality value encoding score.
	The quality characters will be converted to quality values by taking their 
	ASCII value and subtracting the phred_value. Tuples of (name, description, sequence, quality)
	will then be yielded, where quality is a list of integer quality values."""
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
					print("Warning: Names or descriptions do not match.", file = sys.stderr)
				mode = "quality"
			else:
				if seq is None: seq = ''
				seq += line
		elif mode == "quality":
			if quality is None: quality = []
			quality += list(line)
			if len(quality) == len(seq):
				mode = "nucleotide"
				decoded_quality = [ord(character) - phred_value for character in quality]
				yield (name, description, seq, decoded_quality)
				seq = None
				description = None
				name = None
				quality = None
def read_fasta_with_quality(sequencefile, qualityfile):
	"""Reads sequences from a FASTA sequence file and the corresponding quality
	strings from a FASTA quality file, and returns tuples of (name, description, sequence, quality)
	where quality is a list of integer quality values."""

	#Iterate through the sequences and their corresponding quality strings
	#by zipping together the yielded FASTA sequences from both files.
	for (name, description, sequence), (qualname, qualdescription, quality) in itertools.izip(read_fasta(sequencefile), read_fasta(qualityfile, reading_quality_values = True)):
		if qualname != name or qualdescription != description:
			print("Warning: Names or descriptions do not match.", sys.stderr)
		quality = quality.split()
		quality = [int(q) for q in quality] #convert the quality to integers
		yield (name, description, sequence, quality)

	
