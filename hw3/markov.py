import collections, itertools, math

def get_counts(infile, parser, k, alphabet):
	"""Uses the provided parser to read sequences from infile, and counts the number of occurrences 
	of each k-mer in the file. The parser should yield tuples of (name, description, sequence). A k-mer 
	is defined as a string of length k, containing characters in the provided alphabet. Characters 
	not in the alphabet will be ignored when counting. K-mers can also
	contain the characters '^' and '$', which occur k times at the beginning and end of every sequence respectively,
	meaning that "^^A", "$$$", "^^^", and "CA$" are valid 3-mers in the sequence "ACA". The counts will be returned in a collections.Counter
	object, so that if "counts" is the name of the returned object, counts["AAA"] will be the number of occurrences 
	of "AAA" in the input file (assuming that k = 3).
	"""

	#A counter object to store the number of occurrences of each k-mer.
	#This means that counts["AAA"] is the number of occurences
	#of "AAA" in the input file.
	counts = collections.Counter()



	for name, description, seq in parser(infile):

		fixed_seq = fix_sequence(seq, alphabet, k) 

		#Step through the sequence and increment the count for 
		#the substring of the sequence starting at the current position and 
		#continuing for k characters, where k is the length of the k-mers being counted.
		#The one added to the range of the loop causes the last character to be counted.
		for kmer in read_kmer(fixed_seq, k):
			counts[kmer] += 1

	return counts

def read_kmer(seq, k):
	"""Reads kmers from seq and yields one at a time."""

	for start in xrange(len(seq) - k + 1):
		yield seq[start:start + k]

def fix_sequence(sequence, alphabet, k):
	"""Adds k start and stop characters (^ and $) to the beginning 
	and end of sequence. Deletes all characters in sequence that are not in alphabet.
	Converts all characters in sequence to uppercase.
	"""


	#the characters to be added to the beginning and end of every sequence,
	#allowing k-mers to go off the end of a sequence or to start before the beginning.
	stop_chars = ['$' for i in xrange(k)]
	start_chars = ['^' for i in xrange(k)]

	good_characters = []
	for character in sequence:
		if character in alphabet:
			good_characters.append(character.upper())
	good_characters = start_chars + good_characters + stop_chars
	fixed_sequence = ''.join(good_characters)
	return fixed_sequence

def iter_kmers(alphabet, k):
	"""Iterates through every valid kmer (string of length k) over the 
	provided alphabet. The validity of kmers is determined by the 
	is_valid_kmer() function.
	"""

	#make sure the start and stop characters are in the alphabet
	alphabet = alphabet.union(['$'])
	alphabet = alphabet.union(['^'])

	alphabet_list = [alphabet for i in xrange(k)]
	
	#Iterate through the cartesian product AxAx... where A is the alphabet
	#over which kmers are being listed, and the cross product occurs k times.
	#This will list every kmer over the alphabet.
	for kmer in itertools.product(*alphabet_list):
		
		#check if this kmer contains any $ or ^ characters in the wrong place, since
		#the cross product generates all kmers and some are invalid.
		if is_valid_kmer(kmer):
			yield ''.join(kmer)

def is_valid_kmer(kmer):
	"""Checks whether the provided kmer meets the following conditions:
	1) The '$' character does not appear before any character except '$'
	2) The '^' character does not appear after any character other than '^'
	If the conditions are met, returns True. Otherwise returns false.
	"""

	#In mode 0, ^ characters are expected, and other characters trigger a switch to the
	#appropriate mode. In mode 1, sequence characters are expected, and ^ characters indicate
	#an invalid kmer since the sequence has already started. In mode 2, $ characters are expected 
	#and any other character indicates an invalid kmer, since the sequence has ended.
	mode = 0
	valid_kmer = True
	for character in kmer:
		if mode == 0:
			if character == '$':
				mode = 2
			elif character != '^':
				mode = 1
		if mode == 1:
			if character == '$':
				mode = 2
			elif character == '^':
				valid_kmer = False
		if mode == 2:
			if character != '$':
				valid_kmer = False
	return valid_kmer

	

def make_probability_dict(counts, alphabet, order):
	"""Returns a dict containing the base 2 logarithm of the estimated probabilities of
	each valid kmer over the provided alphabet, given it's condition. The condition of a kmer
	is defined as the string containing its first k-1 characters. The probabilities are estimated
	from a provided dictionary of counts of various kmers, and a non-zero probability is guaranteed for
	every valid kmer over the given alphabet.
	"""

	k = order + 1

	pseudocounts = collections.defaultdict(float)
	for kmer, count in counts.items():
		pseudocounts[kmer] = float(count)

	for kmer in iter_kmers(alphabet, k):
		pseudocounts[kmer] += 1.0

	conditional_probabilities = collections.defaultdict(float)
	for kmer in iter_kmers(alphabet, k - 1):
		sum = 0
		for character in alphabet:
			sum += pseudocounts[kmer + character]
		for character in alphabet:
			sub_kmer = kmer + character
			if is_valid_kmer(sub_kmer):
				conditional_probabilities[sub_kmer] = math.log(float(pseudocounts[sub_kmer]), 2) - math.log(sum, 2)

	return conditional_probabilities
				
