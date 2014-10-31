import collections, itertools, math

#Alden Deran (adderan), BME 205 hw3

"""Functions for counting kmers in a sequence, iterating through kmers over an alphabet, checking the 
validity of kmers, and estimating the probability of a kmer by its number of occurrences in a file."""

def get_counts(parser, k, alphabet):
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
		

	for name, description, seq in parser:

		fixed_seq = fix_sequence(seq, alphabet, k)

		#Step through the sequence and increment the count for 
		#the substring of the sequence starting at the current position and 
		#continuing for k characters, where k is the length of the k-mers being counted.
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
		character_upper = character.upper()
		if character_upper in alphabet:
			good_characters.append(character_upper)

	full_sequence = start_chars + good_characters + stop_chars
	full_sequence = ''.join(full_sequence)
	return full_sequence

def iter_kmers(alphabet, k, include_start_kmer):
	"""Iterates through every valid kmer (string of length k) over the 
	provided alphabet. The validity of kmers is determined by the 
	is_valid_kmer() function. If include_start_kmer is True, the kmer containing
	only ^ characters will be yielded as one possible kmer. This may be undesirable
	when making conditional probability tables, because that kmer always occurs
	at the beginning of a sequence, but will be given a less-than-one probability based
	on counts, and will affect the counts of other kmers beginning with k-1 ^ characters. 
	The kmer consisting of only $ characters will always be yielded, because all characters
	after $ in a valid kmer must also be $. This means that, for example, the kmer $$$ will always have conditional
	probability 1 with condition $$, which is not true for the kmer ^^^ and condition ^^.
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
		if is_valid_kmer(kmer, include_start_kmer):
			yield ''.join(kmer)

def is_valid_kmer(kmer, include_start_kmer):
	"""Checks whether the provided kmer meets the following conditions:
	1) The '$' character does not appear before any character except '$'
	2) The '^' character does not appear after any character except '^'
	3) If include_start_kmer is False, the kmer does not consist of only ^ characters.
	If all the conditions are met, returns True. Otherwise returns false.
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
	#if mode is still zero, the kmer only has ^ characters
	if not include_start_kmer and mode == 0:
		valid_kmer = False
	return valid_kmer

	

def make_probability_dict(counts, alphabet, k):
	"""Returns a dict containing the base 2 logarithm of the estimated probabilities of
	each valid kmer over the provided alphabet, given it's condition. The alphabet should be a python 
	set. The condition of a kmer is defined as the string containing its first k-1 characters. The probabilities are estimated
	from a provided dictionary of counts of various kmers, obtained from a training dataset.
	A non-zero probability is guaranteed for every valid kmer over the given alphabet.
	"""
	#pseudocounts will store a count for every allowed kmer over the alphabet.
	#Valid kmers that aren't obesrved in the counts dict will be given a pseudocount of 1.0
	#so that the probability of every valid kmer is non-zero.
	pseudocounts = collections.defaultdict(float)

	#add the counts of all the observed kmers to the pseudocounts 
	for kmer, count in counts.items():
		if is_valid_kmer(kmer, False): #don't want to include start kmers here because they shouldn't affect other probabilities.
			pseudocounts[kmer] = float(count)

	#add 1.0 to the pseudocount of every valid kmer
	for kmer in iter_kmers(alphabet, k, False):
		pseudocounts[kmer] += 1.0
	
	#create a dictionary to hold the logarithm base 2 of the probabilities of each valid kmer.
	conditional_probabilities = collections.defaultdict(float)


	#Iterate through every valid condition. All valid conditions are generated by iterating through
	#each n-mer, where n = k - 1. Then each character from the alphabet can be appended to the condition
	#to yield every kmer with that condition (the resulting kmers are then checked for validity, because some
	#will be invalid). A sum of the pseudocounts of every kmer with that condition can be calculated, allowing the
	#conditonal probability of each kmer given its condition to be calculated by dividing its pseudocount by the sum
	#for its condition.

	#Should include start kmers here, because another character will be appended to the condition, and
	#kmers such as ^^A are valid.
	for condition in iter_kmers(alphabet, k - 1, True):
		sum = 0

		#Find the sum of pseudocounts over every kmer with this condition
		for character in alphabet:
			sum += pseudocounts[condition + character]
		
		#fill the conditional probability table for this condition
		for character in alphabet:
			kmer = condition + character

			#Since all characters in the alphabet are being tried, the validity of the kmer must be checked
			#because invalid kmers (and start-only kmers) will have zero pseudocount, causing a math error when taking the 
			#logarithm.
			if is_valid_kmer(kmer, False):
				conditional_probabilities[kmer] = math.log(pseudocounts[kmer], 2) - math.log(sum, 2)
	return conditional_probabilities
				
