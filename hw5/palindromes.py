#Alden Deran (adderan) BME 205

import string, itertools, collections, math, sys

complement_table = string.maketrans("ACGT", "TGCA")
def reverse_comp(dna):
	"""Returns the given sequence reversed, with each
	character replaced by its opposite strand 
	complement (A -> T, C -> G, T -> A, G -> C)."""
	
	#This function was directly copied from the BME 205 homework
	#5 assignment description.
	

	return dna[::-1].translate(complement_table)

def iter_kmers(alphabet, k):
	"""Generator function that yields every kmer (substring of length k) over an
	alphabet, which should be given as a Python set."""

	alphabets = [alphabet for i in xrange(k)]

	for kmer in itertools.product(*alphabets):
		yield ''.join(kmer)


def iter_palindromes(k, allowed_middle_characters, alphabet):
	"""Generator function that yields every DNA reverse-complement palindrome
	of length k, including odd palindromes with center characters determined by
	allowed_middle_characters."""

	for kmer in iter_kmers(alphabet, k/2):
		comp = reverse_comp(kmer)
		if k % 2 != 0:
			for character in allowed_middle_characters:
				yield kmer + character + comp
		else:
			yield kmer + comp

def iter_palindrome_range(k1, k2, allowed_middle_characters, alphabet):
	"""Generator function that yields all DNA reverse 
	complement palindromes from length k1 to k2, including
	palindromes of length k1 and k2."""

	for k in xrange(k1, k2 + 1):
		for palindrome in iter_palindromes(k, allowed_middle_characters, alphabet):
			yield palindrome
			
		
def count_kmers(seq, k):
	"""Counts all k-mers in seq and returns a
	collections.Counter() object containing
	the counts."""

	counts = collections.Counter()
	for start in xrange(len(seq) - k + 1):
		counts[seq[start:start + k]] += 1
	return counts

def remove_non_alphabetical(alphabet, seq):
	"""Removes all characters from seq that are not in the set alphabet."""
	fixed_seq = []
	for character in seq:
		if character in alphabet:
			fixed_seq.append(character)
	return ''.join(fixed_seq)
			
def count_kmers_in_range(seq, k1, k2):
	"""Returns a counter containing all kmers of length
	k1 through k2 in seq."""

	counts = collections.Counter()
	for k in range(k1, k2 + 1):
		counts = counts + count_kmers(seq, k)
	return counts


def find_expected_counts(palindrome, counts):
	"""Calculates the number of expected counts for a length-n palindrome
	by extending its center n-1 characters in both directions according to 
	a Markov chain. Takes a dictionary of counts for maximum-likelihood 
	estimation of the Markov chain conditional probabilities. Returns the
	number of expected counts"""
	suffix = palindrome[1:]
	prefix = palindrome[0:len(palindrome)-1]
	center = palindrome[1:len(palindrome)-1]
	if counts[center] == 0:
		print("Warning: Zero count for:", palindrome, sys.stderr)
		return 0
	else:
		expected_counts = float(counts[prefix])*float(counts[suffix])/counts[center]
		return expected_counts

def palindrome_statistics(palindrome, counts, N, N_hypotheses):
	"""Takes a palindrome as a string, kmer counts as a 
	collections.Counter(), and integers containing the number of
	characters that were read to measure the counts and the number
	of palindromes that are being studied in this experiment. Computes
	the expected number of occurrences of the palindrome in N characters
	assuming that palindromes are formed by a Markov Model in both directions.""" 
	
	expected_counts = find_expected_counts(palindrome, counts)
	observed_counts = counts[palindrome]
	if len(palindrome) %2 != 0:
		equiv_palindrome = reverse_comp(palindrome)
		expected_counts += find_expected_counts(equiv_palindrome, counts)
		observed_counts += counts[equiv_palindrome]
		
	sigma = math.sqrt(expected_counts * (1 - expected_counts/N))

	Z = (observed_counts - expected_counts)/sigma

	sign = 1
	if Z < 0:
		sign = -1
	P = math.erfc(abs(Z)/math.sqrt(2))/2
	return (palindrome, observed_counts, expected_counts, Z, P, N_hypotheses*P)


