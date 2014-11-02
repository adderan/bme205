#Alden Deran (adderan) BME 205

import string, itertools, collections, math, sys

complement_table = string.maketrans("ACGT", "TGCA")
def reverse_comp(dna):
	return dna[::-1].translate(complement_table)

def iter_kmers(alphabet, k):
	alphabets = [alphabet for i in xrange(k)]

	for kmer in itertools.product(*alphabets):
		yield ''.join(kmer)

allowed_middle_characters = ['G', 'C']
def iter_palindromes(alphabet, k):
	for kmer in iter_kmers(alphabet, k/2):
		comp = reverse_comp(kmer)
		if k % 2 != 0:
			for character in allowed_middle_characters:
				yield kmer + character + comp
		else:
			yield kmer + comp

def iter_palindrome_range(alphabet, k1, k2):
	for k in xrange(k1, k2 + 1):
		for palindrome in iter_palindromes(alphabet, k):
			yield palindrome
			
		
def count_kmers(seq, k, counts):
	for start in xrange(len(seq) - k + 1):
		counts[seq[start:start + k]] += 1
	return counts

def remove_non_alphabetical(alphabet, seq):
	fixed_seq = []
	for character in seq:
		if character in alphabet:
			fixed_seq.append(character)
	return ''.join(fixed_seq)
			
def count_kmers_in_range(alphabet, seq, k1, k2, counts):
	fixed_seq = remove_non_alphabetical(alphabet, seq)
	for k in range(k1, k2 + 1):
		counts = count_kmers(fixed_seq, k, counts)
	return counts

def palindrome_statistics(palindrome, counts, N):
	prefix = palindrome[1:]
	suffix = palindrome[0:len(palindrome)-1]
	center = palindrome[1:len(palindrome)-1]

	#print(center, counts[center])
	if counts[center] == 0:
		return (palindrome, 0, 0, 0, 0, sys.maxint)
	else:
		expected_counts = float(counts[prefix])*float(counts[suffix])/counts[center]

	sigma = math.sqrt(expected_counts * (1 - expected_counts/N))

	Z = (counts[palindrome] - expected_counts)/sigma

	sign = 1
	if Z < 0:
		sign = -1
	P = math.erfc(abs(Z)/math.sqrt(2))/2
	return (palindrome, counts[palindrome], expected_counts, Z, P, N*P)


