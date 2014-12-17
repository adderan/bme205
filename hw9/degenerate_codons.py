#Alden Deran (adderan) bme 205 homework 9

from __future__ import print_function
import operator, itertools, string, collections
class DegenerateCodons:
	
	"""Dictionary that specifies the nucleotides
	that each degenerate base can represent in a degenerate codon."""
	degenerate_bases = {
		'A': ['A'],
		'G': ['G'],
		'C': ['C'],
		'T': ['T'],
		'U': ['T'],
		'R': ['G', 'A'],
		'Y': ['T', 'C'],
		'K': ['G', 'T'],
		'M': ['A', 'C'],
		'S': ['G', 'C'],
		'W': ['A', 'T'],
		'B': ['G', 'T', 'C'],
		'D': ['G', 'A', 'T'],
		'H': ['A', 'C', 'T'],
		'V': ['G', 'C', 'A'],
		'N': ['G', 'C', 'A', 'T']
	}

	def __init__(self, codon_table, codon_frequency):
		"""Creates an object for determining optimal degenerate codons. The codon_table contains
		each of the 64 codons and the amino acid it represents, and the frequency table contains the 
		frequency of usage for each codon (relative to other codons that represent that amino acid) for a 
		particular organism."""
		self.codon_table = codon_table
		self.codon_frequency = codon_frequency
		self.codon_expand_table = self.make_codon_expand_table()
		self.amino_counts = self.make_amino_counts()
		self.possible_dcodons = self.make_possible_dcodons()
	
	def generate_all_codons(self):
		"""Generates all possible degenerate codons using the standard
		degenerate base alphabet."""
		for codon in itertools.product(self.degenerate_bases.keys(), repeat = 3):
			yield ''.join(codon)
		
	def make_codon_expand_table(self):
		"""Makes a dictionary that maps each degenerate codon to a list of 
		codons that it can represent."""
		codon_expand_table = {}
		for codon in self.generate_all_codons():

			#create a list of the sets of bases that 
			#each character in the degenerate codon can translate
			#to. This will be used to create a cartesian product of 
			#these sets, which contains all codons the degenerate codon
			#can represent.
			degenerate_options = [self.degenerate_bases[base] for base in codon]
			codon_expand_table[codon] = []
			for expand_codon in itertools.product(*degenerate_options):
				expand_codon = ''.join(expand_codon)
				codon_expand_table[codon].append(expand_codon)
		return codon_expand_table
			
	def make_amino_counts(self):
		"""Makes a dictionary of the amino acids that each degenerate codon
		can code for. Each key in the dictionary is a degenerate codon, and each
		value is a dictionary of amino acids mapped to the number of times that acid is 
		coded for by the degenerate codon. An amino acid can be encoded by a degenerate codon
		more than once if the degenerate codon can represent two different codons that encode that
		amino acid."""
		amino_counts = {}
		for degenerate_codon in self.codon_expand_table.keys():
			amino_counts[degenerate_codon] = collections.Counter()

			#iterate through all codons that this codon can expand to
			#and check which amino acid that codon represents.
			for expand_codon in self.codon_expand_table[degenerate_codon]:
				acid = self.codon_table[expand_codon]
				amino_counts[degenerate_codon][acid] += 1
		return amino_counts
	def make_possible_dcodons(self):
		"""Makes a dictionary that translates sets of amino acids into lists
		of degenerate codons that can represent them. The amino acid sets are 
		represented by strings of the amino acid letters in alphabetical order,
		with stop codons represented by '*'. The values in the dictionary are
		lists of degenerate codons."""
		possible_dcodons = {}
		for degenerate_codon, multiset in self.amino_counts.items():
			
			#create a string of all the amino acids represented by a 
			#degenerate codon in alphabetical order. This string represents
			#a set of amio acids, and can be used as a key in a dictionary
			#to map sets of acids to lists of degenerate codons that represent
			#that set.
			keystring = ''.join(sorted(multiset.keys()))
			if keystring not in possible_dcodons:
				possible_dcodons[keystring] = []
			possible_dcodons[keystring].append(degenerate_codon)
		return possible_dcodons

	def imbalance(self, dcodon):
		"""Returns a score of the difference between the maximum and minimum
		multiplicity of the amino acids encoded by a degenerate codon. The codon
		might encode one amino acid in many different ways, and another in only one way,
		which would result in a high imbalance. The score is max(multiplicity) - min(multiplicity)/(total acids encoded)."""
		counts = self.amino_counts[dcodon].values()
		return (max(counts) - min(counts))/float(sum(counts))
		
	def dcodon_avg_frequency(self, dcodon):
		"""Finds the average codon usage frequency of a degenerate codon,
		by averaging the usage frequencies of all codons it represents."""
		possible_codons = self.codon_expand_table[dcodon]
		frequencies = [self.codon_frequency[c] for c in possible_codons]
		return float(sum(frequencies))/len(frequencies)

	def make_dcodon_info(self, dcodons):
		"""Takes a list of degenerate codons and creates a list of tuples (degenerate codon, imbalance, frequency)
		sorted by ascending imbalance, with ties broken by descending frequency."""
		codon_info = [(dc, self.imbalance(dc), self.dcodon_avg_frequency(dc)) for dc in dcodons]

		#sort by usage frequency descending, so that the highest frequency is first
		codon_info_freq = sorted(codon_info, key = operator.itemgetter(2), reverse = True)

		#sort by imbalance
		codon_info_imbalance = sorted(codon_info, key = operator.itemgetter(1))
		return codon_info_imbalance


	def print_minimal_table(self):
		"""Print a three-column table of sets of amino acids, the minimum imbalance for that set,
		and the optimal codon for that set, determined by the minimum imbalance and the maximum frequency,
		with imbalance taking higher priority. The rows are in alphabetical order by the amino acid sets."""
		for acid_set, dcodons in sorted(self.possible_dcodons.items(), key = operator.itemgetter(0)):
			dcodon_info = self.make_dcodon_info(dcodons)
			min_imbalance = dcodon_info[0][1]
			best_codon = dcodon_info[0][0]
			print("{:15s} {:.3f} {:10s}".format(acid_set, min_imbalance, best_codon))

	def print_min_codons_table(self):
		"""Print an alphabetical table of amino acid sets, the minimum imbalance for that set,
		and all the codons with minimum imbalance, sorted by usage frequency in the organism."""
		for acid_set, dcodons in sorted(self.possible_dcodons.items(), key = operator.itemgetter(0)):
			dcodon_info = self.make_dcodon_info(dcodons)
			min_imbalance = dcodon_info[0][1]
			best_codons = []
			for dcodon, imbalance, freq in dcodon_info:
				
				#only print the degenerate codons with minimum imbalance
				if imbalance == min_imbalance:
					best_codons.append(dcodon)
			best_codon_str = ','.join(best_codons)
			print("{:15s} {:.3f} {:s}".format(acid_set, min_imbalance, best_codon_str))
		
			
