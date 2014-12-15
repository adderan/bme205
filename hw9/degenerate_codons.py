from __future__ import print_function
import operator, itertools, string, collections
class DegenerateCodons:
	
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
		self.codon_table = codon_table
		self.codon_frequency = codon_frequency
		self.codon_expand_table = self.make_codon_expand_table()
		self.amino_counts = self.make_amino_counts()
		self.possible_dcodons = self.make_possible_dcodons()
	
	def generate_all_codons(self):
		for codon in itertools.product(self.degenerate_bases.keys(), repeat = 3):
			yield ''.join(codon)
		
	def make_codon_expand_table(self):
		codon_expand_table = {}
		for codon in self.generate_all_codons():
			degenerate_options = [self.degenerate_bases[base] for base in codon]
			codon_expand_table[codon] = []
			for expand_codon in itertools.product(*degenerate_options):
				expand_codon = ''.join(expand_codon)
				codon_expand_table[codon].append(expand_codon)
		return codon_expand_table
			
	def make_amino_counts(self):
		amino_counts = {}
		for degenerate_codon in self.codon_expand_table.keys():
			amino_counts[degenerate_codon] = collections.Counter()
			for expand_codon in self.codon_expand_table[degenerate_codon]:
				acid = self.codon_table[expand_codon]
				amino_counts[degenerate_codon][acid] += 1
		return amino_counts
	def make_possible_dcodons(self):
		possible_dcodons = {}
		for degenerate_codon, multiset in self.amino_counts.items():
			keystring = ''.join(sorted(multiset.keys()))
			if keystring not in possible_dcodons:
				possible_dcodons[keystring] = []
			possible_dcodons[keystring].append(degenerate_codon)
		return possible_dcodons

	def imbalance(self, dcodon):
		counts = self.amino_counts[dcodon].values()
		return (max(counts) - min(counts))/float(sum(counts))
		
	def dcodon_avg_frequency(self, dcodon):
		possible_codons = self.codon_expand_table[dcodon]
		frequencies = [self.codon_frequency[c] for c in possible_codons]
		return float(sum(frequencies))/len(frequencies)

	def make_dcodon_info(self, dcodons):
		codon_info = [(dc, self.imbalance(dc), self.dcodon_avg_frequency(dc)) for dc in dcodons]
		codon_info = sorted(codon_info, key = operator.itemgetter(2), reverse = True)
		cdoon_info = sorted(codon_info, key = operator.itemgetter(1), reverse = True)
		return codon_info


	def print_minimal_table(self):
		for acid_set, dcodons in sorted(self.possible_dcodons.items(), key = operator.itemgetter(0)):
			dcodon_info = self.make_dcodon_info(dcodons)
			min_imbalance = dcodon_info[0][1]
			best_codon = dcodon_info[0][0]
			print("{:15s} {:.3f} {:10s}".format(acid_set, min_imbalance, best_codon))

	def print_min_codons_table(self):
		for acid_set, dcodons in sorted(self.possible_dcodons.items(), key = operator.itemgetter(0)):
			dcodon_info = self.make_dcodon_info(dcodons)
			min_imbalance = dcodon_info[0][1]
			best_codons = []
			for dcodon, imbalance, freq in dcodon_info:
				if imbalance == min_imbalance:
					best_codons.append(dcodon)
			best_codon_str = ','.join(best_codons)
			print("{:15s} {:.3f} {:s}".format(acid_set, min_imbalance, best_codon_str))
		
			
