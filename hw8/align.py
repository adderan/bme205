def score_a2m_global(s1, s2, subst, open_gap, extend_gap, double_gap, alphabet, alphabet_lower):
	s1_gap = False
	s2_gap = False
	score = 0
	s1_index = 0

	#no gap 0, s1 gap 1, s2 gap 2
	gap = 0

	for character in s2:

			if character in alphabet:
				score += subst[s1[s1_index] + character]
				s1_index += 1
				gap = 0
			elif character == "-":
				if gap == 0:
					score -= open_gap
					gap = 2
				elif gap == 1:
					score -= double_gap
					gap = 2
				elif gap == 2:
					score -= extend_gap
				s1_index += 1
			elif character in alphabet_lower:
				if gap == 0:
					score -= open_gap
					gap = 1
				elif gap == 1:
					score -= extend_gap
				elif gap == 2:
					score -= double_gap
					gap = 1
			else:
				pass
	return score

class local_aligner:
	def __init__(self, subst, open, extend, double):
		self.subst = subst
		self.open = open
		self.extend = extend
		self.double = double
		
		self.alphabet = set([])
		self.alphabet_lower = set([])
		for pair, score in self.subst.items():
			self.alphabet = self.alphabet.union(pair[0].upper())
			self.alphabet_lower = self.alphabet_lower.union(pair[0].lower())
	

	def score_a2m(self, s1, s2):
		alignment_columns = []
		for index in xrange(len(s2)):
			if s2[index] in self.alphabet:
				alignment_columns.append(index)
		alignment_start = min(alignment_columns)
		alignment_stop = max(alignment_columns)
		#global_score = score_a2m_global(s1, s2, self.subst, self.open, self.extend, self.double, self.alphabet, self.alphabet_lower)
		trimmed_s2 = s2[alignment_start:alignment_stop + 1]
		trimmed_s1 = s1[alignment_start:alignment_stop + 1]
		local_score = score_a2m_global(trimmed_s1, trimmed_s2, self.subst, self.open, self.extend, self.double, self.alphabet, self.alphabet_lower)
		return local_score

	def update_value(i, j, Ir, M, Ic):
		index = str(i) + "," + str(j)
		new_value = self.align_matrices[index]
		if Ir is not None:
			new_value[0] = Ir
		if M is not None:
			new_value[1] = M
		if Ic is not None:
			new_value[2] = Ic
		self.align_matrices[index] = new_value
	
	def get_value(i, j):
		index = str(i) + "," + str(j)
		return self.align_matrices[index]

	def	align(row_seq, col_seq):
		self.align_matrices = {}
		nrow = len(row_seq)
		ncol = len(col_seq)

		minus_inf = -100000


		m00 = self.subst(row_seq[0], col_seq[0])
		update_value(0, 0, minus_inf, m00, minus_inf)

		#Fill in first row and first column using the boundary conditions
		for i in xrange(1, nrow):
			j = 0
			M = self.subst(row_seq[i], col_seq[j])
			Ic = minus_inf

			prev = get_value(i-1, j)
			prevIr = prev[0]
			prevM = prev[1]
			prevIc = prev[2]
			Ir = max(prevM - self.open, prevIr - self.extend, prevIc - self.double)
			
			update_value(i, j, Ir, M, Ic)
		for j in xrange(1, ncol):
			i = 0
			M = self.subst(row_seq[i], col_seq[j])
			Ir = minus_inf
			prev = get_value(i, j - 1)
			prevIr = prev[0]
			prevM = prev[1]
			prevIc = prev[2]

			Ic = max(prevM - self.open, prevIr - self.double, prevIc - self.extend)
			
			update_value(i, j, Ir, M, Ic) 

			
		
		


class global_aligner:
	pass
