import operator

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
	def __init__(self, subst, alphabet, alphabet_lower, open, extend, double):
		self.subst = subst
		self.open = open
		self.extend = extend
		self.double = double
		
		self.alphabet = alphabet
		self.alphabet_lower = alphabet_lower
			

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


	def update_value(self, i, j, Ir, M, Ic):
		index = str(i) + "," + str(j)
		old_value = (None, None, None)
		if index in self.align_matrices:
			old_value = self.align_matrices[index]
		IrNew = old_value[0]
		MNew = old_value[1]
		IcNew = old_value[2]
		if Ir is not None:
			IrNew = Ir
		if M is not None:
			MNew = M
		if Ic is not None:
			IcNew = Ic
		new_value = IrNew, MNew, IcNew
		self.align_matrices[index] = new_value
	
	def get_value(self, i, j):
		index = str(i) + "," + str(j)
		return self.align_matrices[index]
	
	def print_matrix(self):
		for i in range(self.row_seq_length):
			print([get_value(i, j)[0] for j in range(self.col_seq_length)])
				

	def	align(self, row_seq, col_seq):
		self.align_matrices = {}
		self.row_seq = row_seq
		self.col_seq = col_seq

		nrow = len(row_seq)
		ncol = len(col_seq)

		minus_inf = -100000


		m00 = self.subst[row_seq[0] + col_seq[0]]
		self.update_value(0, 0, minus_inf, m00, minus_inf)

		#Fill in first row and first column using the boundary conditions
		for i in xrange(1, nrow):
			j = 0
			M = self.subst[row_seq[i] + col_seq[j]]
			Ic = minus_inf

			upIr, upM, upIc = self.get_value(i-1, j)
			Ir = max(upM - self.open, upIr - self.extend, upIc - self.double)
			
			self.update_value(i, j, Ir, M, Ic)
		for j in xrange(1, ncol):
			i = 0
			M = self.subst[row_seq[i] + col_seq[j]]
			Ir = minus_inf
			leftIr, leftM, leftIc = self.get_value(i, j - 1)

			Ic = max(leftM - self.open, leftIr - self.double, leftIc - self.extend)
			
			self.update_value(i, j, Ir, M, Ic)

		#Fill in rest of matrix using recurrence relations
		for i in xrange(1, nrow):
			for j in xrange(1, ncol):
				upIr, upM, upIc = self.get_value(i-1, j)
				leftIr, leftM, leftIc = self.get_value(i, j-1)
				diagIr, diagM, diagIc = self.get_value(i-1, j-1)

				Sij = self.subst[row_seq[i] + col_seq[j]]

				Ir = max(upM - self.open, upIr - self.extend, upIc - self.double)
				Ic = max(leftM - self.open, leftIr - self.double, leftIc - self.extend)
				M = Sij + max(0, diagIc, diagIr, diagM)

				self.update_value(i, j, Ir, M, Ic)

		scores = [(cell, Ir, M, Ic) for cell, (Ir, M, IC) in self.align_matrices.items()]
		scores = sorted(scores, key = operator.itemgetter(2))
		return scores[0]

	def traceback_col_seq(self):
		res_seq = ""

		cell_scores = [(cell, Ir, M, Ic) for cell, (Ir, M, Ic) in self.align_matrices.items()]
		cell_scores = sorted(cell_scores, key = operator.itemgetter(2), reverse = True)
		start_cell = cell_scores[0][0]
		
		i = len(self.row_seq) - 1
		j = len(self.col_seq) - 1
		start_row, start_col = start_cell.split(sep = ",")

		res_seq += '-' * (len(self.row_seq) - start_row)
		i -= len(self.row_seq) - start_row
		assert(i == start_row)
		
		while j > start_col:
			res_seq += self.col_seq[j].lower()
			j -= 1
		assert(j == start_col)



		



class global_aligner:
	pass
