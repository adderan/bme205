from __future__ import print_function
import operator

#Alden Deran (adderan)

def score_a2m_global(s1, s2, subst, open_gap, extend_gap, double_gap, alphabet, alphabet_lower):
	"""Scores the global alignment of two sequences s1 and s2. S1 should be the master
	sequence from an a2m file, and s2 should be one of the other sequences. A substitution matrix and
	opening, extension, and double-gap penalties must be provided."""

	#scoring of local alignments can be done by trimming the beginning and ending
	#insertions and deletions and then scoring the remaining global alignment, which is why
	#this function is shared by both aligner classes

	score = 0
	s1_index = 0

	#state variable that is zero if no gap is currently open,
	#one if a gap is open in s1, and two if a gap is open in
	#s2.
	gap = 0

	for character in s2:
		if character in alphabet and character != "-" and character not in alphabet_lower:
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

def argmax(t):
	"""Returns the index of the maximum value in a tuple."""
	max_value = max(t)
	for i in range(len(t)):
		if t[i] == max_value: return i
	
class local_aligner:
	"""Contains functions for local alignment of two sequences using affine gap costs."""
	def __init__(self, subst, alphabet, open, extend, double):
		self.subst = subst
		self.open = open
		self.extend = extend
		self.double = double
		
		self.alphabet = alphabet

		#useful in scoring a2m files
		self.alphabet_lower = set([element.lower() for element in alphabet])
			

	def score_a2m(self, s1, s2):
		"""s1 and s2 should be the first and second sequences from
		an a2m-formatted file. The score of the local alignment represented
		by this file is calculated
		using the stored substitiution matrix and open, extend, and double
		gap penalties. Since this is a local alignment, deletions and insertions
		at the end of the sequence will not be included in the score."""

		start_deletions = 0
		start_insertions = 0
		end_deletions = 0
		end_insertions = 0
		alignment_started = False

		#count the insertions and deletions at the beginning
		#of the second sequence. Insertions in the
		#second sequence are represented in a2m by
		#lowercase characters, and deletions in the second sequence
		#(relative to the first) are represented by dashes.
		for character in s2:
			if character in self.alphabet:
				break
			elif character in self.alphabet_lower:
				start_insertions += 1
			elif character == "-":
				start_deletions += 1

		#count the insertions and deletions at the end
		#of s2
		for character in s2[::-1]:
			if character in self.alphabet:
				break
			elif character in self.alphabet_lower:
				end_insertions += 1
			elif character == "-":
				end_deletions += 1
		
		#trim s1 and s2 so that end deletions and end insertions aren't
		#included.
		trimmed_s2 = s2[start_insertions + start_deletions:len(s2) - end_insertions - end_deletions]
		trimmed_s1 = s1[start_deletions:len(s1) - end_deletions]

		#score the local alignment within s1 and s2 (not including
		#the end deletions and end insertions)
		local_score = score_a2m_global(trimmed_s1, trimmed_s2, self.subst, self.open, self.extend, self.double, self.alphabet, self.alphabet_lower)
		return local_score


	def update_value(self, i, j, Ir, M, Ic):
		"""Updates the stored value of the tuple (Ir, M, Ic) at the cell i,j in
		the alignment matrices. The cell i,j represents either an insertion, deletion, 
		or match at the i-th character from the row_seq and the j-th character from the col_seq.
		Each of the values are only updated if the given replacement value
		is not None."""


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
		"""Gets the value of all three alignment matrices (Ir, M, IC)
		as a tuple, for the cell i,j."""

		index = str(i) + "," + str(j)
		return self.align_matrices[index]		

	def	align(self, row_seq, col_seq):
		"""Fllls in the local alignment matrices (Ir, M, Ic) for row_seq and 
		col_seq using affine gap costs."""
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
		scores = sorted(scores, key = operator.itemgetter(2), reverse = True)
		return scores[0][2]
	
	
	def traceback_col_seq(self):
		"""Uses the prcomputed alignment matrices to discover one of the optimal
		local alignments for row_seq and col_seq."""
		res_seq = ""
		
		#find the score of the optimal local alignment ending with a match between two characters.
		cell_scores = [(cell, Ir, M, Ic) for cell, (Ir, M, Ic) in self.align_matrices.items()]
		cell_scores = sorted(cell_scores, key = operator.itemgetter(2), reverse = True)

		#cell coordinates of the end of the local alignment
		start_cell = cell_scores[0][0]

		#If the best score is less than zero,
		#an empty alignment is optimal
		max_score = cell_scores[0][2]
		if max_score <= 0:
			empty_alignment = "-"*len(self.row_seq) + self.col_seq.lower()
			return empty_alignment

		
		#start at the bottom-right cell in the alignment matrix, and
		#print out dashes and lowercase characters for each sequence
		#until the start of the alignment is reached. These insertions
		#and deletions will not contribute to the score of the local alignment.
		i = len(self.row_seq) - 1
		j = len(self.col_seq) - 1

		start_row, start_col = start_cell.split(",")
		
		start_row = int(start_row)
		start_col = int(start_col)

		while i > start_row:
			res_seq += '-'
			i = i - 1
		
		while j > start_col:
			res_seq += self.col_seq[j].lower()
			j -= 1

		#State representing which subproblem is currently 
		#being evaluated. The alignment starts with a match, 
		#so the subproblem is initially one. If a gap is opened in row_seq,
		#the state will transition to 0, or to 2 if a gap is opened in
		#col_seq
		subproblem = 1
		
		#the penalties are different depending on whether a gap 
		#is already open, is not open, or is open on the opposite sequence
		penaltiesIr = (self.extend, self.open, self.double)
		penaltiesM = (0, 0, 0)
		penaltiesIc = (self.double, self.open, self.extend)

		#traverse the matrix until reaching the first row or first column.
		#At that point, the alignment has terminated. It can also terminate 
		#within the matrix (on a match) if all options for continuing provide
		#a negative score.
		while i > 0 and j > 0:
			#output the character for the current square based
			#on the subproblem.
			if subproblem == 1:
				res_seq += self.col_seq[j].upper()
				j = j - 1
				i = i - 1
			elif subproblem == 0:
				res_seq += "-"
				i = i - 1
			elif subproblem == 2:
				res_seq += self.col_seq[j].lower()
				j = j - 1

			#decide which section of the next square to go to
			#next
			nextIr, nextM, nextIc = self.get_value(i, j)
			penalties = None

			if subproblem == 0: penalties = penaltiesIr
			elif subproblem == 1: penalties = penaltiesM
			elif subproblem == 2: penalties = penaltiesIc
			adjustedNextScores = (nextIr - penalties[0], nextM - penalties[1], nextIc - penalties[2])
			
			#if at a match, and all the scores for the next cell are negative,
			#end the alignment.
			if subproblem == 1 and adjustedNextScores[0] < 0 and adjustedNextScores[1] < 0 and adjustedNextScores[2] < 0:
				break

			#assign the next state to Ir, M, or Ic depending on which
			#has the maximum score minus penalty for the next cell.
			subproblem = argmax(adjustedNextScores)

		#add the last character if its score is greater than
		#zero. 
		if self.subst[self.row_seq[i] + self.col_seq[j]] > 0:
			res_seq += self.col_seq[j].upper()
			i = i - 1
			j = j - 1

		#print beginning gaps (deletions) and lowercase characters (insertions)
		while j >= 0:
			res_seq += self.col_seq[j].lower()
			j = j - 1
		while i >= 0:
			res_seq += "-"
			i = i - 1

		res_seq_forward = res_seq[::-1]
		#print("score should be: ", cell_scores[0][2])
		#print("score is: ", self.score_a2m(self.row_seq, res_seq_forward))
		return res_seq_forward


class global_aligner:
	"""Class containing functions for global alignment of sequences
	with affine gap costs. The align() function fills in the Ir, M, and Ic 
	matrices, which are then used by the traceback_col_seq() function to recover
	one of the optimal global alignments."""

	def __init__(self, subst, alphabet, open_penalty, extend_penalty, double_gap_penalty):
		"""Creates a new aligner object with the specified gap opening, gap extension,
		and double gap penalties. Stores the alphabet of possible characters and a substitution
		matrix that gives scores for pairing each character with each other character in an alignment."""
		self.subst = subst
		self.alphabet = alphabet

		#also store the lowercase alphabet, which is used for recognizing
		#insertion characters in a2m-formatted sequences.
		self.alphabet_lower = set([element.lower() for element in alphabet])

		self.open = open_penalty
		self.extend = extend_penalty
		self.double = double_gap_penalty

		#a dict 
		self.align_matrices = {}
	
	def score_a2m(self, s1, s2):
		"""Scores two sequences in a2m format using affine gap costs. s1 should be the master sequence
		and s2 a sequence with possible insertions, deletions, or substitutions relative to s1."""
		return score_a2m_global(s1, s2, self.subst, self.open, self.extend, self.double, self.alphabet, self.alphabet_lower)

	def update_value(self, i, j, Ir, M, Ic):
		"""Updates the stored value of the tuple (Ir, M, Ic) at the cell i,j in
		the alignment matrices. The cell i,j represents either an insertion, deletion, 
		or match at the i-th character from the row_seq and the j-th character from the col_seq.
		Each of the values are only updated if the given replacement value
		is not None."""

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
		"""Gets the value of all three alignment matrices (Ir, M, IC)
		as a tuple, for the cell i,j."""

		index = str(i) + "," + str(j)
		return self.align_matrices[index]		

	def	align(self, row_seq, col_seq):
		"""Fills in the Ir, M, and Ic matrices for global alignment  of
		row_seq and col_seq with affine gap costs."""

		self.row_seq = row_seq
		self.col_seq = col_seq

		nrow = len(row_seq)
		ncol = len(col_seq)

		minus_inf = -100000
		
		#fill in the start-start cell with 0 in the match state,
		#minus infinity in the other states. The alignment will always start
		#in this cell, so matching the start character with any other character
		#is prohibited by the minus infinite score.
		self.update_value(-1, -1, minus_inf, 0, minus_inf)

		#fill in the (0, -1) and (-1, 0) cells that correspond to
		#a gap opening at the beginning of the sequence. These gaps are 
		#charged the full gap opening penalty
		self.update_value(-1, 0, minus_inf, minus_inf, -1*self.open)
		self.update_value(0, -1, -1*self.open, minus_inf, minus_inf)


		#Fill in the start row and start column using the boundary conditions. For the
		#start row, Ir is always minus infinity because it looks one cell up, which is undefined.
		#Ic decreases linearly at each cell because a gap has been opened and is being extended.
		for i in xrange(1, nrow):
			j = -1
			Ic = minus_inf
			M = minus_inf

			upIr, upM, upIc = self.get_value(i-1, j)
			Ir = upIr - self.extend

			self.update_value(i, j, Ir, M, Ic)
		for j in xrange(1, ncol):
			i = -1
			Ir = minus_inf
			M = minus_inf
			leftIr, leftM, leftIc = self.get_value(i, j - 1)

			Ic = leftIc - self.extend			
			self.update_value(i, j, Ir, M, Ic)

		#Fill in the rest of the matrix using the recurrence relations for global alignment.
		for i in xrange(0, nrow):
			for j in xrange(0, ncol):
				upIr, upM, upIc = self.get_value(i-1, j)
				leftIr, leftM, leftIc = self.get_value(i, j-1)
				diagIr, diagM, diagIc = self.get_value(i-1, j-1)

				Sij = self.subst[row_seq[i] + col_seq[j]]

				Ir = max(upM - self.open, upIr - self.extend, upIc - self.double)
				Ic = max(leftM - self.open, leftIr - self.double, leftIc - self.extend)

				#this is different from local alignment because the alignment
				#can't start anywhere but the start character. There is no option
				#to stop at i,j even if all the next cells are negative.
				M = Sij + max(diagIc, diagIr, diagM)

				self.update_value(i, j, Ir, M, Ic)

		final_scores = self.get_value(nrow-1, ncol-1)

		return max(final_scores)


	def print_matrix(self):
		"""Prints the Ir, M, and Ic matrices for debugging."""

		for i in range(-1, len(self.row_seq)):
			for j in range(-1, len(self.col_seq)):
				print(self.get_value(i,j)[0], " ", end = '')
			print("")
		print("")
		for i in range(-1, len(self.row_seq)):
			for j in range(-1, len(self.col_seq)):
				print(self.get_value(i,j)[1], " ", end = '')
			print("")
		print("")
		for i in range(-1, len(self.row_seq)):
			for j in range(-1, len(self.col_seq)):
				print(self.get_value(i,j)[2], " ", end = '')
			print("")


	def traceback_col_seq(self):
		"""Uses the pre-computed alignment matrices to find one 
		of the optimal global alignments of row_seq and col_seq."""

		res_seq = ""
		
		final_scores = self.get_value(len(self.row_seq) -1, len(self.col_seq) -1)

		#Start at the maximum score from the cell in the bottom right corner, 
		#corresponding to a match, insertion, or deletion at
		#the last character of both sequences. One must start at this 
		#cell since this is a global alignment and every character 
		#has to be aligned to either a gap or a character from the 
		#opposite sequence.
		i = len(self.row_seq) - 1
		j = len(self.col_seq) - 1

		subproblem = argmax(final_scores)#State representing which subproblem is currently 
		#being evaluated. The alignment starts with a match, 
		#so the subproblem is initially one. If a gap is opened in row_seq,
		#the state will transition to 0, or to 2 if a gap is opened in
		#col_seq

		penaltiesIr = (self.extend, self.open, self.double)
		penaltiesM = (0, 0, 0)
		penaltiesIc = (self.double, self.open, self.extend)

		while True:
			#output the character for the current cell based
			#on the subproblem.
			if subproblem == 1:
				res_seq += self.col_seq[j].upper()

				#since this was a match, go up one cell and left one cell,
				#as determined by the recurrence relations.
				j = j - 1
				i = i - 1
			elif subproblem == 0:

				#since this was a gap in col_seq, go up
				#one row but stay at the same column, because
				#the current character from col_seq still must
				#be added.
				res_seq += "-"
				i = i - 1
			elif subproblem == 2:
				res_seq += self.col_seq[j].lower()
				j = j - 1
			if j < 0 and i < 0: break
			#decide which section of the next square to go to
			#next
			nextIr, nextM, nextIc = self.get_value(i, j)
			penalties = None

			if subproblem == 0: penalties = penaltiesIr
			elif subproblem == 1: penalties = penaltiesM
			elif subproblem == 2: penalties = penaltiesIc
			adjustedNextScores = (nextIr - penalties[0], nextM - penalties[1], nextIc - penalties[2])
			

			#assign the next state to Ir, M, or Ic depending on which
			#has the maximum score minus penalty for the next cell.
			subproblem = argmax(adjustedNextScores)


		res_seq_forward = res_seq[::-1]
		return res_seq_forward


