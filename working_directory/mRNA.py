#!/usr/bin/env python

import re,sys

def get_mRNA(fold_data,speed_scores):

	rna_seq = []
	for datum in fold_data:

		#get rna for segment that was crystalised
		seg_start = int(datum.rna_aligned_start) - (int(datum.protein_aligned_start)-1)*3 - 1
		seg_end = seg_start + len(datum.protein_sequence)*3
		rna = datum.rna_sequence[seg_start:seg_end]

		for i,amino_acid in enumerate(datum.protein_sequence):

			#remove codons which have non-standard nucleotides
			if not rna[i*3:(i+1)*3] in speed_scores[datum.organism].keys():
				rna = rna[:i*3]+'---'+rna[(i+1)*3:]	
		
		#extract fold segments
		for i,pfold in enumerate(datum.protein_domain_fold):

			if pfold == fold:

				seg = datum.protein_domain_residue[i].split(':')[0] # only use first part of fold
				res_seg = map(int,re.findall('\d+', seg))
				if seg[0] == '-': #begins in a negative region
					res_seg[0] *= -1

				rna_start = (res_seg[0]-int(datum.protein_start))*3
				rna_end =  (res_seg[1]-int(datum.protein_start) + 1 )*3
				rna_seq.append(rna[rna_start:rna_end])

	return rna_seq


def align_mRNA(protein_seq,rna_seq):

	aligned_rna = []
	for i,seq in enumerate(protein_seq):

		count = 0
		alignment = ''
		for j,aa in enumerate(seq):

			#skip codons of amino acids not present in crystal structure
			if aa == '-':
				continue
		
			#iterate until you find the next amino acid in alignment
			while alignments[i][count] == '-':
				alignment += '---'
				count += 1
		
			#sanity check
			if alignments[i][count] != aa:
				print "Error: Mismatch when assigning RNA alignment"
				print fnames[i]
				print sum([1 for a in alignments[i] if a != '-'])
				sys.exit()

			alignment += rna_seq[i][j*3:(j+1)*3]
			count += 1
		
		alignment += (len(alignments[0])-len(alignment)/3) * '---'
		aligned_rna.append(alignment)

	return aligned_rna

