#!/usr/bin/env python

from string import maketrans
import os,get_tRNA

#parameter defining strength of interaction between codon and anticodon
def wobble_param(codon,anti):
	a = codon[2]
	b = anti[0]
	if a == codon2anti(b):
		return 0
	if a == 't':
		return 0.41
	if a == 'c':
		return 0.28
	if a == 'a':
		return 0.9999
	if a == 'g':
		return 0.68

#returns the complementary sequence for a given nucleotide sequence 
def codon2anti(codon):
	codon = codon[::-1]
	intab = "atcg"
	outtab = "tagc"
	trantab = maketrans(intab, outtab)
	return codon.translate(trantab)

#output tAI scores
def output_tAI(fname,scores,databases):
	if not os.path.exists(databases+'tAI/'):
    		os.makedirs(databases+'tAI')
	with open(fname,'w') as outfile:
		for codon in scores.keys():
			outfile.write(codon+','+str(scores[codon])+'\n')

#read tAI scores from file
def read_tAI(fname):
	scores = dict()
	for line in open(fname,'r'):
		codon,score = line.strip().split(',')
		scores[codon] = float(score)
	return scores

#get the tRNA gene copy numbers
def get_tRNA_counts(org,databases):
	gene_numbers = dict()
	tRNA_genes = get_tRNA.get_tRNA(org,databases)
	for amino,codons in tRNA_genes:
		for codon,count in codons:
			gene_numbers[codon.lower()] = count
	return gene_numbers

def calculate_tAI(gene_numbers):

	#establish speed scores for each codon
	speeds = dict()
	for first_aa in ['t','c','a','g']:
		for second_aa in ['t','c','a','g']:
			for third_aa in ['t','c','a','g']:
				codon = first_aa + second_aa + third_aa
				anti = codon2anti(codon)
				if third_aa == 't':
					speeds[codon] = gene_numbers[anti] + (1. - wobble_param(codon,'g'+anti[1:3])) * gene_numbers['g'+anti[1:3]]
				elif third_aa == 'c':
					speeds[codon] = gene_numbers[anti] + (1. - wobble_param(codon,'a'+anti[1:3])) * gene_numbers['a'+anti[1:3]]
				elif third_aa == 'a':
					speeds[codon] = gene_numbers[anti] + (1. - wobble_param(codon,'a'+anti[1:3])) * gene_numbers['a'+anti[1:3]]
				else:
					speeds[codon] = gene_numbers[anti] + (1. - wobble_param(codon,'t'+anti[1:3])) * gene_numbers['t'+anti[1:3]]

	#calculate the geometric mean
	total = 0
	count = 0
	for speed in speeds.values():
		if speed != 0:
			total *= speed
			count += 1
	geo_mean = total**(1./count)

	#codons which do not generate scores are set to the geometric mean
	for codon in speeds.keys():
		if speeds[codon] == 0:
			speeds[codon] = geo_mean

	return speeds


def get_tAI(org,databases):

	#stupid bloody names
	if org == 'ECALL':
		org = 'ECOLI'
	#check to see if the tAI scores have already been calculated
	fname = databases+'tAI_scores/'+org+'.txt'

	#if yes --> read them
	if os.path.exists(fname):
		speeds = read_tAI(fname)

	#if no --> get the tRNA data and calculate them
	else:
		#get tRNA data
		gene_numbers = get_tRNA_counts(org,databases)

		#calculate the tAI scores		
		speeds = calculate_tAI(gene_numbers)

		#save the tAI scores
		output_tAI(fname,speeds)
	
	#add the gap as a possible codon --> do not assign a speed
	speeds['---'] = None
	return speeds
