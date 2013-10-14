#!/usr/bin/env python

from string import maketrans
import get_tRNA,os

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

def codon2anti(codon):
	codon = codon[::-1]
	intab = "atcg"
	outtab = "tagc"
	trantab = maketrans(intab, outtab)
	return codon.translate(trantab)

def output_tAI(fname,scores):
	if not os.path.exists('tAI/'):
    		os.makedirs('tAI')
	with open(fname,'w') as outfile:
		for codon in scores.keys():
			outfile.write(codon+','+str(scores[codon])+'\n')
	
def read_tAI(fname):
	scores = dict()
	for line in open(fname,'r'):
		codon,score = line.strip().split(',')
		scores[codon] = float(score)
	return scores

def get_tAI(org):

	#stupid bloody names
	if org == 'ECALL':
		org = 'ECOLI'

	fname = 'tAI/'+org+'.txt'
	if os.path.exists(fname):
		speeds = read_tAI(fname)

	else:
		#get tRNA data
		gene_numbers = dict()
		tRNA_genes = get_tRNA.get_tRNA_data(org)
		for tRNA in tRNA_genes:
			for anti in tRNA[1]:
				gene_numbers[anti[0].lower()] = anti[1]

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
	

		#codons which do not generate scores are set to the arithmetic mean
		#Note that in the orginal text they we're set to the geometric mean 
		#but I see no reason as to why
		max_speed = max(speeds.values())
		mean_speed = 1.
		count = 0
		for codon in speeds.keys():
			speeds[codon] /= max_speed
			if speeds[codon] != 0:
				count += 1
				mean_speed += speeds[codon]
		mean_speed = mean_speed * (1.0/count)
		for codon in speeds.keys():
			if speeds[codon] == 0:
				speeds[codon] = mean_speed
		output_tAI(fname,speeds)

	speeds['---'] = None
	return speeds
