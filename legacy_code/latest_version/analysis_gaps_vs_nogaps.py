#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from subprocess import check_call, STDOUT
from tempfile import NamedTemporaryFile
import Bio, get_pdb_fold, get_tAI, csands_folds_reader, os, sys, random, itertools, re

######################################################
#FUNCTIONS
######################################################

def smoothSpeed(speeds):
	count  = 0
	total = 0
	for speed in speeds:
		if speed:
			count += 1
			total += speed
	if count > (len(speeds)-1)/2:
		return total / count
	return None

def smoothAll(speeds, window):
	output = []
	for i in xrange(len(speeds)):
		if speeds[i]:
			j = max(i - (window - 1) / 2, 0)
			k = min(i + (window + 1) / 2, len(speeds))
			output.append(smoothSpeed(speeds[j:k]))
		else:
			output.append(None)
	return output

######################################################
#VARIABLES
######################################################

input_file = 'csands_folds_xray_40.txt'
window_size = 3
min_group = 7
alignment_switch = 0

print "Input file:",input_file
print "Window size:", window_size
print "Minimum group size:",min_group
print "Alignment switch:",alignment_switch

######################################################
#GET LIST OF EXCLUDED SPECIES (no tRNA match)
######################################################

excluded = []
for line in open('excluded_species.csv','r'):
	excluded.append(line.strip().split(',')[0])

######################################################
#GET DATA & COUNT FOLDS
######################################################

f = open(input_file,'r')
line = f.readline()
fold_count = dict()
org_count = dict()
speed_gap = dict()
speed_nogap = dict()
data = []

while line:
	datum = csands_folds_reader.read_entry(f,line[3:].strip())
	if datum.organism not in excluded:
		data.append(datum)
		for i,fold in enumerate(datum.protein_domain_fold):		
			if fold in fold_count.keys():
				fold_count[fold] += 1
			else:
				fold_count[fold] = 1
	line = f.readline()
f.close()
	
#sort folds so that the largest group is processed first
folds = sorted(fold_count, key=fold_count.get, reverse=True)
folds = [fold for fold in folds if fold_count[fold] >= min_group]

print "Total folds:",sum([fold_count[fold] for fold in folds])
print "Fold types:",len(folds),'\n'

for count,fold in enumerate(folds):

	print fold,':\t',count,"\tout of\t",len(folds)

######################################################
#GET SPECIFIC FOLD DATA
######################################################

	fold_data = []
	for datum in data:
		switch = 1
		for i,pfold in enumerate(datum.protein_domain_fold):
			if switch and pfold == fold:
				fold_data.append(datum)
				switch = 0

######################################################
#GET ORGANISM TYPES & SET UP DICTS
######################################################

	organisms = []
	for datum in fold_data:
		for i,pfold in enumerate(datum.protein_domain_fold):
			if pfold == fold:
				organisms.append(datum.organism)
				if datum.organism in org_count.keys():
					org_count[datum.organism] += 1
				else:
					org_count[datum.organism] = 1
					speed_gap[datum.organism] = []
					speed_nogap[datum.organism] = []

######################################################
#GET SPEED SCORES
######################################################

	speed_scores = dict()
	for org in organisms:
		if org not in speed_scores.keys():
			speed_scores[org] = get_tAI.get_tAI(org)
	
######################################################
#STRUCTURAL ALIGNMENT
######################################################

	#get pdb segments for each instance of fold
	fnames = []
	protein_seq = []
	for datum in fold_data:
		for i,pfold in enumerate(datum.protein_domain_fold):
			if pfold == fold:

				seg = datum.protein_domain_residue[i].split(':')[0] # only use first part of fold
				res_seg = map(int,re.findall('\d+', seg))
				if seg[0] == '-': #begins in a negative region
					res_seg[0] *= -1
				
				fname,seq = get_pdb_fold.get_pdb_fold(fold,datum.protein,datum.protein_chain,res_seg)
				fnames.append( fname )
				protein_seq.append( seq )
				
				
	#concatenate domains
	if not os.path.exists('pdb_cat/'):
    		os.makedirs('pdb_cat')
	with open('pdb_cat/'+fold+'.pdb', 'w') as outfile:
	    for fname in fnames:
		with open(fname) as infile:
		    outfile.write(infile.read())

	#run MAMMOTHmulti
	if alignment_switch or not os.path.isfile(fold+'-FINAL.aln'):
		with NamedTemporaryFile() as f:
		    check_call(['mmult',  'pdb_cat/'+fold+'.pdb'], stdout=f, stderr=STDOUT)
		    f.seek(0)
		    output = f.read()

		#check for warnings
		if re.search('WARNING',output):
			print '\nWarning detected when performing structural alignment\n'
			print output
			print 'Discarding fold group:', fold
			print 'Folds removed:', fold_count[fold],'\n'
			for org in organisms:
				org_count[org] -= 1
			[os.remove(fold+'-FINAL.'+ext) for ext in ['aln','pdb','log'] ]
			continue

		#delete chaff
		[os.remove(fold+'-FINAL.'+ext) for ext in ['pdb','log'] ]

	#Get alignmnets
	with open(fold+'-FINAL.aln', 'r') as infile:
		alignments = list(Bio.SeqIO.parse(infile, "clustal"))

######################################################
#Sanity Checks and cleaning up
######################################################

	#remove residues that failed to align
	seqs = []
	for i,a_seq,p_seq in zip(xrange(len(alignments)),alignments,protein_seq):
		tmp = [aa for aa in a_seq if aa != '-']
		tmp_2 = [aa for aa in p_seq if aa != '-']
		
		#same amount of residues ==> Good news
		if len(tmp) == len(tmp_2):
			seqs.append(p_seq)
			continue

		#more residues ==> very bad news
		if len(tmp) > len(tmp_2):
			print "\nError: Additional residues detected in alignment than inputed"
			print "Fold:",fold
			print "Protein",fnames[i]
			print "#",i
			print a_seq.seq,'\n'
			print p_seq
			sys.exit(1)

		#less residues ==> a wish for mmult to be documented
		print "\nWarning: Detected missing residues in alignment"
		print "Fold:",fold
		print "Protein",fnames[i]
		print "#",i,'\n'
		j = 0
		seq = ''
		for aa in tmp:
			while aa != p_seq[j]:
				j += 1
				seq += '-'
			seq += aa
			j += 1
		seq += (len(p_seq)-len(seq))*'-'
		seqs.append(seq)

	protein_seq = seqs

######################################################
#GET RNA SEGMENTS
######################################################

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

######################################################
#ALLIGN RNA SEGMENTS
######################################################

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
		
#######################################################
#CALCULATE SPEED PROFILE
#######################################################

	speed_profiles = []
	for rna_seq,org in zip(aligned_rna,organisms): 
		speed_profile = []
		for i in xrange(len(rna_seq)/3):
			codon = rna_seq[i*3:(i+1)*3]
			speed_profile.append(speed_scores[org][codon])
		
		speed_profiles.append(smoothAll(speed_profile,window_size))

########################################################
#Analysis
########################################################

	for i in xrange(len(speed_profiles[0])):
		gap_count = 0
		for j in xrange(len(speed_profiles)):
			if not speed_profiles[j][i]:
				gap_count +=  1
		if float(gap_count)/len(speed_profiles) >= 0.9:
			for j in xrange(len(speed_profiles)):
				speed_gap[organisms[j]].append( speed_profiles[j][i] )
		elif float(gap_count)/len(speed_profiles) <= 0.1:
			for j in xrange(len(speed_profiles)):
				speed_nogap[organisms[j]].append( speed_profiles[j][i] )
				
if not os.path.exists('plots/'):
	os.makedirs('plots')

for org in speed_gap.keys():
	if org_count[org] > 9:
		print org,org_count[org],len(speed_gap[org]),len(speed_nogap[org])
		plt.figure()
		n, bins, patches = plt.hist([speed_gap[org],speed_nogap[org]],25,range=(0,1),color=['0.75','crimson'],normed=1)
		plt.bar(bins[1:]-(bins[1]-bins[0])/2 - 1./50/2, n[0]-n[1] , 1./50)
		plt.legend(['With gaps (N='+str(len(speed_gap[org]))+')','No gaps (N='+str(len(speed_nogap[org]))+')','Difference'])
		plt.title('Distrubrution of speeds - '+org)
		plt.xlabel('tAI Score')
		plt.ylabel('Relative Frequency')
		plt.savefig('plots/speed_dist_'+org+'.png')
		#plt.show()
