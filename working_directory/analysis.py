#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from subprocess import check_call, STDOUT
from tempfile import NamedTemporaryFile
import Bio, get_pdb_fold, get_tAI, csands_folds_reader, os, sys, random, itertools, re, scipy, math

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

#CSandS database to use
input_file = 'csands_folds_xray_40.txt'

#location of databases
databases = '/local/databases/'

#window size to use for smoothing
window_size = 3

#minimum size of family
min_group = 7

#use to make the alignments be repeated or use those from a previous run
#==> if none found they will be automatically done
alignment_switch = 0

#thresholds for defining codons optimality
high_thresh = 0.9
low_thresh = 0.1

home = os.path.expanduser("~")
databases = home + databases
input_file = databases + input_file

print "\nInput file:",input_file
print "Database directory:",databases
print "Window size:", window_size
print "Minimum group size:",min_group
print "Alignment switch:",alignment_switch
print "High cut threshold:",high_thresh
print "Low cut threshold:",low_thresh

######################################################
#GET LIST OF EXCLUDED SPECIES (no tRNA match)
######################################################

excluded = []
for line in open(databases+'excluded_species_tRNA_40.csv','r'):
	excluded.append(line.strip().split(',')[0])

######################################################
#GET DATA & COUNT FOLDS
######################################################

f = open(input_file,'r')
line = f.readline()
fold_count = dict()
org_count = dict()
speed_scores = dict()
speed_obs = dict()
speed_high = dict()
speed_low = dict()
data = []

if not os.path.exists('plots/'):
	os.makedirs('plots')

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

######################################################
#speed dist
######################################################

for datum in data:
	tmp_scores = get_tAI.get_tAI(datum.organism,databases)
	tmp_speed = []
	
	speed_profile = []
	for i in xrange(len(datum.rna_sequence)/3):
			codon = datum.rna_sequence[i*3:(i+1)*3]
			try:
				speed_profile.append(tmp_scores[codon])
			except:
				continue
		
	if datum.organism in speed_obs.keys():
		speed_obs[datum.organism].extend(smoothAll(speed_profile,window_size))
	else:
		speed_obs[datum.organism] = smoothAll(speed_profile,window_size)

for org in speed_obs.keys():
	speed_obs[org].sort()
	speed_low[org] = speed_obs[org][int(math.floor(low_thresh*len(speed_obs[org])))]
	speed_high[org] = speed_obs[org][int(math.ceil(high_thresh*len(speed_obs[org])))]


######################################################
#fold loop (main)
######################################################

for iteration,fold in enumerate(folds):

	print fold,':\t',iteration+1,"\tout of\t",len(folds)

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
					speed_scores[datum.organism] = get_tAI.get_tAI(datum.organism,databases)			
	
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
				
				fname,seq = get_pdb_fold.get_pdb_fold(fold,datum.protein,datum.protein_chain,res_seg,databases)
				fnames.append( fname )
				protein_seq.append( seq )
				
				
	#concatenate domains
	if not os.path.exists(databases+'pdb_fold_cat/'):
    		os.makedirs(databases+'pdb_fold_cat')
	with open(databases+'pdb_fold_cat/'+fold+'.pdb', 'w') as outfile:
	    for fname in fnames:
		with open(fname) as infile:
		    outfile.write(infile.read())

	#run MAMMOTHmulti
	if alignment_switch or not os.path.isfile(databases+'structural_alignments/'+fold+'-FINAL.aln'):
		with NamedTemporaryFile() as f:
		    check_call(['mmult',  databases+'pdb_fold_cat/'+fold+'.pdb'], stdout=f, stderr=STDOUT)
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
		os.rename(fold+'-FINAL.aln',databases+'structural_alignments/'+fold+'-FINAL.aln')

	#Get alignmnets
	with open(databases+'structural_alignments/'+fold+'-FINAL.aln', 'r') as infile:
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

		speed_profiles.append( smoothAll(speed_profile,window_size))


########################################################
#Analysis
########################################################

	

	codon_classes = []	
	for speed,org in zip(speed_profiles,organisms):

		#generate randomized data and get distribution stats
#		copy_data = copy.deepcopy( speed_profile )
#		random_data = []
#		for i in xrange(10000):
#			random.shuffle(test_case)
#			random_data.extend( smoothAll(test_case,3) )
#		random_data = [x for x in random_data if x is not None]
#		mean = np.mean(np.array(data))
#		sigma = np.std(data)
#		low_limit = mean - 1.96*sigma
#		high_limit = mean + 1.96*sigma
#		n,bins,patches = plt.hist(data,25,range=(0,1),normed=1)
#		plt.show()
	
		codon_class = []
		for i in xrange(len(speed)):
			if not speed[i]:
				codon_class.append('-')
			elif speed[i] <= speed_low[org]:
				codon_class.append(-1)
			elif speed[i] >= speed_high[org]:
				codon_class.append(1)
			else:
				codon_class.append(0)
		codon_classes.append(codon_class)


	#calculate conservation score:
	conservation_score = []
	for i in xrange(len(codon_classes[-1])):
		n = 0
		count = 0
		for j in xrange(len(codon_classes)):
			if codon_classes[j][i] != '-':
					count += codon_classes[j][i]
					n += 1
		if n > 4:
			conservation_score.append( count / (1.*n) )	


	#randomize data
	duplicate_data =  [x for x in codon_classes]
	random_scores = []
	for i in xrange(10000):
		for j in xrange(len(duplicate_data)):
			random.shuffle(duplicate_data[j])
		for j in xrange(len(duplicate_data[0])):
			n = 0
			count = 0
			for k in xrange(len(duplicate_data)):
				if duplicate_data[k][j] != '-':
					count += duplicate_data[k][j]
					n += 1
			if n > 4:
				random_scores.append( count / (1. * n) )

	random_scores.sort()
	high_cut = random_scores[int(39*len(random_scores)/40.)]
	low_cut = random_scores[int(len(random_scores)/40.)]

	count = 0
	for score in conservation_score:
		if score <= low_cut or score >= high_cut:
			count += 1
	
	if iteration == 0:
		if count > 0.05*len(conservation_score):
			sig_conserved = 1
			print "Fold is under selection pressure"
		else:
			sig_conserved = 0
	else:
		if count > 0.05*len(conservation_score):
			sig_conserved += 1
			print "Fold is under selection pressure"


	plt.figure()
	n, bins, patches = plt.hist([random_scores,conservation_score],11,range=(-1,1),color=['0.75','crimson'],normed=1)
	plt.bar(bins[1:]-(bins[1]-bins[0])/2 - 1./25/2, n[0]-n[1] , 1./25)
	plt.legend(['Random','Observed (N='+str(len([a for a in conservation_score if a]))+')','Difference'],'upper left'  )
	plt.title('Distrubrution of scores - '+fold)
	plt.xlabel('Conservation Score')
	plt.ylabel('Relative Frequency')
	plt.axvline(x=low_cut,ls='--',color='k',lw=3)
	plt.axvline(x=high_cut,ls='--',color='k',lw=3)
	plt.savefig('plots/conservation_score_'+fold+'.png')

	if iteration == 0:
		cumulative = [n[0]-n[1]]
	else:
		cumulative.append(n[0]-n[1])

print "Folds under selection pressure: ",sig_conserved


import pickle
pickle.dump( cumulative, open( "plots/cumulative.p", "wb" ) )

totals = np.sum(cumulative,axis=0)
std = np.std(cumulative,axis=0)

plt.figure()
plt.bar(bins[1:]-(bins[1]-bins[0])/2 - 2./11/2, totals,
	width=2./11, yerr=std, color='0.75')
plt.title('Distrubrution of scores - Total')
plt.xlabel('Conservation Score')
plt.ylabel('Relative Frequency')
plt.savefig('plots/conservation_score_total.png')


		
		
