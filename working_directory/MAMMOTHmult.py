#!/usr/bin/env python

import re,get_pdb_fold,os,Bio,sys
from tempfile import NamedTemporaryFile

def get_pdb_domains(data,fold,databases):
	fnames = []
	protein_seq = []
	for datum in data:
		for i,pfold in enumerate(datum.protein_domain_fold):
			if pfold == fold:
				seg = datum.protein_domain_residue[i].split(':')[0] #only use first part of fold
				res_seg = map(int,re.findall('\d+', seg))
				if seg[0] == '-': #begins in a negative region
					res_seg[0] *= -1
				
				fname,seq = get_pdb_fold.get_pdb_fold(fold,datum.protein,datum.protein_chain,res_seg,databases) #must extract sequence due to errors in CSandS
				fnames.append( fname )
				protein_seq.append( seq )
	return fnames,protein_seq

##############################
#combine multiple pdb files
def combine_domains(fnames,fold,databases):
	fname = databases+'pdb_fold_cat/'+fold+'.pdb'
	#check whether the combined file already exists
	if not os.path.exists(fname):
		if not os.path.exists(databases+'pdb_fold_cat/'):
	    		os.makedirs(databases+'pdb_fold_cat')
		with open(databases+'pdb_fold_cat/'+fold+'.pdb', 'w') as outfile:
		    for fname in fnames:
			with open(fname) as infile:
			    outfile.write(infile.read())
	return fname

####################

#run MMAMMOTHmult on a file of pdb structures
def MAMMOTHmult(fname)

	#run MAMMOTHmult
	with NamedTemporaryFile() as f:
		    check_call(['mmult', fname], stdout=f, stderr=STDOUT)
		    f.seek(0)
		    output = f.read()

	#check for warnings
	if re.search('WARNING',output):
		print '\nWarning detected when performing structural alignment\n'
		print output
		print 'Discarding fold group:', fold
		print 'Folds removed:', len(fname)
		print 'Please add this fold to the excluded fold list\n'
		[os.remove(fold+'-FINAL.'+ext) for ext in ['aln','pdb','log'] ]
		return None,1

	#delete chaff
	[os.remove(fold+'-FINAL.'+ext) for ext in ['pdb','log'] ]
	
	#move alignment
	outname = databases+'structural_alignments/'+fold+'-FINAL.aln'
	os.rename(fold+'-FINAL.aln',outname)
	return outname,0

####################
	
#read the MAMMOTHmult alignment
def read_alignment(fname):
	#Get alignmnets
	with open(fname, 'r') as infile:
		alignments = list(Bio.SeqIO.parse(infile, "clustal"))
	return alignments

########################

def  sanity_checks(alignments,protein_seq):

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
			sys.exit()

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

	return seqs

###############
#main function to interface with MAMMOTHmult
def structural_alignment(data,fold,databases):

	fnames,seq = get_pdb_domains(data,fold,databases)
	
	fname = combine_domains(fnames,fold,databases)

	fname,flag = MAMMOTHmult(fname)
	if flag:
		return None,None,1

	alignments = read_alignment(fname)

	seq = sanity_check(alignmnets,seq)

	return alignments,seq,0

