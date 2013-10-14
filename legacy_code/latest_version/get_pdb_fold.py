#!/usr/bin/env python
#
# Downloads pdb files and extracts the specified residues (fold region) from a specific 
# chain. Returns the path to the fold file.

import Bio, Bio.PDB, get_pdb, os, sys

class Selector(Bio.PDB.Select):
 
	def __init__(self, chain_letter, res_ranges):
		self.chain_letter = chain_letter
		self.res_ranges = res_ranges

	def accept_model(self, model):
		return (model.get_id() == 0)

	def accept_chain(self, chain):
		return (chain.get_id() in self.chain_letter)

	def accept_residue(self, residue):
		#return any( start <= residue.get_id()[1] <= end for start,end in self.res_ranges )
		if residue.get_id()[1] in range(self.res_ranges[0],self.res_ranges[1]+1):
			if residue.get_id()[0] == ' ' and residue.get_id()[2] == ' ' and residue.has_id('CA'):
				return  True
		return False	

def get_pdb_seq(pdb_id,fname,res_range):

	#get pdb info
	s = Bio.PDB.PDBParser().get_structure(pdb_id,fname)
	ids = [res.get_id()[1] for res in s.get_residues()]
	seq = [Bio.PDB.Polypeptide.three_to_one(res.get_resname()) for res in s.get_residues()]
	#seq = [pp.get_sequence() for pp in Bio.PDB.PPBuilder().build_peptides(s)]
	#seq = str(reduce(lambda x, y: x+y, seq))

	#sanity check
	if len(ids) != len(seq):
		print "Error 1: mismatch between ids and seq",pdb_id
		print "File:",fname 
		print "Length of ids",len(ids)
		print "Length of seq",len(seq)
		print ids,'\n'
		print seq
		sys.exit()

	#map onto domain range given:
	seg = ''
	count = 0
	for i in xrange(res_range[0],res_range[1]+1):
		if count != len(seq) and ids[count] == i :
			seg += seq[count]
			count += 1
		else:
			seg += '-'
	return seg


def get_pdb_fold(fold,pdb_id,chain,res_ranges):
	fname = 'pdb_fold/'+fold+'/'+pdb_id+'_'+chain+'_'+str(res_ranges[0])+'_'+str(res_ranges[1])+'.pdb'

	#check if fold is present
	if os.path.exists(fname):
		seq = get_pdb_seq(pdb_id,fname,res_ranges)
		return fname,seq

	#check for general fold dir
	if not os.path.exists('pdb_fold/'):
    		os.makedirs('pdb_fold')

	#check for specific fold dir
	if not os.path.exists('pdb_fold/'+fold+'/'):
		os.makedirs('pdb_fold/'+fold+'/')

	gname = get_pdb.get_pdb(pdb_id)
	io=Bio.PDB.PDBIO()
	s=Bio.PDB.PDBParser().get_structure(pdb_id,gname)
	io.set_structure(s)
	io.save(fname, select=Selector(chain, res_ranges))
	seq = get_pdb_seq(pdb_id,fname,res_ranges)

	return fname,seq
