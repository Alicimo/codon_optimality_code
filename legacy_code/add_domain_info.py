#!/usr/bin/env python

import re

###########################

class datum_csands:
	
	def __init__(self,ID,DBS,ORG,IDE,POS,COV,EVA,
			PRO,pCHA,DOM,STA,END,SSE,pSEQ,pASE,pALS,pALE,
			RNA,nCHA,nSEQ,PSE,nASE,nALS,nALE):

		self.identity = ID
		
		#Information
		self.databases = DBS
		self.organism_databases = ORG[1:]
		self.organism = ORG[0]
		self.indentical_hits = IDE[0]
		self.positive_hits = POS[0]
		self.coverage = COV[0]
		self.evalue = EVA[0]
		
		#protein information
		self.protein = PRO[0]
		self.protein_chain = pCHA[0]
		#self.protein_domain_class = DOM[1].split(':')
		#self.protein_domain_start = DOM[2].split(':')
		#self.protein_domain_end = DOM[3].split(':')
		self.protein_start = STA[0]
		self.protein_end = END[0]
		self.protein_secondary_struct = SSE[0]
		self.protein_sequence = pSEQ[0]
		self.protein_aligned_sequence = pASE[0]
		self.protein_aligned_start = pALS[0]
		self.protein_aligned_end = pALE[0]


		#RNA information
		self.rna = RNA[0]
		self.rna_chain = nCHA[0]
		self.rna_sequence = nSEQ[0]
		self.rna_protein_sequence = PSE[0]
		self.rna_aligned_sequence = nASE[0]
		self.rna_aligned_start = nALS[0]
		self.rna_aligned_end = nALE[0]

	def get_domains(self):
		f = open('dir.des.scop.1.75B-2013-06-20.txt')
		line = f.readline()
		tmp = []
		while line:
			if re.search('d'+self.protein.lower()+self.protein_chain.lower(), line):
				tmp.append(line.strip().split())
			line = f.readline()

		if len(tmp) == 0:
			self.protein_domain_fold = None
			self.protein_domain_residue = None
		else:
			if tmp[0][3][6] == '_':
				self.protein_domain_fold = [tmp[0][2]]
				self.protein_domain_residue = [[self.protein_start+'-'+self.protein_end]]
			else:
				folds = [entry[2] for entry in tmp]
				self.protein_domain_fold = folds

				residues = []
				for entry in tmp:
					res  = entry[5].split(',')
					res  =  [a[2:] for a in res]
					residues.append(res)
				self.protein_domain_residue = residues	
		f.close()
		


	def output(self,f):
		f.write('###'+self.identity+'\n')

		f.write('##INF\n')
		f.write('#DBS,'+','.join(self.databases)+'\n')
		f.write('#ORG,'+self.organism+','+','.join(self.organism_databases)+'\n')
		f.write('#IDE,'+self.indentical_hits+'\n')
		f.write('#POS,'+self.positive_hits+'\n')
		f.write('#COV,'+self.coverage+'\n')
		f.write('#EVA,'+self.evalue+'\n')

		f.write('##PRO,'+self.protein+'\n')
		f.write('#CHA,'+self.protein_chain+'\n')
		f.write('#DFD,'+','.join(self.protein_domain_fold)+'\n')

		tmp = [':'.join(residues) for residues in self.protein_domain_residue]
		f.write('#DRS,'+','.join(tmp)+'\n')

		f.write('#STA,'+self.protein_start+'\n')
		f.write('#END,'+self.protein_end+'\n')
		f.write('#SSE,'+self.protein_secondary_struct+'\n')
		f.write('#SEQ,'+self.protein_sequence+'\n')
		f.write('#ASE,'+self.protein_aligned_sequence+'\n')
		f.write('#ALS,'+self.protein_aligned_start+'\n')
		f.write('#ALE,'+self.protein_aligned_end+'\n')

		f.write('##RNA,'+self.rna+'\n')
		f.write('#CHA,'+self.rna_chain+'\n')
		f.write('#SEQ,'+self.rna_sequence+'\n')
		f.write('#PSE,'+self.rna_protein_sequence+'\n')
		f.write('#ASE,'+self.rna_aligned_sequence+'\n')
		f.write('#ALS,'+self.rna_aligned_start+'\n')
		f.write('#ALE,'+self.rna_aligned_end+'\n')
		


###########################

def read_entry_csands(f,ID):

	#INFORMATION SECTION
	f.readline()
	DBS = read_csands(f) #databases it's included
	ORG = read_csands(f) #organism (& databases)
	IDE = read_csands(f) #identical hits mRNA->protein (blastx)
	POS = read_csands(f) #positive hits mRNA->protein (blastx)
	COV = read_csands(f) #coverage of protein w.r.t mRNA
	EVA = read_csands(f) #e-value (presume w.r.t CSandS?)

	#AMINO ACID SECTION
	PRO = read_csands(f) #protein identifier
	pCHA = read_csands(f) #protein chain
	DOM = read_csands(f) #domains and domain identifiers
	STA = read_csands(f) #start residue
	END = read_csands(f) #end residue
	SSE = read_csands(f) #JOY assigned secondary struct
	pSEQ = read_csands(f) #protein seq
	pASE = read_csands(f) #aligned protein seq
	pALS = read_csands(f) #aligned start
	pALE = read_csands(f) #aligned end

	#NUCLEOTIDE SECTION
	RNA = read_csands(f) #RNA identifier
	nCHA = read_csands(f) #mRNA chain
	nSEQ = read_csands(f) #mRNA seq
	PSE = read_csands(f) #translation into protein
	nASE = read_csands(f) #aligned protein seq
	nALS = read_csands(f) #aligned start
	nALE = read_csands(f) #aligned end

	return datum_csands(ID,DBS,ORG,IDE,POS,COV,EVA,
			PRO,pCHA,DOM,STA,END,SSE,pSEQ,pASE,pALS,pALE,
			RNA,nCHA,nSEQ,PSE,nASE,nALS,nALE)


def read_csands(f):
	return f.readline().strip().split(',')[1:]


###########################
if __name__ == "__main__":

	f = open('csands.txt','r')
	g = open('csands_folds.txt','w')
	line = f.readline()
	count_1 = 0
	count_2 = 0
	while line:
		if re.match('^###\d+',line):
			count_1 += 1
			print "Entry:\t"+str(count_1)
			datum = read_entry_csands(f,line[3:].strip())
			datum.get_domains()
			if datum.protein_domain_fold != None:
				count_2 += 1
				print "Found fold"
				datum.output(g)
			else:
				print "No fold / Multidomain"
		line = f.readline()

	f.close()
	g.close()

	print "Entries found:\t"+str(count_1)
	print "Folds found:\t"+str(count_2)
