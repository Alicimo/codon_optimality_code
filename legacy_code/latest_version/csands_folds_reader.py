#!/usr/bin/env python

class Datum_csands_folds:
	
	def __init__(self,ID,DBS,ORG,IDE,POS,COV,EVA,
			PRO,pCHA,DFD,DRS,STA,END,SSE,pSEQ,pASE,pALS,pALE,
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
		self.protein_domain_fold = DFD
		self.protein_domain_residue = DRS
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
		f.write('#DRS,'+','.join(self.protein_domain_residue)+'\n')
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

def read_entry(f,ID):

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
	DFD = read_csands(f) #domain folds
	DRS = read_csands(f) #domain residues
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

	return Datum_csands_folds(ID,DBS,ORG,IDE,POS,COV,EVA,
			PRO,pCHA,DFD,DRS,STA,END,SSE,pSEQ,pASE,pALS,pALE,
			RNA,nCHA,nSEQ,PSE,nASE,nALS,nALE)

def read_csands(f):
	return f.readline().strip().split(',')[1:]

if __name__ == "__main__":
	print "This is the reader class for interaction with the csands_folds database"
