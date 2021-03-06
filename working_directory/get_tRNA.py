#!/usr/bin/env python

#module to get tRNA gene counts from GtRNAdb of species contained within CSandS

import BeautifulSoup as BS
import urllib2, os, sys, difflib

#given the CSandS abbreviation for a species, get the corresponding full name
def get_full_name(abrv_name,databases):
 	full_name = None
        for line in open(databases+'tRNA_abbreviations.csv'):
                if line.split(',')[0].strip() == abrv_name:
                        full_name = line.split(',')[1][:-1]
			break
	return full_name

#heuristicly match a full name to a url from the GtRNAdb
def get_url(full_name):

	best_match = ''
	best_score = 0
	best_url = ''

	#parse the species corresponding to each link
	for line in open(databases+'tRNA_links.csv'):
		species = line.split(',')[1].strip()

		#get the matching blocks between full name and species
		matches = difflib.SequenceMatcher(None,full_name,species).get_matching_blocks()#

		#score the match
		score = 0
		for match in matches:
			if match[2] > 2:
				#remove trailing blank space match
				if full_name[match[0]+match[2]-1] == ' ':
					score += match[2] - 1
				else:
					score += match[2]

		#test to see if best match
		if (score > best_score) or (score == best_score and len(species) < len(best_match)):
			best_match = species
			best_score = score
			best_url = line.split(',')[0].strip()

	return best_match,best_url,best_score


#download the gene count data from GtRNAdb and extract content
def get_data(url):
	
	#download html in packets
        f = urllib2.urlopen(url)
	html = ''
        while 1:
                packet = f.read()
                if not packet:
                        break
                html += packet
        f.close()

	#extract the four meanignful tables
	soup = BS.BeautifulSoup(html)
	four_box = soup('table')[10]
	six_box = soup('table')[11]
	two_box = soup('table')[12]
	other_box = soup('table')[13]

	#extract the tRNA gene count data from these tables
	data = []
	data.extend( parse_html_table( four_box ))
	data.extend( parse_html_table( six_box ))
	data.extend( parse_html_table( two_box ))
	data.extend( parse_html_table( other_box ))

	return data

#extract tRNA gene count data from these tables
def parse_html_table(html_table):
	data = []
	for row in html_table('tr')[1:]:
		amino = str(row.th.text)
		tmp = [str(entry.text) for entry in row.findAll('td')]
		codons = []
		for codon in tmp[:-1]:
			#check for a codon entry
			if len(codon.split(';')) != 3:
				#check to see if zero
				if codon[3:] == '&nbsp;':
					codons.append((codon[:3],0))
				else:
					codons.append((codon[:3],int(codon[3:])))
		data.append( (amino,codons) )
	return data

#save the gene counts
def save_data(fname,data):
	f = open(fname,'w')
	for amino in data:
		f.write('##'+amino[0]+'\n')
		for codon in amino[1]:
			f.write('#'+codon[0]+','+str(codon[1])+'\n')
	return

#read gene count data
def read_data(fname):
	f = open(fname,'r')
	data = []
	line = f.readline()
	while line:
		amino = line[2:-1]
		line = f.readline()
		codons = []
		while line and line[:2] != '##':
			line = line[1:-1].split(',')
			codons.append( (line[0],int(line[1])) )
			line = f.readline()
		data.append( (amino,codons) )
	return data

#get the tRNA gene counts
def get_tRNA(abrv_name,databases):

	#check to see if tRNA data present
	fname = databases+'tRNA_gene_counts/'+abrv_name+'.txt'

	#if yes: read it
	if os.path.exists(fname):
		data = read_data(fname)
	
	#if no:
	else:
		#make tRNA directory if it doesn't exist
		if not os.path.exists(databases+'tRNA_gene_counts/'):
    			os.makedirs(databases+'tRNA_gene_counts')

		#get the full name of species
		full_name = get_full_name(abrv_name,databases)

		#error check to see if full name was found
		if not full_name:
			print '\nERROR: No full name found for species',abbr_name
			print 'Please consider appending this species to the list of those excluded\n'
			sys.exit(1)

		#heuristicly match that name to a url from the GtRNAdb
		match,url,score = get_url(full_name)

		#error check to see if the match was good
		if (float(score) / min(len(match),len(full_name)) ) < 0.95:
			print '\nWARNING: Match for organism',abrv_name,'(',full_name,')','found to be below threshold'
			print 'Best macth was found to be',match
			print 'Please consider appending this species to the list of those excluded\n'

		#download the data and extract content
		data = get_data(url)

		#save the gene counts
		save_data(fname,data)

	return data
