#!/usr/bin/env python

import BeautifulSoup as BS
import urllib2, re, os, sys, difflib

def get_full_name(abrv_name):

 	full_name = None
        for line in open('tRNA_species.csv'):
                if line.split(',')[0].strip() == abrv_name:
                        full_name = line.split(',')[1][:-1]
			break
	return full_name

def get_url(full_name):
	
	best_match = ''
	best_score = 0
	best_url = ''
	for line in open('tRNA_links.csv'):
		species = line.split(',')[1].strip()
		matches = difflib.SequenceMatcher(None,full_name,species).get_matching_blocks()
		score = 0
		for match in matches:
			if match[2] > 2:
				if full_name[match[0]+match[2]-1] == ' ':
					score += match[2] - 1
				else:
					score += match[2]
		if (score > best_score) or (score == best_score and len(species) < len(best_match)):
			best_match = species
			best_score = score
			best_url = line.split(',')[0].strip()
	return best_match,best_url,best_score


def get_data(url):
	
        f = urllib2.urlopen(url)
	html = ''
        while 1:
                packet = f.read()
                if not packet:
                        break
                html += packet
        f.close()

	soup = BS.BeautifulSoup(html)
	
	four_box = soup('table')[10]
	six_box = soup('table')[11]
	two_box = soup('table')[12]
	other_box = soup('table')[13]

	data = []
	data.extend( parse_html_table( four_box ))
	data.extend( parse_html_table( six_box ))
	data.extend( parse_html_table( two_box ))
	data.extend( parse_html_table( other_box ))

	return data


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

def save_data(fname,data):
	f = open(fname,'w')
	for amino in data:
		f.write('##'+amino[0]+'\n')
		for codon in amino[1]:
			f.write('#'+codon[0]+','+str(codon[1])+'\n')
	return

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

def get_tRNA_data(abrv_name):

	#check to see if data present
	fname = 'tRNA/'+abrv_name+'.txt'
	if not os.path.exists(fname):
		#make tRNA directory
		if not os.path.exists('tRNA/'):
    			os.makedirs('tRNA')
		full_name = get_full_name(abrv_name)
		match,url,score = get_url(full_name)
		data = get_data(url)
		save_data(fname,data)
		if (float(score) / min(len(match),len(full_name)) ) < 0.95:
			print '\nMatch for organism',abrv_name,'(',full_name,')','found to be below threshold'
			print 'Best macth was found to be',match
			print 'Please consider appending this species to the list of those excluded\n'
			#print ','.join([abrv_name,full_name,match,url,str(min(len(match),len(full_name))),str(score)])
	else:
		data = read_data(fname)
	return data


if __name__ == '__main__':
	get_tRNA_data("HUMAN")
