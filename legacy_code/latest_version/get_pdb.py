#!/usr/bin/env python
#
# Provides simple functionallity to download pdb files using python.
# Returns the path to the downloaded file

import os, urllib2, gzip

def get_pdb(pdb_id):
	fname = 'pdb/'+pdb_id+'.pdb'
	#check if pdb is present
	if os.path.exists(fname):
		return fname

	#check for pdb dir
	if not os.path.exists('pdb/'):
    		os.makedirs('pdb')	

	#download pbd.gz
	f = urllib2.urlopen("http://www.rcsb.org/pdb/files/"+pdb_id+".pdb.gz")
	g = open(fname+'.gz','w')
	while 1:
		packet = f.read()
		if not packet:
			break
		g.write(packet)
	f.close()
	g.close()

	#unzip
	f = gzip.open(fname+'.gz', 'r')
	g = open(fname,'w')
	g.write(f.read())
	f.close()
	g.close()

	os.remove('pdb/'+pdb_id+'.pdb.gz')

	return fname
