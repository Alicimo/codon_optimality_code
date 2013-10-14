#!/usr/bin/env python

import BeautifulSoup as BS

soup = BS.BeautifulSoup(open('tRNA_Database.html'))
for link in soup.findAll('a',href=True):
	print link.get('href'),',', link.text
