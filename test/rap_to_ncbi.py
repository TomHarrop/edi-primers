#!/usr/bin/env python3

import csv
from Bio import Entrez
from Bio import SeqIO

rap_id_file = 'test/rap_ids.csv'
Entrez.email = my_email_address

# get entrez GI for a RAP ID
def get_gi(rap_id):
    term = rap_id
    handle = Entrez.esearch(db='nuccore', term=term)
    record = Entrez.read(handle)
    gi_list = record['IdList']
    return(gi_list)

# get genbank record from GI
def get_gb_records(gi_list):
    handle = Entrez.efetch(
        db='nuccore',
        id=gi_list,
        rettype='gb',
        retmode='text')
    return(SeqIO.parse(handle, 'gb'))

# read RAP ID from file
loc_to_rap = {}
with open(rap_id_file, 'r') as f:
    csvreader = csv.reader(f)
    for line in csvreader:
        loc_to_rap[line[0]] = line[1]

# search entrez for gene ID
test_gis = get_gi('OS01G0848400')

test_handle = Entrez.efetch(
    db='nuccore',
    id=test_gis[0],
    rettype='gb',
    retmode='text')



test_gbk = get_gb_records(test_gis)


