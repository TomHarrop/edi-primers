#!/usr/bin/env python3

import csv
import re
from Bio import Entrez

rap_id_file = 'test/rap_ids.csv'
Entrez.email = my_email_address

#############
# FUNCTIONS #
#############

# get the GIs from RAP
def get_gis_from_rap(rap):
    handle = Entrez.esearch(
        db='gene',
        term=rap + ' AND txid39947[ORGN]')
    record = Entrez.read(handle)
    return(record['IdList'])


# download the record for a gi
def download_gi_record(gi):
    handle = Entrez.efetch(
        db='gene',
        id=gi,
        retmode='xml')
    return(handle.read())


# parse the XM
def parse_xm_from_record(record):
    xm_search = re.compile(r'XM_\d+')
    if xm_search.search(record):
        return(xm_search.search(record).group(0))

# parse the NM
def parse_nm_from_record(record):
    nm_search = re.compile(r'NM_\d+')
    if nm_search.search(record):
        return(nm_search.search(record).group(0))


# csv function
def write_dict_to_csv(result_list, file_name):
    with open(file_name, 'w') as csvfile:
        headers = result_list[0].keys()
        dict_writer = csv.DictWriter(csvfile, fieldnames=headers)
        dict_writer.writeheader()
        dict_writer.writerows(result_list)


########
# CODE # 
########

# read RAP ID from file
with open(rap_id_file, 'r') as f:
    csvreader = csv.reader(f)
    rap_ids = [x[0] for x in csvreader]

rap_to_gi_list = []
for rap in rap_ids:
    gi_list = get_gis_from_rap(rap)
    for gi in gi_list:
        rap_to_gi_list.append({'rap': rap, 'gi': gi})
    
# get XMs and NMs for each gi
gi_list = [x['gi'] for x in rap_to_gi_list]
gi_to_nm = []
gi_to_xm = []
for gi in gi_list:
    record = download_gi_record(gi)
    gi_to_xm.append({'gi': gi, 'transcript': parse_xm_from_record(record)})
    gi_to_nm.append({'gi': gi, 'transcript': parse_nm_from_record(record)})

# write to csv
write_dict_to_csv(rap_to_gi_list, 'test/rap_to_gi.csv')
write_dict_to_csv(gi_to_nm, 'test/gi_to_nm.csv')
write_dict_to_csv(gi_to_xm, 'test/gi_to_xm.csv')

