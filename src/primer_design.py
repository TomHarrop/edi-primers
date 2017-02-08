#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################
# DEPENDENCIES #
################

import rt_primer_design.__init__
import sys
import os

##############
# PARAMETERS #
##############

# initial blast settings
strict_parameters = {
    'PRIMER_PRODUCT_MIN': '70',
    'PRIMER_PRODUCT_MAX': '180',
    'PRIMER_NUM_RETURN': '10',
    'PRIMER_MIN_TM': '55.0',
    'PRIMER_OPT_TM': '60.0',
    'PRIMER_MAX_TM': '65.0',
    'PRIMER_MAX_DIFF_TM': '5',
    'MIN_INTRON_SIZE': '0',
    'MAX_INTRON_SIZE': '1000000',
    'PRIMER_SPECIFICITY_DATABASE': 'refseq_mrna',
    'EXCLUDE_ENV': 'on',
    'ORGANISM': 'Oryza sativa Japonica Group (taxid:39947)',
    'TOTAL_MISMATCH_IGNORE': '7',
    'ALLOW_TRANSCRIPT_VARIANTS': 'on',
    'MAX_CANDIDATE_PRIMER': '1000',
    'PRIMER_MIN_GC': '45.0',
    'PRIMER_MAX_GC': '55.0',
    'GC_CLAMP': '2',
    'POLYX': '3',
    'SELF_ANY': '3.00',
    'SELF_END': '1.00',
    'SEARCH_SPECIFIC_PRIMER': 'on',
    'SHOW_SVIEWER': 'on',
    'UNGAPPED_BLAST': 'on',
    'LOW_COMPLEXITY_FILTER': 'on',
    'SHOW_SVIEWER': 'on',
    'SPAN_INTRON': 'on'
}

# use NCBI usage guidelines:

# Do not poll for any single RID more often than once a minute.
wait_seconds = 60

# Do not contact the server more often than once every three seconds.
max_jobs = int(wait_seconds/3)

# have to increase the recursion limit for joblib
sys.setrecursionlimit(10000)

# tell NCBI who we are
my_email = 'tom.harrop@otago.ac.nz'
strict_parameters['EMAIL'] = my_email

# load the list of refseq IDs
with open('data/refseq.txt', 'r') as file:
    refseq_list = [x.strip() for x in file.readlines()]

# run the blast jobs
refseq_blast_results = rt_primer_design.multiple_primer_blast(
    ref_seq_list=refseq_list,
    starting_parameters=strict_parameters,
    wait_seconds=60,
    n_jobs=max_jobs)

# NM_001050201 has multiple similar seqs, test?

# write output
if not os.path.isdir('output'):
    os.mkdir('output')
if not os.path.isdir('output/html'):
    os.mkdir('output/html')


with open('output/primers.csv', 'w') as file:
    file.write('ref_seq,status,F,TM_F,R,TM_R,product_size,intron_size\n')
    for line in [x.csv_line() for x in refseq_blast_results]:
        file.write(line + '\n')

for primerset in refseq_blast_results:
    primerset.print_file(output_subdirectory='output/html')

html_header = """
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" \n
    "http://www.w3.org/TR/html4/strict.dtd">\n
    <html lang="en">\n
    <head>\n
    <meta http-equiv="content-type" content="text/html; charset=utf-8">\n
    <title>Links for primer design</title>\n
    <link rel="stylesheet" type="text/css" href="style.css">\n
    <script type="text/javascript" src="script.js"></script>\n
    </head>\n
    <body>\n
    <p>\n"""

html_footer = '\n</p>\n</body>\n</html>'

with open('output/html/links.html', 'w') as file:
    file.write(html_header)
    for primer_set in refseq_blast_results:
        result_link = (
            """<a href="{}">{}</a><br />\n""".format(
                primer_set.url, primer_set.ref_seq))
        file.write(result_link)
    file.write(html_footer)