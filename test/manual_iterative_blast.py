#!/usr/bin/env python3

import rt_primer_design


refseq = 'XM_015786920'
strict_parameters = 'defined elsewhere'

test = rt_primer_design.iterate_primer_blast(
    ref_seq=refseq,
    starting_parameters=strict_parameters,
    verbose=True)

refseq_list = ["XM_015776733",
               "XM_015778389",
               "XM_015769405",
               "XM_015786920",
               "XM_015771259",
               "XM_015795930"]

# something strange about XM_015771259

refseq_blast_results = []
for rs in refseq_list:
    attempt = 0
    while attempt < 3:
        attempt += 1
        try:
            print("%s attempt %i" % (rs, attempt))
            result = rt_primer_design.iterate_primer_blast(
                ref_seq=rs,
                starting_parameters=strict_parameters,
                verbose=True)
            refseq_blast_results.append(result)
            attempt = 3
        except Exception as e:
            print("%s attempt %i FAILED" % (rs, attempt))

os.mkdir('output/2017-05-10_reruns')

with open('output/2017-05-10_reruns/primers.csv', 'w') as file:
    file.write('ref_seq,status,F,TM_F,R,TM_R,product_size,intron_size\n')
    for line in [x.csv_line() for x in refseq_blast_results]:
        file.write(line + '\n')

for primerset in refseq_blast_results:
    primerset.print_file(output_subdirectory='output/2017-05-10_reruns/html')

with open('output/2017-05-10_reruns/html/links.html', 'w') as file:
    file.write(html_header)
    for primer_set in refseq_blast_results:
        result_link = (
            """<a href="{}">{}</a><br />\n""".format(
                primer_set.url, primer_set.ref_seq))
        file.write(result_link)
    file.write(html_footer)