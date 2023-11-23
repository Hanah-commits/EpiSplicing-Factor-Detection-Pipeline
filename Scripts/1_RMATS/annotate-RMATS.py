import os

exon_types = ['AS', 'CS']
event_types = ['SE', 'MXE']
flanks = ['50','100','200']

for exon_type in exon_types:
    for type in event_types:

        for flank in flanks:

            os.system(f'bedtools intersect -loj -s -a 0_Files/flanks{flank}.bed -b 0_Files/{type}_{exon_type}.bed | sort | uniq > 0_Files/{type}_flanks{flank}_{exon_type}.bed')
