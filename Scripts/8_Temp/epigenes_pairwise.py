import pandas as pd
import os
import sys
import json

mode = 'flanks' #sys.argv[1]

with open('paths_multi.json') as f:
        data = json.load(f)

op_dirs = []
for process in data['list_of_processes']:
    dir = data[process]['Output directory']
      
    os.system(f'cp -r {dir}0_Files ./0_Files')

    with open('paths.json', 'w') as out:
        json.dump(data[process], out, indent=4)

    # get epigenes for current analysis
    os.system('python 8_Temp/get_epi_nonepi_flanks.py')

    # pre-rbp
    os.system('mkdir ../RBPmap')
    os.system('python 8_Temp/pre_rbp.py 0')

    # run rbpmap
    os.system('python 7_Post/run_rbpmap.py')

    # post-rbpmap
    os.system('python 8_Temp/post-rbp.py 0')

    # get nonepi flanks
    os.system(f'cp ./Post-processing/nonepi_flanks.bed ./0_Files/Post-processing/')
    os.system(f'cp ./Post-processing/FilteredZscores_nonepi.csv ./0_Files/Post-processing/')

    # features

    # classifier

    #DEA
