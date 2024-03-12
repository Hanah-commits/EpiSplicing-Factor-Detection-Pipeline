import time
import os
import json
from pathlib import Path
import shutil


def check_args():

    with open('paths.json') as f:
        d = json.load(f)

    args = {
    "RNASeq files" : d['RNASeq files'],
    "Reference GFF3" : d['Reference GFF3'],
    "Reference GTF" : d['Reference GTF'],
    "Reference fasta": d['Reference fasta'],
    "tissue1" : d["tissue1"],
    "tissue2" : d["tissue2"],
    "Histone modifications" : d["Histone modifications"],
    "ChIPSeq files" : d["ChIPSeq files"],
    "MAJIQ config" : d["MAJIQ config"],
    "RBPmap directory" : d["RBPmap directory"],
    "threads" : d['threads'],
    "read_length":  d['read_length'],
    "Output directory" : d['Output directory'],
    "DEXSEQ directory": d['DEXSEQ directory'],
    "RMATS directory": d['RMATS directory'],
    }

    #check histone modifications & tissue names
    args["Histone modifications"] = [x for x in args["Histone modifications"] if x]
    if len(args["Histone modifications"]) == 0:
        raise ValueError('No Histone Modifications')

    if len(args["tissue1"].strip()) == 0 or len(args["tissue2"].strip()) == 0:
        raise ValueError('No Tissue Name(s)')
        
    #check datatypes
    for val in ['threads', 'read_length']:
        try:
            int(args[val])
        except:
            raise ValueError('Invalid Input ' + args[val])


    # check path validity of directories
    dirs = ["RNASeq files", "ChIPSeq files", "RBPmap directory", "DEXSEQ directory", "RMATS directory"]
    if len(args['Output directory']) != 0:
        dirs.append('Output directory')
    dir_paths = []
    for dir in dirs:

        if args[dir][-1] != '/' and len(args[dir]) > 0:
            args[dir] += '/'

        dir_paths.append(args[dir])
        

    for path in dir_paths:
        
        if not os.path.exists(os.path.dirname(path)):
            raise ValueError('Path does not exist ' + path)

    # check path validity of files
    file_paths = [args["Reference GFF3"], args["Reference GTF"], args["Reference fasta"], args["MAJIQ config"]]
    for file in file_paths:
        
        if not os.path.isfile(file):
            raise ValueError('File does not exist ' + file)

        
    # check if temp directories already exist
    temp_dirs = ['0_Files/', '../RBPmap/']
    for dir in temp_dirs:
        if os.path.exists(dir):
            # temp dir not empty
            if len(os.listdir(dir)) != 0:
                raise ValueError('Delete or move directory to another location ' + dir)
        else: 
            #create temp dir
            Path(dir).mkdir(parents=True, exist_ok=True)

    # check if log files already exist
    log_file_name = 'output.log'
    try:
        f = open(log_file_name, 'x')
    except FileExistsError:
        raise ValueError('Delete or move the log file to another location ' + log_file_name)

    output_dir = args['Output directory']
    if len(output_dir) == 0:
        # create custome output directory tissue1_tissue2_timestamp
        output_dir = str(Path(os.getcwd()).parent.absolute()) + "/Output/"+ args["tissue1"]+ "_" + args["tissue2"] + "_" + str(time.time()) +"/"
        Path(output_dir).mkdir(parents=True, exist_ok=True)


    with open('paths.json', 'w') as fp:
        json.dump(args, fp, indent=4)

    # copy input arguments (paths.json) to output_dir
    shutil.copyfile('paths.json', output_dir+'paths.json')

    return output_dir


def move_dirs(output_dir):
    shutil.move('0_Files/', output_dir)
    shutil.move('../RBPmap/', output_dir)
    shutil.move('output.log', output_dir)
