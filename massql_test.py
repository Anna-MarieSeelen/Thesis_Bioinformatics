#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: alignment of protein sequences to a reference protein sequence
and returning stats
Usage: python3 *name_of_script* *name_input_file* [*output_file*] [gap_penalty]
    name_of_script:
    name_input_file:
    output_file: (optional)
    gap_penalty: default is 8 (optional)
"""

# import statements
from massql import msql_engine
from sys import argv
import sys
#TODO: fix that it imports the file, because the massql function sucksss
#import ALL_GNPS_210409_positive_processed_annotated_CF_NPC_classes.txt
import ntpath
from pathlib import Path
import pyarrow.feather as feather
import os

# functions
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def try_massql(query, file):
    #sys.path.insert(1, '/lustre/BIF/nobackup/seele006/MS2Query_search_libraries')
    results=msql_engine.process_query(query,file)
    return results

def from_feather_to_df(results_feather):
    read_df = feather.read_feather(results_feather)
    print(read_df)
    return read_df

def main():
    """Main function of this module"""
    # step 1: parse the protein sequence from file 1 and file 2 into dict
    path=argv[1]
    name_file="FractionProfiling_RPOS_ToF10_PreCheck_LTR_01_DDA.mzML"
    print(find(path, name_file))
    #head, tail = ntpath.split(argv[1])
    #path_to_file_with_GNPS_mass_library_txt = argv[1]
    query=("QUERY scaninfo(MS2DATA) WHERE MS2NL=163")
    file = "FractionProfiling_RPOS_ToF10_PreCheck_LTR_01_DDA.mzML"
    results_feather=try_massql(query, file)
    from_feather_to_df(results_feather)


if __name__ == "__main__":
    main()