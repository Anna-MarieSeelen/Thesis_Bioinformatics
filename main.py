

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
from sys import argv
import re
import subprocess
import os.path
from ms2query.run_ms2query import download_default_models, default_library_file_base_names, run_complete_folder
from ms2query.ms2library import create_library_object_from_one_dir

# functions
def ms2query():
    # Set the location where all your downloaded model files are stored
    ms2query_library_files_directory = "./ms2query_library_files"
    # Define the folder in which your query spectra are stored.
    # Accepted formats are: "mzML", "json", "mgf", "msp", "mzxml", "usi" or a pickled matchms object.
    ms2_spectra_directory = "specify_directory"

    # Downloads pretrained models and files for MS2Query (>10GB download)
    download_default_models(ms2query_library_files_directory, default_library_file_base_names())

    # Create a MS2Library object
    ms2library = create_library_object_from_one_dir(ms2query_library_files_directory, default_library_file_base_names())

    # Run library search and analog search on your files.
    run_complete_folder(ms2library, ms2_spectra_directory)

def run_needle(reference, related, gapopen, gapextend=0.5):
    """Running the command line tool needle

    reference:
    related:
    gapopen:
    gapextend:
    return:
    """
    out_fn = "out.needle"
    if os.path.exists(out_fn):
        return out_fn
    cmd = 'needle {0} {1} -gapopen {2} -gapextend {3} -outfile {4}' \
        .format(reference, related, gapopen, gapextend, out_fn)
    e = subprocess.check_call(cmd, shell=True)
    print("EXIT STATUS AND TYPE", e, type(e))
    return out_fn

def run_tool2_when_you_want_output_of_command_line():
    """Run jellyfish program on fasta file

    input_fn: string, filename of input FASTA file
    kmer_size: int, size of k-mers used by jellyfish
    """
    out_fn = 'tomato{}'.format(kmer_size)
    cmd = 'jellyfish count -m {} -s 1000000 -o {} {}'\
        .format(kmer_size, out_fn, input_fn)
    e = subprocess.check_output(cmd, shell=True)
    if os.path.exists(out_fn):
        cmd = 'jellyfish stats {}'.format(out_fn)
    else:
        cmd = 'jellyfish stats {}_0'.format(out_fn)
    res = subprocess.check_output(cmd, shell=True)
    return res

def main():
    """Main function of this module"""
    # step 1: parse the protein sequence from file 1 and file 2 into dict
    reference = argv[1]
    related = argv[2]
    dict = parse_input(reference, related)
    # step 2: determine the lenght of the sequences
    lenght_dict = lenght_seq(dict)
    # step 4: align protein sequences from file 1 to the other species in file 2
    try:
        gapopen = int(argv[3])  # try to see if gap_open is in the commandline, otherwise gap_open is 8.
    except:
        gapopen = 8
    out_fn = run_needle(reference, related, gapopen, gapextend=0.5)
    # step 5: parse needle output to extract pairwise alignments
    records = record_finder(out_fn)
    print("{}\t{}\t{}\t{}\t{}\t{}".format("NAME REF", "LENGHT", "NAME REL", "LENGHT", "HAMDIST", "IDEN"))
    for i in range(len(records)):
        clean_record = extract_alignments(records[i])
        list_of_string = (make_long_string(clean_record))
        # step 5: calculate hamming distance between pairwise alignments
        h_distance = hamming_distance(list_of_string)
        reference_name = list(lenght_dict)[0]
        related_name = list(lenght_dict)[i + 1]
        per_identity = percent_identity(h_distance, list_of_string)
        # step 7: tab delimited format for with the alignments, containing: sequence
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5:.1f}".format(reference_name, lenght_dict[reference_name], related_name,
                                                        lenght_dict[related_name], h_distance, per_identity))


if __name__ == "__main__":
    main()
