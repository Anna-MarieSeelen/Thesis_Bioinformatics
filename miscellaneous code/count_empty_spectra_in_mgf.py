

#/lustre/BIF/nobackup/seele006/final_thesis_data

import re
from sys import argv

def parse_input(mgf_file: str) -> dict:
    """Parses mgf into strings, where each string is a spectrum and stores those in dict with the spectrum_id as key.

    :param mgf_file: str, name of mgf formatted file containing MS/MS spectra from GNPS
    :return: dictionary with {spectrum_id:record} where each record is a string containing the mgf-style accession of one
    compound
    """

    lines_mgf_file=open(mgf_file)
    spectrum_record_bool = False
    mgf_spectrum_records=[]
    mgf_spectrum_record = ""
    for line in lines_mgf_file:
        if line.startswith("BEGIN IONS"):
            spectrum_record_bool=True
        if spectrum_record_bool:
            mgf_spectrum_record+=line
        if line.startswith("END IONS"):
            mgf_spectrum_records.append(mgf_spectrum_record)
            spectrum_record_bool = False
            mgf_spectrum_record=""

    dict_with_mgf_spectra={}
    for mgf_spectrum_record in mgf_spectrum_records:
        mgf_spectrum_record=mgf_spectrum_record.strip()
        # look for the spectrumid in the string and use it as a key for the dict
        key = str(re.search(r'SCANS=(.*)', mgf_spectrum_record).group(1))
        if key is not None:
            dict_with_mgf_spectra[key] = mgf_spectrum_record
    print(dict_with_mgf_spectra["19"])
    return dict_with_mgf_spectra

def count_empty_spectra(dict_with_mgf_spectra):
    count_empty=0
    count_total=len(dict_with_mgf_spectra)
    for key in dict_with_mgf_spectra.keys():
        count_new_line=dict_with_mgf_spectra[key].count("\n")
        if count_new_line<11:
            count_empty+=1
    print(count_empty)
    print(count_total)

def main():
    """Main function of this module"""
    mgf_file = argv[1]
    dict_mgf=parse_input(mgf_file)
    count_empty_spectra(dict_mgf)

if __name__ == "__main__":
    main()

