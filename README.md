# Annotation of Mass2Motif mass fragments and neutral losses with MAGMa

In the Figure below a complete schematic representation of the whole pipeline is depicted. The red blocks with the rounded edges represent the final input and output. 
The blue rectangular blocks represent tools or websites used in the pipeline, and the yellow blocks with a wavy bottom represent intermediate inputs and outputs.

![Workflow_MSc_thesis - Flowchart](https://user-images.githubusercontent.com/107037630/208470162-1e391fda-d40f-4d01-830c-8cd312124ea9.jpeg)

## Input

MS2LDA was run through the GNPS website on the three .mzML files containing the measured MS2 spectra to generate Mass2Motifs. To run MS2LDA on GNPS, first a classical molecular network was generated on GNPS. After the MS2LDA analysis on GNPS was finished, the .dict file containing the information obtained through the MS2LDA analysis on GNPS (e.g. Mass2Motifs, Mass2Motif fragments or losses, etc.) was uploaded on MS2LDA.org using the upload tab in the create experiment option. From MS2LDA.org the .csv containing the extracted fragment and loss Mass2Motif fragments or losses, and the .csv containing all fragmentation spectra and Mass2Motifs matching details were downloaded. The consensus spectra in .mgf format from the classical molecular network and the two .csv files from MS2LDA with the Mass2Motifs, the Mass2Motif fragments or losses, and the spectrum identifiers of experimental spectra that contained certain Mass2Motif fragments or losses were used as an input for the pipeline.

## MS2Query

see https://github.com/iomega/ms2query for installation and run instructions

## Select Mass2Motifs

### Prepare environment

### Install tools

conda install -c conda-forge rdkit

### Run script
e.g. python3 make_pdf_with_smiles.py 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/MS2Query_output.csv 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/MS2LDA_spectra_and_motif.csv 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/MS2LDA_motif_and_fragments.csv 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/consensus_spectra_from_GNPS_classical_molecular_network.mgf 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/output

## MassQL

### Prepare environment

### Install tools

https://pypi.org/project/massql/

### Run script

e.g. python3 massql.py 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/input/motif_massql_querries.txt 
/mnt/LTR_userdata/hooft001/mass_spectral_embeddings/datasets/GNPS_15_12_21/ALL_GNPS_15_12_2021_positive_annotated.pickle 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/output/out_spectrum 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/output/out_files 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/output/json_enzo/GNPS.mgf 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/output/json_enzo/GNPS.json

## MAGMa

### Prepare environment

conda install -c conda-forge rdkit

### Install tools

see https://github.com/NLeSC/MAGMa/tree/master/job

### Run script

e.g. python3 MAGMa_final.py /lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/output/out_spectrum /lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MAGMa/output/MAGMa_results_database_for_every_spectrum_from_massql /home/seele006/thesis/motif_massql_querries.txt /lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MAGMa/output/pic_mass2Motif_frag

## Output

From the MAGMa
