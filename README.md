# Annotation of Mass2Motif mass fragments and neutral losses with MAGMa

In the Figure below a complete schematic representation of the whole pipeline is depicted. The red blocks with the rounded edges represent the final input and output. 
The blue rectangular blocks represent tools or websites used in the pipeline, and the yellow blocks with a wavy bottom represent intermediate inputs and outputs.

![Workflow_MSc_thesis - Flowchart](https://user-images.githubusercontent.com/107037630/208470162-1e391fda-d40f-4d01-830c-8cd312124ea9.jpeg)

## Input


## MS2Query

see https://github.com/iomega/ms2query for installation and run instructions

### Run script

e.g. 

## Select Mass2Motifs

### Prepare environment

conda install -c conda-forge rdkit

### Install tools

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

e.g. python3 massql_test.py 
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

## Output
