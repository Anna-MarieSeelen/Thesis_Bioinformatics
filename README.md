# Annotation of Mass2Motif mass fragments and neutral losses with MAGMa

In the Figure below a schematic representation of the whole pipeline is depicted. The red blocks with the rounded edges represent the final input and output. 
The blue rectangular blocks represent tools or websites used in the pipeline, and the yellow blocks with a wavy bottom represent intermediate inputs and outputs. The picture files generated by some scripts are not included in this pipeline.

![Workflow_MSc_thesis - Flowchart](https://user-images.githubusercontent.com/107037630/208470162-1e391fda-d40f-4d01-830c-8cd312124ea9.jpeg)

## Input

MS2LDA was run through the GNPS website on the three .mzML files containing the measured MS2 spectra to generate Mass2Motifs. To run MS2LDA on GNPS, first a classical molecular network was generated on GNPS. After the MS2LDA analysis on GNPS was finished, the .dict file containing the information obtained through the MS2LDA analysis on GNPS (e.g. Mass2Motifs, Mass2Motif fragments or losses, etc.) was uploaded on MS2LDA.org using the upload tab in the create experiment option. From MS2LDA.org the .csv containing the extracted fragment and loss Mass2Motif fragments or losses, and the .csv containing all fragmentation spectra and Mass2Motifs matching details were downloaded. The consensus spectra in .mgf format from the classical molecular network and the two .csv files from MS2LDA with the Mass2Motifs, the Mass2Motif fragments or losses, and the spectrum identifiers of experimental spectra that contained certain Mass2Motif fragments or losses were used as an input for the pipeline.

## MS2Query

see https://github.com/iomega/ms2query for installation and run instructions

A separate conda environment was made to run this script. This environment included the following packages:
However, it should be noted that not all these packages are neccessary to run the script!!

Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
absl-py                   1.1.0                    pypi_0    pypi
astunparse                1.6.3                    pypi_0    pypi
ca-certificates           2022.4.26            h06a4308_0  
cachetools                5.2.0                    pypi_0    pypi
certifi                   2022.6.15        py38h06a4308_0  
charset-normalizer        2.1.0                    pypi_0    pypi
cycler                    0.11.0                   pypi_0    pypi
deprecated                1.2.13                   pypi_0    pypi
flatbuffers               1.12                     pypi_0    pypi
fonttools                 4.33.3                   pypi_0    pypi
gast                      0.4.0                    pypi_0    pypi
gensim                    4.2.0                    pypi_0    pypi
google-auth               2.9.0                    pypi_0    pypi
google-auth-oauthlib      0.4.6                    pypi_0    pypi
google-pasta              0.2.0                    pypi_0    pypi
grpcio                    1.47.0                   pypi_0    pypi
h5py                      2.10.0                   pypi_0    pypi
idna                      3.3                      pypi_0    pypi
importlib-metadata        4.12.0                   pypi_0    pypi
joblib                    1.1.0                    pypi_0    pypi
keras                     2.9.0                    pypi_0    pypi
keras-preprocessing       1.1.2                    pypi_0    pypi
kiwisolver                1.4.3                    pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
libclang                  14.0.1                   pypi_0    pypi
libffi                    3.3                  he6710b0_2  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
llvmlite                  0.38.1                   pypi_0    pypi
lxml                      4.9.0                    pypi_0    pypi
markdown                  3.3.7                    pypi_0    pypi
matchms                   0.13.0                   pypi_0    pypi
matplotlib                3.5.2                    pypi_0    pypi
ms2deepscore              0.2.3                    pypi_0    pypi
ms2query                  0.3.2                    pypi_0    pypi
ncurses                   6.3                  h7f8727e_2  
networkx                  2.8.4                    pypi_0    pypi
numba                     0.55.2                   pypi_0    pypi
numpy                     1.22.4                   pypi_0    pypi
oauthlib                  3.2.0                    pypi_0    pypi
openssl                   1.1.1p               h5eee18b_0  
opt-einsum                3.3.0                    pypi_0    pypi
packaging                 21.3                     pypi_0    pypi
pandas                    1.4.3                    pypi_0    pypi
pillow                    9.1.1                    pypi_0    pypi
pip                       21.2.4           py38h06a4308_0  
protobuf                  3.19.4                   pypi_0    pypi
pyasn1                    0.4.8                    pypi_0    pypi
pyasn1-modules            0.2.8                    pypi_0    pypi
pyparsing                 3.0.9                    pypi_0    pypi
pyteomics                 4.5.3                    pypi_0    pypi
python                    3.8.13               h12debd9_0  
python-dateutil           2.8.2                    pypi_0    pypi
pytz                      2022.1                   pypi_0    pypi
readline                  8.1.2                h7f8727e_1  
requests                  2.28.1                   pypi_0    pypi
requests-oauthlib         1.3.1                    pypi_0    pypi
rsa                       4.8                      pypi_0    pypi
scikit-learn              1.1.1                    pypi_0    pypi
scipy                     1.8.1                    pypi_0    pypi
setuptools                61.2.0           py38h06a4308_0  
six                       1.16.0                   pypi_0    pypi
smart-open                6.0.0                    pypi_0    pypi
spec2vec                  0.6.0                    pypi_0    pypi
sqlite                    3.38.5               hc218d9a_0  
tensorboard               2.9.1                    pypi_0    pypi
tensorboard-data-server   0.6.1                    pypi_0    pypi
tensorboard-plugin-wit    1.8.1                    pypi_0    pypi
tensorflow                2.9.1                    pypi_0    pypi
tensorflow-estimator      2.9.0                    pypi_0    pypi
tensorflow-io-gcs-filesystem 0.26.0                   pypi_0    pypi
termcolor                 1.1.0                    pypi_0    pypi
threadpoolctl             3.1.0                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
tqdm                      4.64.0                   pypi_0    pypi
typing-extensions         4.2.0                    pypi_0    pypi
urllib3                   1.26.9                   pypi_0    pypi
werkzeug                  2.1.2                    pypi_0    pypi
wheel                     0.37.1             pyhd3eb1b0_0  
wrapt                     1.14.1                   pypi_0    pypi
xz                        5.2.5                h7f8727e_1  
zipp                      3.8.0                    pypi_0    pypi
zlib                      1.2.12               h7f8727e_2  

## Select Mass2Motifs

### Prepare environment

A separate conda environment was made to run this script. This environment included the following packages:
However, it should be noted that not all these packages are neccessary to run the script!!

Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
ca-certificates           2022.4.26            h06a4308_0  
certifi                   2022.6.15        py39h06a4308_0  
charset-normalizer        2.1.1                    pypi_0    pypi
contourpy                 1.0.5                    pypi_0    pypi
cycler                    0.11.0                   pypi_0    pypi
defusedxml                0.7.1                    pypi_0    pypi
deprecated                1.2.13                   pypi_0    pypi
fonttools                 4.37.2                   pypi_0    pypi
fpdf                      1.7.2                    pypi_0    pypi
fpdf2                     2.5.5                    pypi_0    pypi
idna                      3.4                      pypi_0    pypi
kiwisolver                1.4.4                    pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
llvmlite                  0.39.1                   pypi_0    pypi
lxml                      4.9.1                    pypi_0    pypi
matchms                   0.16.0                   pypi_0    pypi
matplotlib                3.6.0                    pypi_0    pypi
ncurses                   6.3                  h5eee18b_3  
networkx                  2.8.6                    pypi_0    pypi
numba                     0.56.2                   pypi_0    pypi
numpy                     1.23.0                   pypi_0    pypi
openssl                   1.1.1p               h5eee18b_0  
packaging                 21.3                     pypi_0    pypi
pandas                    1.4.3                    pypi_0    pypi
pickydict                 0.4.0                    pypi_0    pypi
pillow                    9.2.0                    pypi_0    pypi
pip                       21.2.4           py39h06a4308_0  
pyparsing                 3.0.9                    pypi_0    pypi
pyteomics                 4.5.5                    pypi_0    pypi
python                    3.9.12               h12debd9_1  
python-dateutil           2.8.2                    pypi_0    pypi
pytz                      2022.1                   pypi_0    pypi
rdkit                     2022.3.4                 pypi_0    pypi
readline                  8.1.2                h7f8727e_1  
requests                  2.28.1                   pypi_0    pypi
scipy                     1.9.1                    pypi_0    pypi
setuptools                59.8.0                   pypi_0    pypi
six                       1.16.0                   pypi_0    pypi
sqlite                    3.38.5               hc218d9a_0  
tk                        8.6.12               h1ccaba5_0  
tzdata                    2022a                hda174b7_0  
urllib3                   1.26.12                  pypi_0    pypi
wheel                     0.37.1             pyhd3eb1b0_0  
wrapt                     1.14.1                   pypi_0    pypi
xz                        5.2.5                h7f8727e_1  
zlib                      1.2.12               h7f8727e_2

### Install tools

conda install -c conda-forge rdkit

### Run script
e.g. python3 select_Mass2Motif_frag_and_loss.py
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/MS2Query_output.csv 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/MS2LDA_spectra_and_motif.csv 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/MS2LDA_motif_and_fragments.csv 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/input/consensus_spectra_from_GNPS_classical_molecular_network.mgf 
/lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/select_Mass2Motifs/output

## MassQL

### Prepare environment

A separate conda environment was made to run this script. This environment included the following packages:
However, it should be noted that not all these packages are neccessary to run the script!!

Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
aiosignal                 1.2.0                    pypi_0    pypi
attrs                     21.4.0                   pypi_0    pypi
ca-certificates           2022.4.26            h06a4308_0  
certifi                   2022.6.15        py38h06a4308_0  
charset-normalizer        2.1.0                    pypi_0    pypi
click                     8.0.4                    pypi_0    pypi
cycler                    0.11.0                   pypi_0    pypi
deprecated                1.2.13                   pypi_0    pypi
distlib                   0.3.5                    pypi_0    pypi
filelock                  3.7.1                    pypi_0    pypi
fonttools                 4.34.4                   pypi_0    pypi
frozenlist                1.3.0                    pypi_0    pypi
greenlet                  1.1.2                    pypi_0    pypi
grpcio                    1.43.0                   pypi_0    pypi
idna                      3.3                      pypi_0    pypi
importlib-resources       5.8.0                    pypi_0    pypi
jsonschema                4.7.2                    pypi_0    pypi
kaleido                   0.2.1                    pypi_0    pypi
kiwisolver                1.4.3                    pypi_0    pypi
lark-parser               0.12.0                   pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
llvmlite                  0.38.1                   pypi_0    pypi
lxml                      4.9.1                    pypi_0    pypi
massql                    0.0.12                   pypi_0    pypi
matchms                   0.13.0             pyh5e36f6f_0    https://anaconda.org/bioconda/matchms/0.13.0/download
matplotlib                3.5.2                    pypi_0    pypi
msgpack                   1.0.4                    pypi_0    pypi
ncurses                   6.3                  h5eee18b_3  
networkx                  2.8.4                    pypi_0    pypi
numba                     0.55.2                   pypi_0    pypi
numpy                     1.22.4                   pypi_0    pypi
openssl                   1.1.1q               h7f8727e_0  
packaging                 21.3                     pypi_0    pypi
pandas                    1.4.3                    pypi_0    pypi
pickydict                 0.4.0                    pypi_0    pypi
pillow                    9.2.0                    pypi_0    pypi
pip                       22.1.2           py38h06a4308_0  
platformdirs              2.5.2                    pypi_0    pypi
plotly                    5.9.0                    pypi_0    pypi
protobuf                  3.20.1                   pypi_0    pypi
psims                     1.0.1                    pypi_0    pypi
py-expression-eval        0.3.14                   pypi_0    pypi
pyarrow                   8.0.0                    pypi_0    pypi
pydot                     1.4.2                    pypi_0    pypi
pymzml                    2.5.1                    pypi_0    pypi
pyparsing                 3.0.9                    pypi_0    pypi
pyrsistent                0.18.1                   pypi_0    pypi
pyteomics                 4.5.3                    pypi_0    pypi
python                    3.8.13               h12debd9_0  
python-dateutil           2.8.2                    pypi_0    pypi
pytz                      2022.1                   pypi_0    pypi
pyyaml                    6.0                      pypi_0    pypi
ray                       1.13.0                   pypi_0    pypi
readline                  8.1.2                h7f8727e_1  
regex                     2022.7.9                 pypi_0    pypi
requests                  2.28.1                   pypi_0    pypi
scipy                     1.8.1                    pypi_0    pypi
setuptools                61.2.0           py38h06a4308_0  
six                       1.16.0                   pypi_0    pypi
sqlalchemy                1.4.39                   pypi_0    pypi
sqlite                    3.38.5               hc218d9a_0  
tenacity                  8.0.1                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
tqdm                      4.64.0                   pypi_0    pypi
urllib3                   1.26.10                  pypi_0    pypi
virtualenv                20.15.1                  pypi_0    pypi
wheel                     0.37.1             pyhd3eb1b0_0  
wrapt                     1.14.1                   pypi_0    pypi
xz                        5.2.5                h7f8727e_1  
zipp                      3.8.1                    pypi_0    pypi
zlib                      1.2.12               h7f8727e_2  

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

A separate conda environment was made to run this script. This environment included the following packages:
However, it should be noted that not all these packages are neccessary to run the script!!

Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
blas                      1.0                         mkl  
bzip2                     1.0.8                h7b6447c_0  
c-ares                    1.18.1               h7f8727e_0  
ca-certificates           2022.9.24            ha878542_0    conda-forge
cairo                     1.16.0               h19f5f5c_2  
cairocffi                 1.3.0              pyhd8ed1ab_0    conda-forge
cairosvg                  2.5.2                    pypi_0    pypi
certifi                   2022.6.15        py38h06a4308_0  
cffi                      1.15.1           py38h74dc2b5_0  
charset-normalizer        2.1.1                    pypi_0    pypi
contourpy                 1.0.5                    pypi_0    pypi
coverage                  6.3.2            py38h7f8727e_0  
cssselect2                0.7.0                    pypi_0    pypi
curl                      7.84.0               h5eee18b_0  
cycler                    0.11.0                   pypi_0    pypi
cython                    0.29.30          py38h6a678d5_0  
defusedxml                0.7.1                    pypi_0    pypi
deprecated                1.2.13                   pypi_0    pypi
expat                     2.4.4                h295c915_0  
fontconfig                2.13.1               h6c09931_0  
fonttools                 4.37.3                   pypi_0    pypi
freetype                  2.10.4               h0708190_1    conda-forge
gettext                   0.21.0               hf68c758_0  
git                       2.34.1          pl5262hc120c5b_0  
glib                      2.69.1               h4ff587b_1  
greenlet                  2.0.0a2                  pypi_0    pypi
icu                       58.2                 he6710b0_3  
idna                      3.3                      pypi_0    pypi
importlib-metadata        4.12.0                   pypi_0    pypi
intel-openmp              2021.4.0          h06a4308_3561  
kiwisolver                1.4.4                    pypi_0    pypi
krb5                      1.19.2               hac12032_0  
ld_impl_linux-64          2.38                 h1181459_1  
libcurl                   7.84.0               h91b91d3_0  
libedit                   3.1.20210910         h7f8727e_0  
libev                     4.33                 h7f8727e_1  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libnghttp2                1.46.0               hce63b2e_0  
libpng                    1.6.37               hbc83047_0  
libssh2                   1.10.0               h8f2d780_0  
libstdcxx-ng              11.2.0               h1234567_1  
libuuid                   1.0.3                h7f8727e_2  
libxcb                    1.15                 h7f8727e_0  
libxml2                   2.9.14               h74e7548_0  
libxslt                   1.1.35               h4e12654_0  
llvmlite                  0.39.1                   pypi_0    pypi
lxml                      4.9.1            py38h1edc446_0  
macauthlib                0.6.0                    pypi_0    pypi
magma                     1.3                       dev_0    <develop>
matchms                   0.16.0                   pypi_0    pypi
matplotlib                3.6.0                    pypi_0    pypi
mkl                       2021.4.0           h06a4308_640  
mkl-service               2.4.0            py38h7f8727e_0  
mkl_fft                   1.3.1            py38hd3c417c_0  
mkl_random                1.2.2            py38h51133e4_0  
mock                      4.0.3                    pypi_0    pypi
ncurses                   6.3                  h5eee18b_3  
networkx                  2.8.6                    pypi_0    pypi
nose                      1.3.7           pyhd3eb1b0_1008  
numba                     0.56.2                   pypi_0    pypi
numpy                     1.23.1           py38h6c91a56_0  
numpy-base                1.23.1           py38ha15fc14_0  
openssl                   1.1.1q               h7f8727e_0  
packaging                 21.3                     pypi_0    pypi
pandas                    1.5.0                    pypi_0    pypi
pcre                      8.45                 h9c3ff4c_0    conda-forge
pcre2                     10.37                he7ceb23_1  
perl                      5.26.2               h14c3975_0  
pickydict                 0.4.0                    pypi_0    pypi
pillow                    9.2.0                    pypi_0    pypi
pip                       22.1.2           py38h06a4308_0  
pixman                    0.40.0               h36c2ea0_0    conda-forge
pp                        1.6.4.4                  pypi_0    pypi
pycparser                 2.21               pyhd8ed1ab_0    conda-forge
pyparsing                 3.0.9                    pypi_0    pypi
pyteomics                 4.5.5                    pypi_0    pypi
python                    3.8.13               h12debd9_0  
python-dateutil           2.8.2                    pypi_0    pypi
pytz                      2022.4                   pypi_0    pypi
rdkit                     2022.3.5                 pypi_0    pypi
readline                  8.1.2                h7f8727e_1  
requests                  2.28.1                   pypi_0    pypi
scipy                     1.9.1                    pypi_0    pypi
setuptools                59.8.0                   pypi_0    pypi
six                       1.16.0             pyhd3eb1b0_1  
sqlalchemy                1.4.41                   pypi_0    pypi
sqlite                    3.39.2               h5082296_0  
tinycss2                  1.2.1                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
urllib3                   1.26.12                  pypi_0    pypi
webencodings              0.5.1                    pypi_0    pypi
webob                     1.8.7                    pypi_0    pypi
wheel                     0.37.1             pyhd3eb1b0_0  
wrapt                     1.14.1                   pypi_0    pypi
xz                        5.2.5                h7f8727e_1  
zipp                      3.8.1                    pypi_0    pypi
zlib                      1.2.12               h7f8727e_2  

![image](https://user-images.githubusercontent.com/107037630/208644077-c91e8a1e-e96b-4de5-a505-b6de630c91b4.png)
![image](https://user-images.githubusercontent.com/107037630/208644146-935b61b0-6aa4-4d20-9f6e-3305899de913.png)
![image](https://user-images.githubusercontent.com/107037630/208644254-02dcd1f5-fd77-4502-870a-1fab0422c9c5.png)

### Install tools

see https://github.com/NLeSC/MAGMa/tree/master/job
conda install -c conda-forge rdkit

### Run script

e.g. python3 MAGMa_final.py /lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MassQL/output/out_spectrum /lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MAGMa/output/MAGMa_results_database_for_every_spectrum_from_massql /home/seele006/thesis/motif_massql_querries.txt /lustre/BIF/nobackup/seele006/MSc_thesis_annotation_Mass2Motif_fragments_data/MAGMa/output/pic_mass2Motif_frag

## Output

The annotations for each Mass2Motif mass fragment and neutral loss were combined in a tsv-formatted output file. The frequency that each SMILES annotation for a Mass2Motif fragment or loss was obtained, was also tracked in the .tsv file. If the molecular weight of the SMILES annotated to the neutral loss was not similar 1 decimal after the comma to the weight of the Mass2Motif neutral loss, the molecular weight of the SMILES structure was also written to the .tsv file.
