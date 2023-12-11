# TIGTOG: Taxonomic Information of Giant viruses using Trademark Orthologous Groups
TIGTOG is a command-line tool for assigning taxonomy to giant virus metagenome-assembled genomes (GVMAGs) using protein family trademarks. 

This tool will search predicted proteins of input MAGs against a set of curated Hidden Markov Models (HMM) for protein families prevalent in different groups of giant viruses. It utilizes the protein profiles and other sequence features (e.g. GC content) to predict the taxonomic classification of the input MAGs. The tool may also report the Average Amino Acid Identity (AAI) between input MAGs and taxonomically identified GVs.

TIGTOG is intended to be used on viral MAGs. We suggest that users first confirm that their MAGs are viral using tools such as ViralRecall (https://github.com/faylward/viralrecall). TIGTOG will issue warnings if the input sequence has an unusually large genome size and/or low coding density, which are typically observed in non-viral genomes.

## How to use
### Clone the Repository
```bash
git clone https://github.com/anhd-ha/TIGTOG.git
```
### Requirements
* All required dependencies are specified in the file `tigtog.yml` in this repository. Please download this file and use `conda` (see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)) to recreate the tigtog environment by
```shell
conda env create -f tigtog.yml
```
* If you are using conda on MacOS, please use the file `tigtog_macos.yml` to install the environment.
  
* Input for this tool is a directory that contains single contigs or MAGs as nucleic acid sequences (fasta format).


## Running 
Activate the tigtog environment by
```shell
conda activate tigtog
```
### MINIMAL USAGE: 
```shell
python tigtog.py -i <directory of fasta files> -n <project_name>
```

### Options

**-i, --input**
Folder containing input genome FASTA file. 

**-n, --project**
The project name for outputs. Please ensure it differs from the input folder name.

**-t, --cpus**
Number of CPUs to use for the HMMER3 search and AAI calculation (if requested).

**-e, --evalue**
e-value that is passed to HMMER3 for the GVOG hmmsearch (default=1e-10, recommended).

**-a, --aai**
If you would like to calculate one-way Average Amino Acid Identity (AAI) between your MAGs and a database of reference giant virus sequences (containing representatives of every genus as delineated in [Aylward et al. 2021](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001430)), you can use the -a flag.

**-r, --redo**
If you have already run the script and you want to re-run it with different parameters, you can use the -r flag to avoid re-running prodigal and HMMER.


### Examples
To get the taxonomic prediction of MAGs using 8 threads
```shell
python tigtog.py -i test_MAGs -n test_run -t 8
```
To get taxonomic prediction and run AAI calculation without re-running previous prodigal and HMMER search
```shell
python tigtog.py -i test_MAGs -n test_run -t 8 -a -r
```

## Result
* Taxonomic prediction of input MAGs (provided at the Order and Family levels) and confidence of predictions are stored in the file `<project_name>.prediction_result.tsv`.

* If you requested AAI calculation using the -a flag, the reference giant virus with the best AAI to your MAGs will also be reported in the result file, along with the AAI value, alignment fraction (AF), and the taxonomic assignment of that reference. Only AAI hits with an AF>20 will be reported.


### Contact
Please contact us if you have any questions: 
Anh Ha: anhdha@vt.edu
Frank Aylward: faylward@vt.edu


### Citation
Coming soon!


