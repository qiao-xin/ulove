# ULOVE v1.0.0
A software used to assess gene annotation completeness with universal low-copy orthologs Viridiplantae (ULOVE).

| | |
| --- | --- |
| Authors | Xin Qiao ([Xin Qiao](https://github.com/qiao-xin)) |
| | Qionghou Li ([Qionghou Li](https://github.com/LQHHHHH)) |
| Email   | <qiaoxin@njau.edu.cn> |

## Contents
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Running](#running)
* [Result Files](#result-files)
* [Citation](#citation)

## Dependencies

- [Perl](https://www.perl.org)
- [HMMER](http://hmmer.org)
- [Python](https://www.python.org)
- [Matplotlib](https://matplotlib.org)

## Installation

```bash
cd /path/to/your_software_directory  # or any directory of your choice
git clone https://github.com/qiao-xin/ulove.git
cd ulove
chmod 750 ulove.pl generate_figure.py set_PATH.sh decompress_hmm_files.sh
sh decompress_files.sh
source set_PATH.sh
```

Test you can run ULOVE:
```bash
ulove.pl
```
If install successfully, ULOVE will print its 'help' text.
```
Usage: ulove.pl -i input_directory -o output_directory -l lineage_dataset -s species_code
#####################
-i the directory storing protein sequence file, please name the file like this Arath.pep (*species_code*.pep)
-o the directory storing output data
-l the lineages dataset used for completeness assessement, please choose an appropriate lineage for your species:
    viridiplantae
    chlorophyta
    streptophyta
    embryophyta
    tracheophyta
    spermatophyta
    angiosperms
    monocots
    eudicots
-s species code, for example, Arath can be used as the species code of Arabidopsis thaliana
```

## Running

A typical command to run ULOVE could look like this:
```bash
ulove.pl -i ./data -o . -l viridiplantae -s Arath
```

## Result Files
### 1 - short_summary.specific.viridiplantae.Arath.ulove.txt
```
# ULOVE version is: 1.0.0
# The lineage dataset is: viridiplantae (Creation date: 2025-06-06, number of genomes: 1163, number of ULOVEs: 371)
# Summarized benchmarking in ULOVE notation for file ./Arath.pep
# ULOVE was run in mode: protein

	***** Results: *****

	C:99.5%[S:98.7%,D:0.8%],F:0.5%,M:0.0%,n:371
	369	Complete ULOVEs (C)
	366	Complete and single-copy ULOVEs (S)
	3	Complete and duplicated ULOVEs (D)
	2	Fragmented BUSCOs (F)
	0	Missing ULOVEs (M)
	371	Total ULOVE groups searched

Dependencies and versions:
	hmmsearch: 3.3.2
	ulove: 1.0.0
```

### 2 - Arath_ulove_assessment_result.pdf
![Arath_ulove_assessment_result](/data/Arath_ulove_assessment_result.png)
*Figure 1: The visualization of gene annotation completeness assessment of Arabidopsis thaliana.*

### 3 - A directory "hmmsearch_results", including output files from hmmsearch
```
OG0004472_Arath.domtblout
OG0004472_Arath.out
OG0004571_Arath.domtblout
OG0004571_Arath.out
OG0004791_Arath.domtblout
OG0004791_Arath.out
...
```

## Citation
Preparation and submission.
