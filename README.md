# Enterovirus D68 VP1 Nextstrain Analysis

Performs a full Nextstrain analysis on the VP1 segment of Enterovirus D68.

## Quickstart 
### Setup
You will need to install local BLAST for this analysis to work. Do this with:
`sudo apt-get install ncbi-blast+`
_(is this right?)_

Download in _tab delimited format_ all samples that are Enterovirus -> Enterovirus D -> Enterovirus D68 using ViPR's search function, without restriction on sequence length or dates. (There should be over 3,000.) 

Place this file in `enterovirus_vp1/data` (you may need to create `data`), and include the name of that file on line 9 of the Snakefile (replacing `data/allEntero-20Nov18.tsv` or similar). 

Place sequences and metadata from Swedish sequences in the `data` folder, and ensure the filenames on lines 11 and 12 of the Snakefile match your own. 

If you have other sequences & metadata (not on GenBank) you'd like to add, you can include these on lines 15 and 16.

### Regions
This script will allow you to look at sequences by region as well as country. As using VP1 means many countries are included, regions can be easier to look at for general patterns. The Snakefile is already set up for this kind of analysis, and region will be automatically generated for all downloaded sequences.

*However*, you should ensure the Swedish metadata file, and any additional 'manual' files, have an additional column called 'region' with an entry for each sample ('europe' - lowercase). Otherwise, no Swedish/manual sequences will have a region. 

### Running
Navigate to the `enterovirus_vp1` folder and run `snakemake`. Initial runs will take some time, as downloading all sequences from GenBank is slow.

All accession numbers are compared, so a sequence already included in 'Swedish' or 'manual' files will not be downloaded from GenBank.

## Reruns
This Snakefile is written to make adding new data from ViPR easier. Simply download the latest full collection of samples from ViPR (using the same instructions as above), place the new file in `data`, and replace the filename on line 9 of the Snakefile. Run `snakemake`, and the script should automatically only download and BLAST sequences with accesssion numbers that have not previously been checked (even if they were not included in the analysis). 

After adding any new sequences, the a new full Nextstrain analysis will proceed. 


# Technical Notes
## Strain names
In ViPR downloads as specified above, `strain` is not a unique identifier, as multiple segments may come from the same `strain`. This causes problems unique to VP1 analysis (with full-genome, this is not an issue). To handle this, the `vipr_parse.py` script generates new `strain` identifiers by combiing the original `strain` column with the accession number, separated by a double-underscore. 

## Blasting
ViPR sequences are not reliably labelled with the segment(s) they include (excepting whole-genome, it seems). In order to decide which sequences contain VP1, this script creates a local BLAST database against an EV-D68 reference genome VP1 sequence, then BLASTs all downloaded sequences against it.  

Sequences with matches of at least 700bp are included. This was chosen because in initial runs, using >=600bp added only 47 sequences more and >=800bp lost 289 sequences. Only the matching sequence segment is taken for analysis.

## Reruns
This Snakefile saves a copy of the most recently run parsed, downloaded ViPR file, and uses this to decide whether an accession number is 'new.' This means that even sequences that do not include VP1 (and are not included in later stages of the analysis) will not be checked again. 

Here are some graphics showing the snakemake run diagram for a new run of the pipeline, and a re-run, adding new sequences:

New run:
![alt text][run1]

Re-run:
![alt text][run2]

[run1]: https://github.com/emmahodcroft/enterovirus_vp1/tree/master/images/newrun-snakemake-DAG.PNG "new run"
[run2]: https://github.com/emmahodcroft/enterovirus_vp1/tree/master/images/rerun-snakemake-DAG.PNG "rerun"
