# soYouWantTheLongestTranscript

This tool prepares genomes/transcriptomes for tasks such as orthology prediction, where multiple isoforms per gene are undesireable. It ingests a fasta file (a genome/transcriptome in nucleotide space) and a gff file, and filters the genome/transcriptome down to the longest transcript for each gene. It can process multiple data sources at once, and can accommodate variable formatting of the id field in the gff file using regex (examples provided below for common formats). 

# help:

```
usage: sywtlt.py [-h] [-i IN_DIR] [-o OUT_DIR]
                 [-s SAMPLE_NAMES [SAMPLE_NAMES ...]] [-f FASTA_EXT]
                 [-g GFF_EXT] [-r REGEX]

returns the longest transcript per gene in both nucleotide space and peptide
space along with a simplified gff with only the relevant lines.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_DIR, --in_dir IN_DIR
                        path to directory with input gff and fasta
  -o OUT_DIR, --out_dir OUT_DIR
                        path to output directory (created if missing, unique
                        created if non-empty)
  -s SAMPLE_NAMES [SAMPLE_NAMES ...], --sample_names SAMPLE_NAMES [SAMPLE_NAMES ...]
                        one or more space separated sample_name(s)
  -f FASTA_EXT, --fasta_ext FASTA_EXT
                        optionally set non-default (.fna) file extension
  -g GFF_EXT, --gff_ext GFF_EXT
                        optionally set non-default (.gff) file extension
  -r REGEX, --regex REGEX
                        optionally group isoforms on id field using regex (by
                        default finds competing isoforms using parent field)

```

-----

when using regex, quote the string to avoide problems with special characters:

```
python regex_test.py -r "^(.+?)\|(.+?)_(.+?)_"
```

for more information on how to use regex, check out this [example](https://regex101.com/r/F561kR/4)
