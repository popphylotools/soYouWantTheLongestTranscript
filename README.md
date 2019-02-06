help:

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
                        optionally set isoform grouping on id field (defaults
                        to parent field)
```

-----

when using regex, quote the string like so:

```
python regex_test.py -r "^(.+?)::.+?::(.+?)::.+$"
```

for more information on how to use regex, check out this [example](https://regex101.com/r/0C1B0B/1)