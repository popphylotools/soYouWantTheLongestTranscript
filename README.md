# soYouWantTheLongestTranscript

This tool prepares genomes/transcriptomes for tasks such as orthology prediction, where multiple isoforms per gene are undesireable. It ingests a fasta file (a genome/transcriptome in nucleotide space) and a gff file, and filters the genome/transcriptome down to the longest transcript for each gene. It can process multiple data sources at once, and can accommodate variable formatting of the id field in the gff file using regex (examples provided below for common formats). 

# help:

```
usage: sywtlt.py [-h] [-i IN_DIR] [-o OUT_DIR]
                 [-s SAMPLE_NAMES [SAMPLE_NAMES ...]] [-f FASTA_EXT]
                 [-g GFF_EXT] [-r REGEX]

Returns the longest transcript per gene in both nucleotide space and peptide
space along with a simplified gff with only the relevant lines.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_DIR, --in_dir IN_DIR
                        Path to directory with input gff and fasta.
  -o OUT_DIR, --out_dir OUT_DIR
                        Path to output directory (created if missing, unique
                        created if non-empty).
  -s SAMPLE_NAMES [SAMPLE_NAMES ...], --sample_names SAMPLE_NAMES [SAMPLE_NAMES ...]
                        One or more space separated sample_name(s).
  -f FASTA_EXT, --fasta_ext FASTA_EXT
                        Optionally set non-default (.fna) file extension.
  -g GFF_EXT, --gff_ext GFF_EXT
                        Optionally set non-default (.gff) file extension.
  -r REGEX, --regex REGEX
                        Optionally group isoforms on id field using regex. By
                        default parent field is used to find competing
                        isoforms. Quote this parameter to avoid problems with
                        special characters. For help with regex, visit:
                        https://regex101.com/r/F561kR/4

```
