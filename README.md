# soYouWantTheLongestTranscript

This tool prepares genomes/transcriptomes for tasks such as orthology prediction, where multiple isoforms per gene are undesirable. It ingests a fasta file (a genome/transcriptome in nucleotide space) and a gff file, and filters the genome/transcriptome down to the longest transcript for each gene. It supports processing multiple genomes in parallel, and in the absence of a parent attribute, it supports deduplication of isoforms using regex on a hierarchical id attribute.

# regex examples

[simple example](https://regex101.com/r/bXMjIj/2)

[multi-part example](https://regex101.com/r/F561kR/4)

# usage:

```
usage: longest_transcript.py [-h] [-i IN_DIR] [-o OUT_DIR]
                 [-s SAMPLE_NAMES [SAMPLE_NAMES ...]] [-f FASTA_EXT]
                 [-g GFF_EXT] [-r REGEX] [--freeze_attribute_order]

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
  --freeze_attribute_order
                        Due to implementation, a slight speed boost can be
                        obtained if you don't care about maintaining attribute
                        order. This option exists for regression testing, or
                        for workflows for which this order matters.

```
