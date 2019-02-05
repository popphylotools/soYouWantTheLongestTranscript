#!/usr/bin/env python
# coding=utf-8
""""""

import argparse
import logging
import multiprocessing as mp
import os
import shutil
import sys
import re

import Bio
import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta


# logging ##############################################################################################################

class OneLineExceptionFormatter(logging.Formatter):
    """custom logging configuration for exceptions to
    format on a single line, but with <newline> markers
    for easy reformating in a text editor"""
    def formatException(self, exc_info):
        result = super().formatException(exc_info)
        return repr(result)

    def format(self, record):
        """Format error message to fit on a single line.
        """
        result = super().format(record)
        if record.exc_text:
            result = result.replace(r"\n", "<newline>")
        return result


# initialize logging
handler = logging.StreamHandler()
formatter = OneLineExceptionFormatter(logging.BASIC_FORMAT)
handler.setFormatter(formatter)
root = logging.getLogger()
root.setLevel(os.environ.get("LOGLEVEL", "DEBUG"))
root.addHandler(handler)
log = logging.getLogger("extractor.py")


# directory management #################################################################################################

def create_unique_dir(path, limit=99):
    """Return a path to an empty directory. Either the dir at path, or a dir of the form 'path + _01'
    :param path: The initial path to use
    :param limit: The maximum number of directory variations this function will attempt to create.
    :return: A path to an empty directory.
    """
    if not os.path.isdir(path):
        log.debug("Creating output directory - {}".format(path))
        os.makedirs(path)
    else:
        width = len(str(limit))
        original = path.rstrip(os.sep)
        if len(os.listdir(original)) == 0:
            log.debug("Using output directory - {}".format(path))
            return original  # folder empty, let's use it
        count = 1
        while count < limit:
            try:
                os.mkdir(path)
                log.debug("Creating output directory - {}".format(path))
                return path
            except OSError as path_error:
                if path_error.errno == 17:  # file exists
                    path = "{0}_{1:0>{2}}".format(original, count, width)
                    count += 1
                else:
                    raise
        else:
            msg = "could not uniquely create directory {0}: limit `{1}` reached"
            raise Exception(msg.format(original, limit))


# helper functions #####################################################################################################

def group_isos_by(iso, regex):
    if regex is None:
        key = iso.attributes["Parent"][0]
        return key
    else:
        try:
            text = iso.attributes["ID"][0]
            search_obj = re.search(regex, text)
            key = search_obj.groups()
            return key
        except Exception as e:
            log.exception("regex failed: {}".format(e))


def rank_isos_by(iso, db):
    iso_len = db.children_bp(iso, child_featuretype='CDS', merge=False, ignore_strand=False)
    return -iso_len, iso.attributes["ID"]  # biggest first, alphabetical for ties


def get_frame(iso):
    if iso.frame == ".":
        return 0
    else:
        return int(iso.frame)


def get_db(in_gff_fn):
    db_fn = in_gff_fn + ".db"
    if os.path.isfile(db_fn):
        try:
            db = gffutils.FeatureDB(db_fn, keep_order=True)
            return db
        except TypeError as e:
            message = "Error reading db, creating {}.bak and recreating - db_fn:{} - error:{}".format(db_fn, db_fn, e)
            log.debug(message)
            shutil.copy2(db_fn, db_fn + '.bak')

    db = gffutils.create_db(data=in_gff_fn,
                            dbfn=db_fn,
                            force=True,
                            merge_strategy='merge',
                            id_spec=['ID', 'Name'])

    return db


# core logic ###########################################################################################################

def sample_input_generator(args):
    group_regex = args.regex
    out_dir = create_unique_dir(args.out_dir)
    for sample_name in args.sample_names:

        in_gff_fn = os.path.join(args.in_dir, sample_name + args.gff_ext)
        in_fasta_fn = os.path.join(args.in_dir, sample_name + args.fasta_ext)

        yield sample_name, in_gff_fn, in_fasta_fn, out_dir, group_regex


def process_sample(in_tup):
    sample_name, in_gff_fn, in_fasta_fn, out_dir, group_regex = in_tup

    db = get_db(in_gff_fn)

    fasta = Fasta(in_fasta_fn)

    selected_isos = select_isoforms(sample_name, db, group_regex)

    extract_and_write_sequences(sample_name, out_dir, selected_isos, db, fasta)

    return sample_name


def select_isoforms(sample_name, db, group_regex):
    iso_groups = dict()
    selected_isos = dict()

    def get_cds_parents(db, level=1):
        cds_list = list(db.features_of_type(featuretype='CDS'))
        parents = set()
        for cds in cds_list:
            parents.update(db.parents(cds, level=level))
        parents = sorted(list(parents), key=lambda p: (p.seqid, p.start, p.id))
        return parents

    # group with group_isos_by function
    for iso in get_cds_parents(db):
        try:
            group_key = group_isos_by(iso, group_regex)
        except Exception as e:
            message = "problem generating iso group key - sample_name:{} - id:{} - error:{}".format(sample_name, iso.attributes["ID"], e)
            log.debug(message)
            continue
        if group_key not in iso_groups:
            iso_groups[group_key] = []
        iso_groups[group_key].append(iso)

    # sort groups with rank_isos_by function and select first from each group
    for iso_list in iso_groups.values():
        iso = sorted(iso_list, key=lambda iso: rank_isos_by(iso, db))[0]
        selected_isos[iso.attributes["ID"][0]] = iso

    return selected_isos


def extract_and_write_sequences(sample_name, out_dir, selected_isos, db, fasta):
    pre_filter_count = len(selected_isos)
    post_filter_count = 0

    sorted_key_list = sorted(selected_isos.keys(), key=lambda k: (selected_isos[k].seqid, selected_isos[k].start))

    with open(os.path.join(out_dir, sample_name + ".fna"), 'w') as nuc_f, \
         open(os.path.join(out_dir, sample_name + ".faa"), 'w') as pep_f, \
         open(os.path.join(out_dir, sample_name + ".gff"), 'w') as gff_f:

        for key in sorted_key_list:
            # remove exons filtering errors and empties
            cds_list = list(db.children(selected_isos[key], order_by="start", featuretype="CDS"))
            try:
                n_seq = Seq("".join([c.sequence(fasta, use_strand=False) for c in cds_list]))
            except KeyError as e:
                message = "KeyError - sp:{} - key:{} - error:{}".format(sample_name, key, e)
                log.info(message)
                continue
            except ValueError as e:
                message = "ValueError - sp:{} - key:{} - error:{}".format(sample_name, key, e)
                log.info(message)
                continue
            if n_seq is "":
                message = "Excluding - sp:{} - key:{} - n_seq is empty".format(sample_name, key)
                log.info(message)
                continue
            # take reverse complement if appropriate
            if selected_isos[key].strand == "-":
                n_seq = n_seq.reverse_complement()
                frame = get_frame(cds_list[-1])
            else:
                frame = get_frame(cds_list[0])

            # trim front according to frame and pad last frame with N if needed
            trans_n_seq = n_seq[frame:]
            pad = (3 - (len(trans_n_seq) % 3)) % 3
            if pad != 0:
                trans_n_seq += "N" * pad

            # translate to peptide sequence filtering errors, internal stops, and empties
            try:
                p_seq = trans_n_seq.translate()
            except Bio.Data.CodonTable.TranslationError as e:
                message = "TranslationError - sp:{} - key:{} - error:{}".format(sample_name, key, e)
                log.info(message)
                continue
            if "*" in p_seq[:-1]:
                message = "Excluding - sp:{} - key:{} - internal stop codon".format(sample_name, key)
                log.info(message)
                continue
            if len(p_seq) == 0:
                message = "Excluding - sp:{} - key:{} - p_seq is empty".format(sample_name, key)
                log.info(message)
                continue

            post_filter_count += 1

            # write nuc_fasta and pep_fasta lines
            nuc_f.write(SeqRecord(n_seq, id=key, description="").format("fasta"))
            pep_f.write(SeqRecord(p_seq, id=key, description="").format("fasta"))

            # write gff lines
            gff_group = list(db.parents(selected_isos[key], level=1)) + [selected_isos[key]] + cds_list
            for gff_line in gff_group:
                gff_f.write(str(gff_line) + "\n")

    message = "sp: {} - isoforms kept: {} - isoforms discarded: {} - percent discarded: {}"
    log.info(message.format(sample_name, post_filter_count, pre_filter_count - post_filter_count,
                            (pre_filter_count - post_filter_count) / pre_filter_count * 100))

    return sample_name, post_filter_count, pre_filter_count


# cli and multiprocessing ##############################################################################################

class HelpAndQuitOnFailParser(argparse.ArgumentParser):
    """custom argparse configuration
    if error parsing, prints help and exits"""
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    parser = HelpAndQuitOnFailParser(description=("returns the longest transcript per gene "
                                                  "in both nucleotide space and peptide space "
                                                  "along with a simplified gff with only the relevant lines."))

    # directories
    parser.add_argument('-i', '--in_dir', type=str, default="input",
                        help='path to directory with input gff and fasta')

    parser.add_argument('-o', '--out_dir', type=str, default="output",
                        help='path to output directory (created if missing, unique created if non-empty)')

    # sample names
    parser.add_argument('-s', '--sample_names', nargs='+', type=str, default=['GCF_000789215.1_ASM78921v2_genomic'],
                        help='one or more space separated sample_name(s)')

    # extensions
    parser.add_argument('-f', '--fasta_ext', type=str, default=".fna",
                        help="optionally set non-default (.fna) file extension")

    parser.add_argument('-g', '--gff_ext', type=str, default=".gff",
                        help='optionally set non-default (.gff) file extension')

    # regex
    parser.add_argument('-r', '--regex', type=str, default=None,
                        help='optionally set isoform grouping on id field (defaults to parent field)')

    args = parser.parse_args()

    processor_count = min(len(args.sample_names), mp.cpu_count())
    with mp.Pool(processor_count) as p:
        data = p.map(process_sample, sample_input_generator(args))

    print(data)


if __name__ == "__main__":
    # run main
    try:
        exit(main())
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
