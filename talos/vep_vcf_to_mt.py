#!/usr/bin/env python3

"""
This is an adapter process to take a VCF annoated by VEP, and re-arrange/re-name/re-type the attributes into
the same format as generated by the Broad seqr-loader pipeline

We don't need all the VEP fields, and the ones that we do need have different names depending on whether they
are generated by VEP's JSON output, or the direct VCF annotation
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import Any

import hail as hl

from talos.static_values import get_logger

MISSING_INT = hl.int32(0)

# get the hail types - this will be used in re-coding the attributes
HAIL_TYPES = {
    'int': {'type': hl.int32, 'default': 0},
    'string': {'type': hl.str, 'default': ''},
    'float': {'type': hl.float64, 'default': 0.0},
}
RELEVANT_FIELDS = relevant_fields = [
    'am_class',
    'am_pathogenicity',
    'biotype',
    'cdna_position',
    'cds_position',
    'consequence',
    'ensp',
    'exon',
    'feature',
    'feature_type',
    'gene',
    'gnomade_af',
    'gnomadg_af',
    'hgvsc',
    'hgvsp',
    'lof',
    'mane_select',
    'mane_plus_clinical',
    'polyphen',
    'protein_position',
    'sift',
    'symbol',
    'variant_class',
]
TYPE_UPDATES: dict[str, dict] = {
    'am_pathogenicity': {'insert': False, 'type': 'float'},
    'am_class': {'insert': False, 'type': 'string'},
    'biotype': {'type': 'string'},
    'cdna_position': {'type': 'int'},
    'cds_position': {'type': 'int'},
    'consequence': {'insert': False, 'name': 'consequence_terms'},
    'gnomade_af': {'type': 'float'},
    'gnomadg_af': {'type': 'float'},
    'sift': {'insert': False, 'type': 'string', 'name': 'sift_prediction'},
    'sift_prediction': {'type': 'string'},
    'polyphen': {'insert': False, 'name': 'polyphen_prediction'},
    'polyphen_prediction': {'type': 'string'},
    'protein_position': {'type': 'int'},
}


def csq_strings_into_hail_structs(csq_strings: list[str], mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Take the list of CSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        mt (hl.MatrixTable): the MatrixTable to annotate

    Returns:
        a MatrixTable with the VEP annotations re-arranged
    """

    # get the CSQ contents as a list of lists of strings, per variant
    split_csqs = mt.info.CSQ.map(lambda csq_entry: csq_entry.split('\|'))  # type: ignore

    # generate a struct limited to the fields of interest
    # use a couple of accessory methods to re-map the names and types for compatibility
    return mt.annotate_rows(
        vep=hl.struct(
            transcript_consequences=split_csqs.map(
                lambda x: hl.struct(
                    **{
                        remap_name(csq_strings[n]): remap_type(csq_strings[n], x[n])
                        for n in range(len(csq_strings))
                        if csq_strings[n] in RELEVANT_FIELDS
                    },
                ),
            ),
        ),
    )


def remap_name(input_name: str) -> str:
    """
    take the current variable name, consider re-coding it

    Args:
        input_name (str): the attribute name in the VEP data

    Returns:
        a new name, or the original if satisfactory
    """

    if input_name not in TYPE_UPDATES:
        return input_name
    return TYPE_UPDATES[input_name].get('name', input_name)


def remap_type(input_name, input_value) -> Any:
    """
    take the field name & value, re-type them based on a lookup

    Args:
        input_name ():
        input_value ():

    Returns:
        the same value, with an appropriate Type
    """
    # special case for consequence_terms - take the `&` delimited String and make it a list/ArrayExpression
    if input_name == 'consequence':
        return input_value.split('&')

    # check if we need to re-map it
    elif input_name not in TYPE_UPDATES:
        return input_value

    # if the value needs to be re-typed, cast it
    if re_type := TYPE_UPDATES[input_name].get('type', False):
        # when the current contents are a missing String, run a type cast on 0 instead
        # this is on the basis that the values default to String, so the only translations are to numeric
        return hl.if_else(
            input_value != '',
            HAIL_TYPES[re_type]['type'](input_value),
            HAIL_TYPES[re_type]['type'](HAIL_TYPES[re_type]['default']),
        )

    return input_value


def extract_and_split_csq_string(vcf_path: str) -> list[str]:
    """
    Extract the CSQ header from the VCF and split it into a list of strings
    Args:
        vcf_path (str): path to the local VCF

    Returns:
        list of strings
    """

    # get the headers from the VCF
    all_headers = hl.get_vcf_metadata(vcf_path)

    # get the '|'-delimited String of all header names
    csq_whole_string = all_headers['info']['CSQ']['Description'].split('Format: ')[-1]

    # split it all on pipes, return the list
    return csq_whole_string.lower().split('|')


def implant_detailed_af(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    VEP by default doesn't provide the AF granularity Talos is built on
    so substitute the general gnomad exome/genome AF as all separate categories
    and something similar for SpliceAI, which we are getting from a separate source

    Args:
        mt ():

    Returns:

    """
    return mt.annotate_rows(
        gnomad_exomes=hl.struct(
            AF=mt.vep.transcript_consequences[0].gnomade_af,
            AN=MISSING_INT,
            AC=MISSING_INT,
            Hom=MISSING_INT,
            Hemi=MISSING_INT,
        ),
        gnomad_genomes=hl.struct(
            AF=mt.vep.transcript_consequences[0].gnomadg_af,
            AN=MISSING_INT,
            AC=MISSING_INT,
            Hom=MISSING_INT,
            Hemi=MISSING_INT,
        ),
    )


def insert_spliceai_annotation(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    only use if spliceAI is absent from the current MT schema

    Args:
        mt ():

    Returns:
        the same MT, but STRONGER
    """

    return mt.annotate_rows(splice_ai=hl.struct(delta_score=hl.float64(0), splice_consequence=hl.str('')))


def insert_am_annotations_if_missing(mt: hl.MatrixTable, am_table: str | None = None) -> hl.MatrixTable:
    """
    Load up a Hail Table of AlphaMissense annotations, and annotate this data unless the AM annotations already exist

    Args:
        mt ():
        am_table (str | None):
    """

    # check if am is missing
    fields_and_types = dict(mt.vep.transcript_consequences[0].items())
    if 'am_class' in fields_and_types:
        get_logger().info('am_class already present, skipping')
        return mt

    get_logger().info(f'Reading AM annotations from {am_table} and applying to MT')

    if am_table is None:
        get_logger().error('AM annotations table is not present, and AM annotations are not in the VCF. Please create')

    # read in the hail table containing alpha missense annotations
    am_ht = hl.read_table(am_table)

    # gross - this needs a conditional application based on the specific transcript_id
    mt = mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=hl.map(
                lambda x: x.annotate(
                    am_class=hl.if_else(
                        x.feature == am_ht[mt.row_key].transcript_id,
                        am_ht[mt.row_key].am_class,
                        hl.str(''),
                    ),
                    am_pathogenicity=hl.if_else(
                        x.feature == am_ht[mt.row_key].transcript_id,
                        am_ht[mt.row_key].am_pathogenicity,
                        hl.float64(0),
                    ),
                ),
                mt.vep.transcript_consequences,
            )
        )
    )

    return mt


def insert_missing_annotations(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Insert any missing annotations into the MatrixTable if they aren't present
    use all the expected field names, and include the relevant expected data type

    Args:
        mt ():

    Returns:
        the same MT, with any missing CSQ entries added
    """
    fields_and_types = dict(mt.vep.transcript_consequences[0].items())
    for fieldname, fieldtype in TYPE_UPDATES.items():
        if fieldname not in fields_and_types:
            # a couple of special cases - don't insert if the field is marked as 'insert: False'
            # these are where we renamed the native VEP fields to the Broad schema, we don't want to also
            # insert the old one with a naff default value
            # alpha missense scores are too core to the category reasoning - we annotate instead of using blank values
            if not fieldtype.get('insert', True):
                continue

            get_logger().info(f'{fieldname} was absent, inserting {fieldtype}')
            mt = mt.annotate_rows(
                vep=mt.vep.annotate(
                    transcript_consequences=hl.map(
                        lambda x: x.annotate(
                            **{
                                fieldname: HAIL_TYPES[fieldtype['type']]['type'](
                                    HAIL_TYPES[fieldtype['type']]['default']
                                )
                            }
                        ),
                        mt.vep.transcript_consequences,
                    )
                )
            )
    return mt


def main():
    """
    take an input VCF and an output MT path
    optionally, also supply the alpha_missense table created by helpers/parse_amissense_into_ht.py
    if AlphaMissense annotations aren't already present in the VCF, this will annotate them in
    """

    parser = ArgumentParser(description='Takes a VEP annotated VCF and makes it a MT')
    parser.add_argument('vcf', help='Path to the annotated VCF')
    parser.add_argument('output', help='output MatrixTable path')
    parser.add_argument('--am', help='Hail Table containing AlphaMissense annotations', default=None)
    args, unknown = parser.parse_known_args()
    if unknown:
        raise ValueError(f'Whats the deal with {unknown}?')

    # maybe this should be a larger local cluster
    # and maybe partitions should be managed/enforced
    hl.init()
    hl.default_reference('GRCh38')

    # pull and split the CSQ header line
    vep_header_elements = extract_and_split_csq_string(args.vcf)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(args.vcf, array_elements_required=False, force_bgz=True)

    # checkpoint it locally to make everything faster
    mt = mt.checkpoint('checkpoint.mt', overwrite=True, _read_if_exists=True)

    # re-shuffle the CSQ elements
    mt = csq_strings_into_hail_structs(vep_header_elements, mt)

    # insert super detailed AF structure - no reannotation, just re-organisation
    mt = implant_detailed_af(mt)

    # if we need AlphaMissense scores to be added, add them
    mt = insert_am_annotations_if_missing(mt, am_table=args.am)

    # check if all required annotations are present - insert if absent
    mt = insert_missing_annotations(mt)

    # check if spliceAI annotations are present - insert if absent
    if 'splice_ai' not in mt.row_value:
        mt = insert_spliceai_annotation(mt)

    # audit all required annotations?
    mt.describe()

    mt.write(args.output)


if __name__ == '__main__':
    main()
