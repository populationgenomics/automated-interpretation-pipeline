"""
runs between classification and publishing results
takes a number of inputs:
    - Classified VCF
    - JSON describing the present compound-het pairs
    - PanelApp data

reads in all compound het pairs
reads in all panelapp details
for each variant in each participant, check MOI in affected
participants relative to the MOI described in PanelApp
"""


import json
import logging
from argparse import ArgumentParser
from collections import defaultdict
from typing import Any, Dict, List, Union

from cloudpathlib import AnyPath
from cyvcf2 import VCFReader
from peddy.peddy import Ped

from reanalysis.moi_tests import MOIRunner
from reanalysis.utils import (
    AbstractVariant,
    canonical_contigs_from_vcf,
    gather_gene_dict_from_contig,
    get_simple_moi,
    read_json_from_path,
    CompHetDict,
    CustomEncoder,
    ReportedVariant,
)


def set_up_inheritance_filters(
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    config: Dict[str, Any],
    pedigree: Ped,
    comp_het_lookup: CompHetDict,
) -> Dict[str, MOIRunner]:
    """
    parse the panelapp data, and find all MOIs in this dataset
    for each unique MOI, set up a MOI filter instance
    save each one to a dictionary

    {MOI_string: MOI_runner (with a .run() method)}

    The MOI_runner uses a MOI string to select appropriate filters

    All logic regarding how MOI is applied, and which MOIs to
    apply to which PanelApp MOI descriptions is partitioned off into
    the MOI classes. All we need here is a Run() method, that returns
    either a list of results, or an empty list

    for every variant, we can then do a simple lookup using this
    dictionary to find the correct MOI runner, and run it
    that will return all matching MOIs for the variant

    This dictionary format means we only have to set up each once
    A billion variants, 6 MOI = 6 test instances, each created once
    :param panelapp_data:
    :param config:
    :param pedigree:
    :param comp_het_lookup:
    :return:
    """

    moi_dictionary = {}

    # iterate over all genes
    for key, gene_data in panelapp_data.items():

        # skip over the stored metadata
        if '_version' in key:
            continue

        # extract the per-gene MOI, and SIMPLIFY
        gene_moi = get_simple_moi(gene_data.get('moi'))

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:

            # get a MOIRunner with the relevant filters
            moi_dictionary[gene_moi] = MOIRunner(
                pedigree=pedigree,
                target_moi=gene_moi,
                config=config['moi_tests'],
                comp_het_lookup=comp_het_lookup,
            )

    return moi_dictionary


def apply_moi_to_variants(
    variant_dict: Dict[str, Dict[str, AbstractVariant]],
    moi_lookup: Dict[str, MOIRunner],
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
) -> List[ReportedVariant]:
    """
    take a collection of all variants on a given contig & MOI filters
    find all variants/compound hets which fit the PanelApp MOI

    :param variant_dict:
    :param moi_lookup:
    :param panelapp_data:
    :return:
    """

    results = []

    # NOTE - if we want details from both sides of a compound het in the
    # MOI check or report details, this lookup gives us a way to access that data
    # For now just iterate over the individual variants
    for gene, variants in variant_dict.items():

        # extract the panel data specific to this gene
        # extract once per gene, not once per variant
        panel_gene_data = panelapp_data.get(gene)

        # variant appears to be in a red gene
        if panel_gene_data is None:
            logging.error(f'How did this gene creep in? {gene}')
            continue

        simple_moi = get_simple_moi(panel_gene_data.get('moi'))

        for variant in variants.values():

            # if this variant is category 1, 2, 3, or 4; evaluate is as a 'primary'
            if variant.category_non_support:

                # this variant is a candidate for MOI checks
                # - find the simplified MOI string
                # - use to get appropriate MOI model
                # - run variant, append relevant classification(s) to the results
                results.extend(
                    moi_lookup[simple_moi].run(
                        principal_var=variant, gene_lookup=variants
                    )
                )

    return results


def clean_initial_results(
    result_list: List[ReportedVariant], samples: List[str]
) -> Dict[str, Dict[str, ReportedVariant]]:
    """
    Possibility 1 variant can be classified multiple ways
    This cleans those to unique for final report
    Join all possible classes for the condensed variants
    :param result_list:
    :param samples: all samples from the VCF
    """

    clean_results = defaultdict(dict)

    for each_event in result_list:
        support_id = (
            ','.join(sorted(each_event.support_vars))
            if each_event.support_vars is not None
            else 'Unsupported'
        )
        var_uid = (
            f'{each_event.var_data.coords.string_format}__'
            f'{each_event.gene}__'
            f'{support_id}'
        )

        # if this variant was already found, combine the selection 'reasons'
        if var_uid in clean_results[each_event.sample]:
            # combine any possible reasons
            clean_results[each_event.sample][var_uid].reasons.update(each_event.reasons)

        # otherwise insert this variant into the dict
        else:
            clean_results[each_event.sample][var_uid] = each_event

    # Empty list for 0 variant samples; explicitly record samples in this joint-call
    # the PED could have more samples than a joint call, due to sub-setting or QC.
    # When presenting results, we want all samples with negative findings, without
    # comparing both VCF and PED files
    for sample in samples:
        if sample not in clean_results:
            clean_results[sample] = {}

    return clean_results


def main(
    labelled_vcf: str,
    comp_het: str,
    config_path: Union[str, Dict[str, Any]],
    out_json: str,
    panelapp: str,
    pedigree: str,
):
    """
    VCFs used here should be small
    These have been pre-filtered to retain only a small number of candidate variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call

    Might be able to use a single output path, just altering the extension
    Depends on how this is handled by Hail, as the object paths are Resource File paths

    Re-working of the comp-het logic means that we only store pairings as strings
    Not needing to reach the annotations attached to variant pairs opens up choices:
        - process each variant in turn (original design)
        - parse each chromosome separately, then process the group of variants
        - parse all variants, then process as a group

    these come with incrementing memory footprints...

    preference is for #2; process an entire contig together (note, still heavily
    filtered, so low variant numbers expected)
        - we can look-up the partner variant's attributes in future if we want
        - this will be required when we do familial checks, e.g. for a compound het,
            we need to check the presence/absence of a pair of variants in unaffected
            family members

    :param labelled_vcf:
    :param comp_het:
    :param config_path:
    :param out_json:
    :param panelapp:
    :param pedigree:
    """

    # parse the pedigree from the file (via write to temp)
    with open('i_am_a_temporary.ped', 'w', encoding='utf-8') as handle:
        handle.write(AnyPath(pedigree).read_text())
    pedigree_digest = Ped('i_am_a_temporary.ped')

    # parse panelapp data from dict
    panelapp_data = read_json_from_path(panelapp)

    # get the runtime configuration
    if isinstance(config_path, dict):
        config_dict = config_path
    elif isinstance(config_path, str):
        config_dict = read_json_from_path(config_path)
    else:
        raise Exception(
            f'What is the conf path then?? "{config_path}": {type(config_path)}'
        )

    # set up the inheritance checks
    moi_lookup = set_up_inheritance_filters(
        panelapp_data=panelapp_data,
        pedigree=pedigree_digest,
        config=config_dict,
        comp_het_lookup=read_json_from_path(comp_het),
    )

    # open the VCF using a cyvcf2 reader
    vcf_opened = VCFReader(labelled_vcf)

    # permit a blacklist to exclude known artefacts
    variant_blacklist = None
    if 'variant_blacklist' in config_dict:
        variant_blacklist = read_json_from_path(config_dict['variant_blacklist'])

    results = []

    # obtain a set of all contigs with variants
    for contig in canonical_contigs_from_vcf(vcf_opened):

        # assemble {gene: [var1, var2, ..]}
        contig_dict = gather_gene_dict_from_contig(
            contig=contig,
            variant_source=vcf_opened,
            config=config_dict,
            panelapp_data=panelapp_data,
            blacklist=variant_blacklist,
        )
        results.extend(
            apply_moi_to_variants(
                variant_dict=contig_dict,
                moi_lookup=moi_lookup,
                panelapp_data=panelapp_data,
            )
        )

    # remove duplicate variants
    cleaned_results = clean_initial_results(results, samples=vcf_opened.samples)

    # dump results using the custom-encoder to transform sets & DataClasses
    with AnyPath(out_json).open('w') as fh:
        json.dump(cleaned_results, fh, cls=CustomEncoder, indent=4)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument(
        '--labelled_vcf', help='Path to VCF resulting from category labelling process'
    )
    parser.add_argument(
        '--comp_het', help='JSON file containing all compound het variants'
    )
    parser.add_argument('--config_path', help='path to the runtime JSON config')
    parser.add_argument('--pedigree', help='Path to joint-call PED file')
    parser.add_argument('--panelapp', help='Path to JSON file of PanelApp data')
    parser.add_argument('--out_json', help='Path to write JSON results to')
    args = parser.parse_args()
    main(
        labelled_vcf=args.labelled_vcf,
        comp_het=args.comp_het,
        config_path=args.config_path,
        out_json=args.out_json,
        panelapp=args.panelapp,
        pedigree=args.pedigree,
    )
