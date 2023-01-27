#!/usr/bin/env python3


"""
Complete revision
"""


import logging
import sys

from datetime import datetime
from dateutil.relativedelta import relativedelta

import click

from cpg_utils.config import get_config

from reanalysis.utils import (
    find_latest_file,
    get_json_response,
    get_simple_moi,
    read_json_from_path,
    save_new_historic,
    write_output_json,
    IRRELEVANT_MOI,
    ORDERED_MOIS,
)


PanelData = dict[str, dict | list[dict]]
PANELAPP_HARD_CODED_DEFAULT = 'https://panelapp.agha.umccr.org/api/v1/panels'
PANELAPP_BASE = get_config()['workflow'].get('panelapp', PANELAPP_HARD_CODED_DEFAULT)

# pylint: disable=no-value-for-parameter,unnecessary-lambda


def request_panel_data(url: str) -> tuple[str, str, list]:
    """
    just takes care of the panelapp query
    Args:
        url ():

    Returns:
        components of the panelapp response
    """

    panel_json = get_json_response(url)

    # steal attributes
    panel_name = panel_json.get('name')
    panel_version = panel_json.get('version')
    panel_genes = panel_json.get('genes')

    # log name and version
    logging.info(f'Current {panel_name} version: {panel_version}')

    return panel_name, panel_version, panel_genes


def get_panel_green(
    gene_dict: PanelData,
    old_data: dict,
    panel_id: int | None = None,
    version: str | None = None,
):
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    Args:
        gene_dict (): dictionary to continue populating
        old_data ():
        panel_id (): specific panel or 'base' (e.g. 137)
        version (): version, optional. Latest panel unless stated
    """

    if panel_id is None:
        panel_id = get_config()['workflow'].get('default_panel', 137)

    panel_url = f'{PANELAPP_BASE}/{panel_id}'

    # include the version if required
    if version:
        panel_url = f'{panel_url}?version={version}'

    panel_name, panel_version, panel_genes = request_panel_data(panel_url)

    # add metadata for this panel & version
    gene_dict['metadata'].append(
        {'name': panel_name, 'version': panel_version, 'id': panel_id}
    )

    # iterate over the genes in this panel result
    for gene in panel_genes:

        # only retain green genes
        if gene['confidence_level'] != '3' or gene['entity_type'] != 'gene':
            continue

        ensg = None
        symbol = gene.get('entity_name')
        chrom = None

        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensembl_data = content[list(content.keys())[0]]
                ensg = ensembl_data['ensembl_id']
                chrom = ensembl_data['location'].split(':')[0]

        if ensg is None:
            logging.info(f'Gene "{symbol} lacks an ENSG ID, so it is being excluded')
            continue

        # check if this is a new gene in this analysis
        # all panels previously containing this gene
        gene_panels_for_this_gene = old_data['genes'].get(ensg, {}).get('panels', [])
        gene_prev_in_panel = panel_id in gene_panels_for_this_gene

        if not gene_prev_in_panel:
            gene_panels_for_this_gene.append(panel_id)

        exact_moi = gene.get('mode_of_inheritance', 'unknown').lower()

        # either update or add a new entry
        if ensg in gene_dict['genes'].keys():
            this_gene = gene_dict['genes'][ensg]

            # add this moi to the set
            if exact_moi not in IRRELEVANT_MOI:
                this_gene['moi'].add(exact_moi)

            # if this is/was new - it's new
            if not gene_prev_in_panel:
                this_gene['panels'].add(panel_id)
                this_gene['new'].append(panel_id)

        else:

            # save the entity into the final dictionary
            gene_dict['genes'][ensg] = {
                'symbol': symbol,
                'moi': {exact_moi} if exact_moi not in IRRELEVANT_MOI else set(),
                'new': [] if gene_prev_in_panel else [panel_id],
                'panels': set(map(int, gene_panels_for_this_gene)),
                'chrom': chrom,
            }


def get_best_moi(gene_dict: dict):
    """
    From the collected set of all MOIs, take the most lenient
    If set was empty (i.e. no specific MOI found) find the default

    Default is found autosome/sex chrom aware

    Args:
        gene_dict (): the 'genes' index of the collected dict
    """

    for content in gene_dict.values():

        # accept the simplest MOI if no exact moi found
        if not content['moi']:
            content['moi'] = get_simple_moi(None, chrom=content['chrom'])
            continue

        # otherwise accept the most lenient valid MOI
        moi_set = {
            get_simple_moi(moi, chrom=content['chrom']) for moi in content['moi']
        }

        # force a combined MOI here
        if 'Biallelic' in moi_set and 'Monoallelic' in moi_set:
            content['moi'] = 'Mono_And_Biallelic'

        else:
            # take the more lenient of the gene MOI options
            content['moi'] = sorted(moi_set, key=lambda x: ORDERED_MOIS.index(x))[0]


def read_panels_from_participant_file(panel_json: str) -> set[int]:
    """
    reads the per-participants panels into a set
    Args:
        panel_json (): path to a per-participant panel dump

    Returns:
        set of all the panels across all participants
    """
    participant_panels = read_json_from_path(panel_json)
    panel_set = set()
    for details in participant_panels.values():
        panel_set.update(details.get('panels', []))

    return panel_set


def find_core_panel_version() -> str | None:
    """
    take the default panel ID from config
    iterate through its associated activities
    return the panel version which was closest to 12 months ago

    or return None, i.e. if the panel is not >= 12 months old

    Returns:
        a version string from 12 months prior
    """

    date_threshold = datetime.today() - relativedelta(years=1)
    panel_id = get_config()['workflow'].get('default_panel', 137)

    # activities URL
    activities_url = f'{PANELAPP_BASE}/{panel_id}/activities'

    # query for data from this endpoint
    activities: list[dict] = get_json_response(activities_url)

    # iterate through all activities on this panel
    for activity in activities:

        # cast the activity datetime to day-resolution
        activity_date = datetime.strptime(activity['created'].split('T')[0], '%Y-%M-%d')

        # keep going until we land on the day, or skip past it
        if activity_date > date_threshold:
            continue

        return activity['panel_version']

    # it's possible we won't find one for some panels
    return None


def get_new_genes(older_version: str) -> set[str]:
    """
    query for two versions of the same panel
    find all genes new on that panel between versions
    Args:
        older_version ():

    Returns:

    """

    # I don't like this implementation

    # this should be a None'able object
    old_data = {'genes': {}}

    current: PanelData = {'metadata': [], 'genes': {}}
    get_panel_green(current, old_data)

    old: PanelData = {'metadata': [], 'genes': {}}
    get_panel_green(old, old_data, version=older_version)

    return set(current['genes'].keys()).difference(set(old['genes'].keys()))


def overwrite_new_status(gene_dict: PanelData, new_genes: set[str]):
    """
    ignores any previous notion of new, replaces it with a manually assigned one

    Args:
        gene_dict ():
        new_genes ():

    Returns:

    """

    panel_id = get_config()['workflow'].get('default_panel', 137)

    for gene, gene_data in gene_dict['genes'].items():
        if gene not in new_genes:
            gene_data['new'] = []
        else:
            gene_data['new'] = [panel_id]


@click.command()
@click.option('--panels', help='JSON of per-participant panels')
@click.option('--out_path', required=True, help='destination for results')
@click.option('--previous', help="previous data for finding 'new' genes")
def main(panels: str | None, out_path: str, previous: str | None):
    """
    if present, reads in any prior reference data
    if present, reads additional panels to use
    queries panelapp for each panel in turn, aggregating results

    Args:
        panels ():
        out_path ():
        previous ():
    """

    logging.info('Starting PanelApp Query Stage')

    # create old_data - needs to be a json in this format
    # absolutely no need for this to be global?
    # {'genes': {'ENSG***': {'panels': [1, 2, 3]}}}
    old_data = {'genes': {}}
    new_genes = None
    if previous:
        logging.info(f'Reading legacy data from {previous}')
        old_data = read_json_from_path(previous)

    elif get_config()['dataset_specific'].get('historic_results'):
        old_file = find_latest_file(start='panel_')
        if old_file is not None:
            logging.info(f'Grabbing legacy panel data from {old_file}')
            old_data = read_json_from_path(old_file)

    else:
        logging.info('No prior data found, running panel diff vs. 12 months ago...')
        older_version = find_core_panel_version()
        if older_version is None:
            raise ValueError('Could not find a version from 12 months ago')
        new_genes = get_new_genes(older_version=older_version)

        logging.info(f'New genes in prev. 12 months: {len(new_genes)}')
        print(new_genes)

    # set up the gene dict
    gene_dict: PanelData = {'metadata': [], 'genes': {}}

    # first add the base content
    get_panel_green(gene_dict, old_data=old_data)

    # if participant panels were provided, add each of those to the gene data
    if panels is not None:
        panel_list = read_panels_from_participant_file(panels)
        logging.info(f'All additional panels: {", ".join(map(str, panel_list))}')
        for panel in panel_list:

            # skip mendeliome
            if panel == get_config()['workflow'].get('default_panel', 137):
                continue

            get_panel_green(gene_dict=gene_dict, old_data=old_data, panel_id=panel)

    # now get the best MOI
    get_best_moi(gene_dict['genes'])

    # if we didn't have prior reference data, scrub down new statuses
    if new_genes is not None:
        overwrite_new_status(gene_dict, new_genes)

    # write the output to long term storage
    write_output_json(output_path=out_path, object_to_write=gene_dict)

    # remove edge case where gene is present, missing for one analysis, then
    # new if it is seen again
    for gene, data in old_data['genes'].items():
        if gene not in gene_dict['genes']:
            gene_dict['genes'][gene] = data

    save_new_historic(gene_dict, prefix='panel_')


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    main()
