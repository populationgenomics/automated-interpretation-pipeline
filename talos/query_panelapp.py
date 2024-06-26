#!/usr/bin/env python3


"""
Complete revision
"""

# mypy: ignore-errors

from argparse import ArgumentParser

import zoneinfo

from cpg_utils.config import config_retrieve

from talos.models import HistoricPanels, PanelApp, PanelDetail, PanelShort, PhenotypeMatchedPanels
from talos.utils import (
    ORDERED_MOIS,
    find_latest_file,
    get_cohort_config,
    get_json_response,
    get_logger,
    get_simple_moi,
    read_json_from_path,
    save_new_historic,
)

PANELAPP_HARD_CODED_DEFAULT = 'https://panelapp.agha.umccr.org/api/v1/panels'
PANELAPP_BASE = config_retrieve(['panels', 'panelapp'], PANELAPP_HARD_CODED_DEFAULT)
DEFAULT_PANEL = config_retrieve(['panels', 'default_panel'], 137)
TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')


def request_panel_data(url: str) -> tuple[str, str, list]:
    """
    takes care of the panelapp query
    Args:
        url ():

    Returns:
        components of the panelapp response
    """

    panel_json = get_json_response(url)
    panel_name = panel_json.get('name')
    panel_version = panel_json.get('version')
    panel_genes = panel_json.get('genes')

    # log name and version
    get_logger().info(f'{panel_name} version: {panel_version}')

    return panel_name, panel_version, panel_genes


def get_panel_green(
    gene_dict: PanelApp,
    old_data: HistoricPanels | None = None,
    panel_id: int = DEFAULT_PANEL,
    version: str | None = None,
    blacklist: list[str] | None = None,
    forbidden_genes: set[str] | None = None,
):
    """
    Takes a panel number, and pulls all GRCh38 gene details from PanelApp
    For each gene, keep the MOI, symbol, ENSG (where present)

    Args:
        gene_dict (): PanelApp obj to continue populating
        old_data (HistoricPanels): dict of sets - panels per gene
        panel_id (): specific panel or 'base' (e.g. 137)
        version (): version, optional. Latest panel unless stated
        blacklist (): list of symbols/ENSG IDs to remove from this panel
        forbidden_genes (set[str]): genes to remove for this cohort
    """

    if blacklist is None:
        blacklist = []

    if forbidden_genes is None:
        forbidden_genes = set()

    # include the version if required
    panel_url = f'{PANELAPP_BASE}/{panel_id}' + (f'?version={version}' if version else '')

    panel_name, panel_version, panel_genes = request_panel_data(panel_url)

    # add metadata for this panel & version
    gene_dict.metadata.append(PanelShort(**{'name': panel_name, 'version': panel_version, 'id': panel_id}))

    # iterate over the genes in this panel result
    for gene in panel_genes:
        symbol = gene.get('entity_name')

        # only retain green genes
        if gene['confidence_level'] != '3' or gene['entity_type'] != 'gene' or symbol in forbidden_genes:
            continue

        ensg = None
        chrom = None

        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene['gene_data']['ensembl_genes'].items():
            if build.lower() == 'grch38':
                # the ensembl version may alter over time, but will be singular
                ensembl_data = content[list(content.keys())[0]]
                ensg = ensembl_data['ensembl_id']
                chrom = ensembl_data['location'].split(':')[0]

        if ensg is None or ensg in blacklist or symbol in blacklist or ensg in forbidden_genes:
            get_logger().info(f'Gene {symbol}/{ensg} removed from {panel_name}')
            continue

        # check if this is a new gene in this analysis
        new_gene = False
        if old_data:
            # new if this gene & panel wasn't in the old data
            if new_gene := (panel_id not in old_data.genes.get(ensg, {})):
                # add this panel to the gene, so it won't be new next time
                old_data.genes.setdefault(ensg, set()).add(panel_id)

        exact_moi = gene.get('mode_of_inheritance', 'unknown').lower()

        # either update or add a new entry
        if ensg in gene_dict.genes:
            this_gene = gene_dict.genes[ensg]

            # now we find it on this panel
            this_gene.panels.add(panel_id)

            # add this moi to the set
            this_gene.all_moi.add(exact_moi)

            # if this is/was new - it's new
            if new_gene:
                this_gene.new.add(panel_id)

        else:
            # save the entity into the final dictionary
            gene_dict.genes[ensg] = PanelDetail(
                **{
                    'symbol': symbol,
                    'all_moi': {exact_moi},
                    'new': {panel_id} if new_gene else set(),
                    'panels': {panel_id},
                    'chrom': chrom,
                },
            )


def get_best_moi(gene_dict: dict):
    """
    From the collected set of all MOIs, take the most lenient
    If set was empty (i.e. no specific MOI found) find the default

    Default is found autosome/sex chrom aware

    Args:
        gene_dict (): the 'genes' index of the collected dict
    """

    for content in gene_dict.values():
        # accept the simplest MOI
        simplified_mois = get_simple_moi(content.all_moi, chrom=content.chrom)

        # force a combined MOI here
        if 'Biallelic' in simplified_mois and 'Monoallelic' in simplified_mois:
            content.moi = 'Mono_And_Biallelic'

        else:
            # take the more lenient of the gene MOI options
            content.moi = sorted(simplified_mois, key=lambda x: ORDERED_MOIS.index(x))[0]


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--panels', help='JSON of per-participant panels')
    parser.add_argument('--out_path', required=True, help='destination for results')
    parser.add_argument('--dataset', default=None, help='dataset to use, optional')
    args = parser.parse_args()
    main(panels=args.panels, out_path=args.out_path, dataset=args.dataset)


def main(panels: str | None, out_path: str, dataset: str | None = None):
    """
    Queries PanelApp for all the gene panels to use in the current analysis
    queries panelapp for each panel in turn, aggregating results

    Args:
        panels (): file containing per-participant panels
        out_path (): where to write the results out to
        dataset (): optional dataset to use
    """

    get_logger().info('Starting PanelApp Query Stage')

    dataset = dataset or config_retrieve(['workflow', 'dataset'])

    # find and extract this dataset's portion of the config file
    # set the Forbidden genes (defaulting to an empty set)
    forbidden_genes = read_json_from_path(get_cohort_config(dataset).get('forbidden'), set())

    # Cat. 2 is greedy - the lower barrier to entry means we should avoid using it unless
    # there is a prior run to bootstrap from. If there's no history file, there are no 'new' genes in this round
    if old_file := find_latest_file(dataset=dataset, start='panel_'):
        get_logger().info(f'Grabbing legacy panel data from {old_file}')
        old_data = read_json_from_path(old_file, return_model=HistoricPanels)  # type: ignore

    else:
        get_logger().info('No prior data found, not treating anything as new')
        old_data = None

    # are there any genes to skip from the Mendeliome? i.e. only report if in a specifically phenotype-matched panel
    remove_from_core: list[str] = config_retrieve(['panels', 'require_pheno_match'], [])
    get_logger().info(f'Genes to remove from Mendeliome: {",".join(remove_from_core)!r}')

    # set up the gene dict
    gene_dict = PanelApp(genes={})

    # first add the base content
    get_logger().info('Getting Base Panel')
    get_panel_green(gene_dict, old_data=old_data, blacklist=remove_from_core, forbidden_genes=forbidden_genes)

    # if participant panels were provided, add each of those to the gene data
    panel_list: set[int] = set()
    if panels is not None:
        get_logger().info('Reading participant panels')
        hpo_panel_object = read_json_from_path(panels, return_model=PhenotypeMatchedPanels)  # type: ignore
        panel_list = hpo_panel_object.all_panels
        get_logger().info(f'Phenotype matched panels: {", ".join(map(str, panel_list))}')

    # now check if there are cohort-wide override panels
    if extra_panels := get_cohort_config(dataset).get('cohort_panels'):
        get_logger().info(f'Cohort-specific panels: {", ".join(map(str, extra_panels))}')
        panel_list.update(extra_panels)

    for panel in panel_list:
        # skip mendeliome - we already queried for it
        if panel == DEFAULT_PANEL:
            continue

        get_logger().info(f'Getting Panel {panel}')
        get_panel_green(gene_dict=gene_dict, panel_id=panel, old_data=old_data, forbidden_genes=forbidden_genes)

    # now get the best MOI, and update the entities in place
    get_best_moi(gene_dict.genes)

    # write the output to long term storage
    with open(out_path, 'w') as out_file:
        out_file.write(PanelApp.model_validate(gene_dict).model_dump_json(indent=4))

    # Only save here if we have a historic location in config
    save_new_historic(old_data, dataset=dataset, prefix='panel_')


if __name__ == '__main__':
    cli_main()
