#!/usr/bin/env python3


"""
track down the latest version of all reports
generate an index HTML page with links to all reports

Generate a second report for the latest variant only report
"""


import logging
import sys
from dataclasses import dataclass
from os.path import join
from pathlib import Path
from typing import Any

import jinja2

from cpg_utils import to_path
from cpg_utils.config import get_config
from metamist.graphql import gql, query

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'
PROJECT_QUERY = gql(
    """
    query MyQuery {
        myProjects {
            dataset
        }
    }
    """
)
REPORT_QUERY = gql(
    """
    query MyQuery($project: String!) {
        project(name: $project) {
            analyses(active: {eq: true}, type:  {eq: "aip-report"}) {
                output
                meta
                timestampCompleted
            }
        }
    }
    """
)


@dataclass
class Report:
    """
    generic object for storing report details
    """

    dataset: str
    address: str
    genome_or_exome: str
    subtype: str
    date: str


def get_my_projects():
    """
    queries metamist for projects I have access to,
    returns the dataset names
    """
    response: dict[str, Any] = query(PROJECT_QUERY)
    return {dataset['dataset'] for dataset in response['myProjects']}


def get_project_analyses(project: str) -> list[dict]:
    """
    find all the active analysis entries for this project
    Args:
        project (str): project to query for
    """

    response: dict[str, Any] = query(REPORT_QUERY, variables={'project': project})
    return response['project']['analyses']


def main(latest: bool = False):
    """
    finds all existing reports, generates an HTML file
    eventually we can latch onto the meta key `type:aip_output_html`
    but that won't be populated until we run more through the pipeline

    Args:
        latest (bool): whether to create the latest-only report
    """

    all_cohorts = {}

    for cohort in get_my_projects():
        if 'test' in cohort:
            continue

        for analysis in get_project_analyses(cohort):
            output_path = analysis['output']
            # mutually exclusive conditional search for 'latest'
            if latest:
                if 'latest' not in output_path:
                    continue
                date = output_path.rstrip('.html').split('_')[-1]
                cohort_key = f'{cohort}_{date}'
            else:
                if 'latest' in output_path:
                    continue
                cohort_key = cohort

            # pull the exome/singleton flags
            exome_output = 'Exome' if 'exome' in output_path else 'Familial'
            singleton_output = 'Singleton' if 'singleton' in output_path else 'Familial'
            try:
                # pipeline runs don't have display_url
                report_address = analysis['output'].replace(
                    get_config()['storage']['default']['web'],
                    get_config()['storage']['default']['web_url'],
                )
                all_cohorts[f'{cohort_key}_{exome_output}_{singleton_output}'] = Report(
                    dataset=cohort,
                    address=report_address,
                    genome_or_exome=exome_output,
                    subtype=singleton_output,
                    date=analysis['timestampCompleted'].split('T')[0],
                )
            except KeyError:
                logging.info(
                    'Failed to construct a Report entry - is this a report HTML entry?'
                )
                logging.info(analysis)

    # if there were no reports, don't bother with the HTML
    if not all_cohorts:
        return

    # smoosh into a list for the report context - all reports sortable by date
    template_context = {'reports': list(all_cohorts.values())}

    # build some HTML
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
    )
    template = env.get_template('index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    to_path(
        join(
            get_config()['storage']['common']['test']['web'],
            'reanalysis',
            'latest_aip_index.html' if latest else 'aip_index.html',
        )
    ).write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    # run once for all main reports, then again for the latest-only reports
    main()
    main(latest=True)
