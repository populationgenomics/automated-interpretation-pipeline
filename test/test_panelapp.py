"""
tests for the PanelApp parser
"""

import json

import pytest

from reanalysis.query_panelapp import (
    get_best_moi,
    get_panel_green,
    read_panels_from_participant_file,
)


@pytest.fixture(name='fake_panelapp')
def fixture_fake_panelapp(requests_mock, latest_mendeliome, latest_incidentalome):
    """
    prepares the web requests mock to serve as stand-in panelapp
    Args:
        requests_mock ():
        latest_mendeliome ():
        latest_incidentalome ():
    """

    requests_mock.register_uri(
        'GET',
        'https://panelapp.agha.umccr.org/api/v1/panels/137',
        json=latest_mendeliome,
    )
    requests_mock.register_uri(
        'GET',
        'https://panelapp.agha.umccr.org/api/v1/panels/126',
        json=latest_incidentalome,
    )


def test_panel_query(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    :param fake_panelapp: fake web hook mock
    """

    gd = {'genes': {}, 'metadata': []}
    get_panel_green(
        gd,
        old_data={
            'genes': {
                'ENSG00ABCD': {'panels': ['pink']},
                'ENSG00EFGH': {'panels': [137]},
            }
        },
    )
    assert gd['genes']['ENSG00ABCD']['moi'] == {'Biallelic'}
    assert gd['genes']['ENSG00ABCD']['panels'] == ['pink', 137]
    assert gd['genes']['ENSG00EFGH']['moi'] == {'Monoallelic'}


def test_panel_query_addition(fake_panelapp):  # pylint: disable=unused-argument
    """
    check that the default parsing delivers correct data
    oof, this was a tricky one
    :param fake_panelapp: fake web hook mock
    """
    # assumed data we already gathered
    gd = {
        'metadata': [{'version': '0.11088', 'name': 'Mendeliome', 'id': 137}],
        'genes': {
            'ENSG00ABCD': {
                'symbol': 'ABCD',
                'moi': {'Monoallelic'},
                'new': [],
                'panels': [137],
            },
            'ENSG00IJKL': {
                'symbol': 'IJKL',
                'moi': {'Mono_And_Biallelic'},
                'new': [137],
                'panels': [123, 137],
            },
        },
    }

    # should query for and integrate the incidentalome content
    get_panel_green(
        gd,
        panel_id=126,
        old_data={
            'genes': {
                'ENSG00EFGH': {'panels': [137, 126]},
                'ENSG00IJKL': {'panels': [137]},
            },
        },
    )
    assert gd['genes']['ENSG00ABCD']['moi'] == {'Monoallelic', 'Biallelic'}
    assert gd['genes']['ENSG00ABCD']['panels'] == [137, 126]
    assert gd['genes']['ENSG00IJKL']['moi'] == {'Mono_And_Biallelic'}
    assert gd['genes']['ENSG00IJKL']['panels'] == [123, 137]
    assert 'ENSG00EFGH' not in gd['genes']


def test_get_list_from_participants(tmp_path):
    """
    tests the unique panel finder
    """
    party_data = {
        'i': {'panels': [1, 2], 'what': 'does'},
        'am': {'panels': [1, 3], 'the': 'fox'},
        'sam': {'panels': [9, 99], 'say?': 'Wa-pa-pa-pa-pa-pa-pow!'},
    }
    tmp_json = tmp_path / 'temp.json'
    with open(tmp_json, 'w', encoding='utf-8') as handle:
        json.dump(party_data, handle)
    assert read_panels_from_participant_file(str(tmp_json)) == {1, 2, 3, 9, 99}


def test_get_best_moi_unknown():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'Unknown'}}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Mono_And_Biallelic'


def test_get_best_moi_mono():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'Monoallelic'}}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Monoallelic'


def test_get_best_moi_mono_and_biallelic():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'Monoallelic', 'Biallelic'}}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Monoallelic'


def test_get_best_moi_1():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'Monoallelic', 'Biallelic', 'both'}}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Mono_And_Biallelic'


def test_get_best_moi_x():
    """
    check that the MOI summary works
    """

    d = {'ensg1': {'moi': {'x-linked biallelic', 'x-linked'}}}
    get_best_moi(d)
    assert d['ensg1']['moi'] == 'Hemi_Mono_In_Female'
