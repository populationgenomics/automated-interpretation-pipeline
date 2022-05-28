"""
test class for AIP comparisons
"""


import logging
from typing import List

import pytest
from peddy import Ped

import hail as hl

from reanalysis.comparison import (
    check_gene_is_green,
    apply_variant_qc_methods,
    check_in_vcf,
    common_format_from_results,
    common_format_from_seqr,
    filter_sample_by_ab,
    # filter_to_well_normalised,
    find_missing,
    find_probands,
    find_variant_in_mt,
    # test_consequences,
    CommonFormatResult,
    Confidence,
)


def test_parse_aip(output_json):
    """
    tests that the AIP output JSON is parsed correctly
    :param output_json:
    :return:
    """

    parsed_result = common_format_from_results(output_json)
    assert list(parsed_result.keys()) == ['SAMPLE_1']

    parsed_variants = parsed_result['SAMPLE_1']
    assert isinstance(parsed_variants, list)
    assert len(parsed_variants) == 1

    parsed_var = parsed_variants[0]
    assert isinstance(parsed_var, CommonFormatResult)

    test_variant = CommonFormatResult('7', 93105286, 'T', 'A', [Confidence.EXPECTED])
    assert test_variant == parsed_var


def test_proband_finder(trio_ped):
    """
    tests function to find probands from a Ped file
    this trio contains PROBAND, MOTHER, and FATHER.
    Only the proband is affected, and the mother and father both have a listed child
    not the toughest test
    :param trio_ped:
    :return:
    """
    # digest that Ped
    ped_parsed = Ped(trio_ped)
    probands = find_probands(ped_parsed)
    assert probands == ['PROBAND']


def test_proband_finder_with_sibling(quad_ped):
    """
    tests function to find probands from a Ped file
    contains the same trio as above ^^ plus an unaffected sibling
    :param quad_ped:
    :return:
    """
    # digest that Ped
    ped_parsed = Ped(quad_ped)
    probands = find_probands(ped_parsed)
    assert probands == ['PROBAND']


def test_seqr_parser(seqr_csv_output):
    """
    First variant is 17:10697288 G>A & tagged as Possible
    Second variant does not have an AIP training tag, so should be ignored
    :param seqr_csv_output:
    :return:
    """

    seqr_results = common_format_from_seqr(seqr_csv_output, probands=['PROBAND'])

    # only results for one sample
    assert len(list(seqr_results.keys())) == 1

    sample_variants = seqr_results.get('PROBAND')
    assert (
        len(sample_variants) == 1
    ), f'Should only be one retained variant: {sample_variants}'

    assert sample_variants == [
        CommonFormatResult('17', 10697288, 'G', 'A', [Confidence.POSSIBLE])
    ]


def test_find_missing_matched(caplog):
    """
    trial the comparison process, check logged results
    one matching variant, and one bonus AIP result
    :param caplog:
    :return:
    """
    caplog.set_level(logging.INFO)

    fake_seqr = {'match': [CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE])]}

    fake_aip = {
        'match': [
            CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE]),
            CommonFormatResult('2', 2, 'A', 'G', [Confidence.POSSIBLE]),
        ]
    }

    discrep = find_missing(fake_aip, fake_seqr)
    assert len(discrep) == 0
    log_records = [rec.message for rec in caplog.records]
    assert 'Sample match - 1 matched variant(s)' in log_records
    assert 'Sample match - 0 missing variant(s)' in log_records


def test_find_missing_fails(caplog):
    """
    trial the comparison process, check logged results
    one mis-matching variant
    :param caplog:
    :return:
    """
    caplog.set_level(logging.INFO)

    fake_seqr = {'match': [CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE])]}

    fake_aip = {
        'match': [
            CommonFormatResult('2', 2, 'A', 'G', [Confidence.POSSIBLE]),
        ]
    }

    discrep = find_missing(fake_aip, fake_seqr)
    assert len(discrep) == 1
    assert len(discrep['match']) == 1
    log_records = [rec.message for rec in caplog.records]
    assert 'Sample match - 0 matched variant(s)' in log_records
    assert 'Sample match - 1 missing variant(s)' in log_records


def test_find_missing_different_sample(caplog):
    """
    trial the comparison process, check logged results
    no common samples
    :param caplog:
    :return:
    """
    caplog.set_level(logging.INFO)

    fake_seqr = {'match': [CommonFormatResult('1', 1, 'A', 'C', [Confidence.POSSIBLE])]}

    fake_aip = {
        'mismatch': [
            CommonFormatResult('2', 2, 'A', 'G', [Confidence.POSSIBLE]),
        ]
    }

    discrep = find_missing(fake_aip, fake_seqr)
    assert len(discrep) == 1
    assert len(discrep['match']) == 1
    log_records = [rec.message for rec in caplog.records]
    assert 'Samples completely missing from AIP results: match' in log_records
    assert 'Sample match: 1 missing variant(s)' in log_records


def test_missing_in_vcf(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    :return:
    """
    common_var = CommonFormatResult('1', 1, 'GC', 'G', [Confidence.EXPECTED])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(in_vcf) == 1
    assert in_vcf['SAMPLE'] == [common_var]
    assert len(not_in_vcf) == 0


def test_missing_in_vcf_confidence_irrelevant(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    :return:
    """
    common_var = CommonFormatResult('1', 1, 'GC', 'G', [Confidence.UNLIKELY])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(in_vcf) == 1
    assert in_vcf['SAMPLE'] == [common_var]
    assert len(not_in_vcf) == 0


def test_missing_in_vcf_allele_mismatch(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    alleles should not match
    :return:
    """
    common_var = CommonFormatResult('1', 1, 'GC', 'A', [Confidence.EXPECTED])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(not_in_vcf) == 1
    assert not_in_vcf['SAMPLE'] == [common_var]
    assert len(in_vcf) == 0


def test_missing_in_vcf_variant_pos_mismatch(single_variant_vcf_path):
    """
    test method which scans VCF for a given variant
    alleles should not match
    :return:
    """
    common_var = CommonFormatResult('11', 11, 'GC', 'G', [Confidence.EXPECTED])
    variant_object = {'SAMPLE': [common_var]}
    in_vcf, not_in_vcf = check_in_vcf(single_variant_vcf_path, variant_object)
    assert len(not_in_vcf) == 1
    assert not_in_vcf['SAMPLE'] == [common_var]
    assert len(in_vcf) == 0


def test_variant_in_mt(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    :return:
    """
    query_variant = CommonFormatResult('1', 1, 'GC', 'G', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 1


def test_variant_in_mt_allele_mismatch(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    :return:
    """
    query_variant = CommonFormatResult('1', 1, 'A', 'T', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 0


def test_variant_in_mt_wrong_chrom(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    :return:
    """
    query_variant = CommonFormatResult('11', 1, 'GC', 'G', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 0


def test_variant_in_mt_wrong_pos(hail_matrix):
    """
    check that a variant can be retrieved from a MT
    :return:
    """
    query_variant = CommonFormatResult('1', 100, 'GC', 'G', [])
    result = find_variant_in_mt(hail_matrix, query_variant)
    assert result.count_rows() == 0


# - - - - - - - - - - - - - - - - - - - #
# tests relating to MT category/quality #
# - - - - - - - - - - - - - - - - - - - #


@pytest.mark.parametrize('gene,rows', (('green', 1), ('other', 0)))
def test_check_gene_is_green(gene, rows, hail_matrix):
    """

    :param gene:
    :param rows:
    :param hail_matrix:
    :return:
    """
    green_genes = hl.set(['green'])
    gene_mt = hail_matrix.annotate_rows(geneIds=gene)
    result = check_gene_is_green(gene_mt, green_genes)
    assert result.count_rows() == rows


@pytest.mark.parametrize(
    'filters,alleles,ac,an,results',
    [
        (
            None,
            ['A', 'C'],
            100,
            100,
            ['QC: AC too high in joint call'],
        )
    ],
)
def test_quality_method(
    filters: List[str],
    alleles: List[str],
    ac: int,
    an: int,
    results: List[str],
    hail_matrix: hl.MatrixTable,
):
    """
    required fields: filters, alleles, AC, AN
    then some extra fluff about AB (GT and AD)

    :param filters:
    :param alleles:
    :param ac:
    :param an:
    :param results:
    :param hail_matrix:
    :return:
    """
    config = {'min_samples_to_ac_filter': 0, 'ac_threshold': 0.1}
    if filters is None:
        filters = hl.empty_set(t=hl.tstr)
    anno_mt = hail_matrix.annotate_rows(
        alleles=alleles,
        filters=filters,
        info=hail_matrix.info.annotate(AC=ac, AN=an),
    )
    assert apply_variant_qc_methods(anno_mt, config) == results


@pytest.mark.parametrize(
    'gt,ad,result',
    [
        ('1/0', [20, 80], ['QC: Variant fails AB ratio']),
        ('1/0', [50, 50], []),
        ('0/0', [80, 20], ['QC: Variant fails AB ratio']),
        ('1/1', [80, 20], ['QC: Variant fails AB ratio']),
        ('1|1', [80, 20], ['QC: Variant fails AB ratio']),
    ],
)
def test_filter_sample_by_ab(gt, ad, result, hail_matrix):
    """

    :param gt:
    :param ad:
    :param result:
    :return:
    """

    anno_mt = hail_matrix.annotate_entries(AD=hl.array(ad), GT=hl.parse_call(gt))
    assert filter_sample_by_ab(anno_mt, 'SAMPLE') == result
