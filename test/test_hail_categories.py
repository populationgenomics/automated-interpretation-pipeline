"""
unit testing collection for the hail MT methods
"""

import json
import os
import pytest

import hail as hl
import pandas as pd

from reanalysis.hail_filter_and_label import (
    annotate_category_1,
    annotate_category_2,
    annotate_category_3,
    annotate_category_4,
    green_and_new_from_panelapp,
    filter_rows_for_rare,
    filter_benign_or_non_genic,
    filter_to_green_genes_and_split,
    filter_by_consequence,
    filter_to_categorised,
    annotate_category_4_only,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

# contains a single variant at chr1:1, with minimal info
HAIL_VCF = os.path.join(INPUT, 'single_hail.vcf.bgz')
PANELAPP_FILE = os.path.join(INPUT, 'panel_changes_expected.json')

category_1_keys = ['locus', 'clinvar_sig', 'clinvar_stars']
category_2_keys = [
    'locus',
    'clinvar_sig',
    'cadd',
    'revel',
    'geneIds',
    'consequence_terms',
]
category_3_keys = ['locus', 'clinvar_sig', 'lof', 'consequence_terms']
category_4_keys = [
    'locus',
    'cadd',
    'revel',
    'mutationtaster',
    'gerp',
    'eigen',
]

category_conf = {
    'critical_csq': ['frameshift_variant'],
    'in_silico': {
        'cadd': 28.0,
        'revel': 0.4,
        'polyphen': 0.85,
        'sift': 0.05,
        'gerp': 1.0,
        'eigen': 0.25,
    },
}

hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'pathogenic', 0], 0),
        ([hl_locus, 'pathogenic', 1], 1),
        ([hl_locus, 'conflicting_interpretations_of_pathogenicity', 1], 0),
        ([hl_locus, 'benign', 1], 0),
        ([hl_locus, 'pathogenic&something&else', 2], 1),
    ],
)
def test_class_1_assignment(values, classified, hail_matrix):
    """
    use some fake annotations, apply to the single fake variant
    check that the classification process works as expected based
    on the provided annotations
    :param values:
    :param classified:
    :param hail_matrix:
    """
    # cast the input as a dictionary
    row_dict = dict(zip(category_1_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
            clinvar_stars=hail_table[hail_matrix.locus].clinvar_stars,
        )
    )

    anno_matrix = annotate_category_1(anno_matrix)
    assert anno_matrix.info.Category1.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'pathogenic', 0.0, 0.0, 'GREEN', 'missense'], 1),
        ([hl_locus, 'benign', 0.0, 0.0, 'GREEN', 'missense'], 0),
        ([hl_locus, 'pathogenic', 99.0, 1.0, 'RED', 'frameshift_variant'], 0),
        ([hl_locus, 'pathogenic', 99.0, 1.0, 'GREEN', 'frameshift_variant'], 1),
        ([hl_locus, 'benign', 30.0, 0.0, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 'meh', 28.01, 0.0, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 'meh', 0, 0.41, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 'meh', 0, 0.0, 'GREEN', 'synonymous'], 0),
    ],
)
def test_class_2_assignment(values, classified, hail_matrix):
    """
    the fields in the input are, respectively:

    :param values: locus clinvar_sig cadd revel geneIds consequence_terms
    :param classified:
    :param hail_matrix:
    """

    csq = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(category_2_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=hail_table[hail_matrix.locus].geneIds,
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
            cadd=hail_table[hail_matrix.locus].cadd,
            revel=hail_table[hail_matrix.locus].revel,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [hl.Struct(consequence_terms=hl.set([csq]))]
            ),
        ),
    )

    anno_matrix = annotate_category_2(
        anno_matrix, config=category_conf, new_genes=hl.set(['GREEN'])
    )
    assert anno_matrix.info.Category2.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'benign', 'hc', 'frameshift_variant'], 0),
        ([hl_locus, 'benign', 'HC', 'frameshift_variant'], 1),
        ([hl_locus, 'pathogenic', 'lc', 'frameshift_variant'], 1),
        ([hl_locus, 'pathogenic', hl.missing(hl.tstr), 'frameshift_variant'], 1),
    ],
)
def test_class_3_assignment(values, classified, hail_matrix):
    """
    the fields in the input are, respectively:

    :param values: locus clinvar_sig loftee consequence_terms
    :param classified:
    :param hail_matrix:
    """

    csq = values.pop()
    lof = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(category_3_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        consequence_terms=hl.set([csq]),
                        lof=lof,
                    )
                ]
            ),
        ),
    )

    anno_matrix = annotate_category_3(anno_matrix, config=category_conf)
    assert anno_matrix.info.Category3.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 0.0, 0.0, 'n', 0.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.9, 'n', 0.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 29.5, 0.9, 'n', 0.0, 0.0, 1.0, 0.0], 1),
        ([hl_locus, 0.0, 0.0, 'D', 0.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 0.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 0.0, 0.9], 1),
    ],
)
def test_class_4_assignment(values, classified, hail_matrix):
    """
    :param values: value order in the class4_keys list
    :param classified: expected classification
    :param hail_matrix: test fixture
    """

    polyphen = values.pop()
    sift = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(category_4_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            cadd=hail_table[hail_matrix.locus].cadd,
            eigen_phred=hail_table[hail_matrix.locus].eigen,
            gerp_rs=hail_table[hail_matrix.locus].gerp,
            mutationtaster=hail_table[hail_matrix.locus].mutationtaster,
            revel=hail_table[hail_matrix.locus].revel,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        polyphen_score=polyphen,
                        sift_score=sift,
                    )
                ]
            ),
        ),
    )

    anno_matrix = annotate_category_4(anno_matrix, config=category_conf)
    assert anno_matrix.info.Category4.collect() == [classified]


def test_green_and_new_from_panelapp():
    """
    check that the set expressions from panelapp data are correct
    this is collection of ENSG names from panelapp
    2 set expressions, one for all genes, one for new genes only
    :return:
    """
    with open(PANELAPP_FILE, 'r', encoding='utf-8') as handle:
        panelapp_data = json.load(handle)
        green_expression, new_expression = green_and_new_from_panelapp(panelapp_data)

        # check types
        assert isinstance(green_expression, hl.SetExpression)
        assert isinstance(new_expression, hl.SetExpression)

        # check content by collecting
        assert sorted(list(green_expression.collect()[0])) == [
            'ENSG00ABCD',
            'ENSG00EFGH',
            'ENSG00IJKL',
        ]
        assert list(new_expression.collect()[0]) == ['ENSG00EFGH']


@pytest.mark.parametrize(
    'exac,gnomad,length',
    [
        (0, 0, 1),
        (1.0, 0, 0),
        (0.04, 0.04, 1),
    ],
)
def test_filter_rows_for_rare(exac, gnomad, length, hail_matrix):
    """

    :param hail_matrix:
    :return:
    """
    conf = {'af_semi_rare': 0.05}
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            exac_af=exac,
            gnomad_af=gnomad,
        )
    )
    matrix = filter_rows_for_rare(anno_matrix, conf)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'clinvar,stars,gene_id,length',
    [
        ('not_benign', 0, '', 1),
        ('not_benign', 1, 'gene', 0),
        ('sth_else', 1, 'gene', 1),
        ('', 3, 'gene', 1),
        ('', 3, hl.missing(t=hl.tstr), 0),
    ],
)
def test_filter_benign_genic(clinvar, stars, gene_id, length, hail_matrix):
    """

    :param hail_matrix:
    :return:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=clinvar,
            clinvar_stars=stars,
        ),
        geneIds=gene_id,
    )
    matrix = filter_benign_or_non_genic(anno_matrix)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'gene_ids,length',
    [
        ({'not_green'}, 0),
        ({'green'}, 1),
        ({'gene'}, 1),
        ({'gene', 'not_green'}, 1),
        ({'green', 'gene'}, 2),
    ],
)
def test_filter_to_green_genes_and_split(gene_ids, length, hail_matrix):
    """

    :param hail_matrix:
    :return:
    """
    green_genes = hl.literal({'green', 'gene'})
    anno_matrix = hail_matrix.annotate_rows(geneIds=hl.literal(gene_ids))
    matrix = filter_to_green_genes_and_split(anno_matrix, green_genes)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'gene_ids,gene_id,consequences,biotype,mane_select,length',
    [
        ('green', 'green', 'frameshift_variant', 'protein_coding', '', 1),
        ('green', 'green', 'frameshift_variant', '', 'NM_relevant', 1),
        ('mis', 'match', 'frameshift_variant', 'protein_coding', 'NM_relevant', 0),
        ('green', 'green', 'frameshift_variant', '', '', 0),
    ],
)
def test_filter_by_consequence(
    gene_ids, gene_id, consequences, biotype, mane_select, length, hail_matrix
):
    """

    :param hail_matrix:
    :return:
    """
    conf = {'useless_csq': ['synonymous']}
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=gene_ids,
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        consequence_terms=hl.set([consequences]),
                        biotype=biotype,
                        gene_id=gene_id,
                        mane_select=mane_select,
                    )
                ]
            ),
        ),
    )
    matrix = filter_by_consequence(anno_matrix, conf)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'one,two,three,four,length',
    [
        (0, 0, 0, 0, 0),
        (0, 1, 0, 0, 1),
        (0, 0, 1, 0, 1),
        (0, 0, 0, 1, 1),
        (0, 1, 1, 0, 1),
        (1, 0, 0, 1, 1),
    ],
)
def test_filter_to_classified(one, two, three, four, length, hail_matrix):
    """

    :param hail_matrix:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            Category1=one, Category2=two, Category3=three, Category4=four
        )
    )
    matrix = filter_to_categorised(anno_matrix)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'one,two,three,four,flag',
    [
        (0, 0, 0, 0, 0),
        (0, 1, 0, 0, 0),
        (0, 0, 1, 0, 0),
        (0, 0, 0, 1, 1),
        (0, 1, 1, 0, 0),
        (1, 0, 0, 1, 0),
    ],
)
def test_c4_only_tag(one, two, three, four, flag, hail_matrix):
    """

    :param hail_matrix:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            Category1=one, Category2=two, Category3=three, Category4=four
        )
    )
    matrix = annotate_category_4_only(anno_matrix)
    assert matrix.category_4_only.collect() == [flag]
