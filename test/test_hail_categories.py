"""
unit testing collection for the hail MT methods
"""


import pytest

import hail as hl
import pandas as pd

from reanalysis.hail_filter_and_label import (
    annotate_aip_clinvar,
    annotate_category_1,
    annotate_category_2,
    annotate_category_3,
    annotate_category_5,
    annotate_category_support,
    green_and_new_from_panelapp,
    filter_to_population_rare,
    split_rows_by_gene_and_filter_to_green,
    filter_to_categorised,
)


category_1_keys = ['locus', 'clinvar_aip_strong']
category_2_keys = [
    'locus',
    'clinvar_aip',
    'cadd',
    'revel',
    'geneIds',
    'consequence_terms',
]
category_3_keys = ['locus', 'clinvar_aip', 'lof', 'consequence_terms']
support_category_keys = [
    'locus',
    'cadd',
    'revel',
    'mutationtaster',
    'gerp',
    'eigen',
]
hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 0], 0),
        ([hl_locus, 1], 1),
        ([hl_locus, 2], 0),
    ],
)
def test_class_1_assignment(values, classified, hail_matrix):
    """
    use some fake annotations, apply to the single fake variant
    check that the classification process works as expected based
    on the provided annotations
    """
    # cast the input as a dictionary
    row_dict = dict(zip(category_1_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_aip_strong=hail_table[hail_matrix.locus].clinvar_aip_strong,
        )
    )

    anno_matrix = annotate_category_1(anno_matrix)
    assert anno_matrix.info.categoryboolean1.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 1, 0.0, 0.0, 'GREEN', 'missense'], 1),
        ([hl_locus, 0, 0.0, 0.0, 'GREEN', 'missense'], 0),
        ([hl_locus, 1, 99.0, 1.0, 'RED', 'frameshift_variant'], 0),
        ([hl_locus, 1, 99.0, 1.0, 'GREEN', 'frameshift_variant'], 1),
        ([hl_locus, 0, 28.11, 0.0, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 0, 0, 0.8, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 0, 0, 0.0, 'GREEN', 'synonymous'], 0),
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
            clinvar_aip=hail_table[hail_matrix.locus].clinvar_aip,
            cadd=hail_table[hail_matrix.locus].cadd,
            revel=hail_table[hail_matrix.locus].revel,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [hl.Struct(consequence_terms=hl.set([csq]))]
            ),
        ),
    )

    anno_matrix = annotate_category_2(anno_matrix, new_genes=hl.set(['GREEN']))
    assert anno_matrix.info.categoryboolean2.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 0, 'hc', 'frameshift_variant'], 0),
        ([hl_locus, 0, 'HC', 'frameshift_variant'], 1),
        ([hl_locus, 1, 'lc', 'frameshift_variant'], 1),
        ([hl_locus, 1, hl.missing(hl.tstr), 'frameshift_variant'], 1),
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
            clinvar_aip=hail_table[hail_matrix.locus].clinvar_aip,
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

    anno_matrix = annotate_category_3(anno_matrix)
    assert anno_matrix.info.categoryboolean3.collect() == [classified]


@pytest.mark.parametrize(
    'spliceai_score,flag',
    [(0.1, 0), (0.11, 0), (0.3, 0), (0.49, 0), (0.5, 1), (0.69, 1), (0.9, 1)],
)
def test_category_5_assignment(spliceai_score: float, flag: int, hail_matrix):
    """
    :param hail_matrix:
    :return:
    """

    matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(splice_ai_delta=spliceai_score)
    )
    matrix = annotate_category_5(matrix)
    assert matrix.info.categoryboolean5.collect() == [flag]


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
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 0.0, 1], 1),
    ],
)
def test_support_assignment(values, classified, hail_matrix):
    """
    :param values: value order in the class4_keys list
    :param classified: expected classification
    :param hail_matrix: test fixture
    """

    polyphen = values.pop()
    sift = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(support_category_keys, values))

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

    anno_matrix = annotate_category_support(anno_matrix)
    assert anno_matrix.info.categorysupport.collect() == [classified]


def test_green_and_new_from_panelapp():
    """
    check that the set expressions from panelapp data are correct
    this is collection of ENSG names from panelapp
    2 set expressions, one for all genes, one for new genes only
    """

    mendeliome = {
        'ENSG00ABCD': {'new': [1]},
        'ENSG00EFGH': {'new': []},
        'ENSG00IJKL': {'new': [2]},
    }
    green_expression, new_expression = green_and_new_from_panelapp(mendeliome)

    # check types
    assert isinstance(green_expression, hl.SetExpression)
    assert isinstance(new_expression, hl.SetExpression)

    # check content by collecting
    assert sorted(green_expression.collect()[0]) == [
        'ENSG00ABCD',
        'ENSG00EFGH',
        'ENSG00IJKL',
    ]
    assert new_expression.collect()[0] == {'ENSG00ABCD', 'ENSG00IJKL'}


@pytest.mark.parametrize(
    'exomes,genomes,clinvar,length',
    [
        (0, 0, 0, 1),
        (1.0, 0, 0, 0),
        (1.0, 0, 1, 1),
        (0.0001, 0.0001, 0, 1),
        (0.0001, 0.0001, 1, 1),
    ],
)
def test_filter_rows_for_rare(exomes, genomes, clinvar, length, hail_matrix):
    """
    :param hail_matrix:
    :return:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            gnomad_ex_af=exomes, gnomad_af=genomes, clinvar_aip=clinvar
        )
    )
    matrix = filter_to_population_rare(anno_matrix)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'gene_ids,length',
    [
        ({'not_green'}, 0),
        ({'green'}, 1),
        ({'gene'}, 1),
        ({'gene', 'not_green'}, 1),
        ({'green', 'gene'}, 2),
        ({hl.missing(t=hl.tstr)}, 0),
    ],
)
def test_filter_to_green_genes_and_split(gene_ids, length, hail_matrix):
    """
    :param hail_matrix:x
    :return:
    """
    green_genes = hl.literal({'green', 'gene'})
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=hl.literal(gene_ids),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [hl.Struct(gene_id='gene', biotype='protein_coding', mane_select='')]
            ),
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == length


def test_filter_to_green_genes_and_split__consequence(hail_matrix):
    """
    :param hail_matrix:
    :return:
    """

    green_genes = hl.literal({'green'})
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=green_genes,
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        gene_id='green', biotype='protein_coding', mane_select=''
                    ),
                    hl.Struct(gene_id='green', biotype='batman', mane_select='NM_Bane'),
                    hl.Struct(gene_id='green', biotype='non_coding', mane_select=''),
                    hl.Struct(
                        gene_id='NOT_GREEN', biotype='protein_coding', mane_select=''
                    ),
                ]
            ),
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == 1
    matrix = matrix.filter_rows(hl.len(matrix.vep.transcript_consequences) == 2)
    assert matrix.count_rows() == 1


@pytest.mark.parametrize(
    'one,two,three,four,five,support,length',
    [
        (0, 0, 0, 'missing', 0, 0, 0),
        (0, 1, 0, 'missing', 0, 0, 1),
        (0, 0, 1, 'missing', 0, 0, 1),
        (0, 0, 0, 'missing', 0, 1, 1),
        (0, 0, 0, 'not_blank', 0, 0, 1),
        (0, 0, 0, 'missing', 1, 0, 1),
        (0, 1, 1, 'missing', 0, 1, 1),
        (1, 0, 0, 'missing', 0, 1, 1),
    ],
)
def test_filter_to_classified(
    one, two, three, four, five, support, length, hail_matrix
):
    """
    :param hail_matrix:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            categoryboolean1=one,
            categoryboolean2=two,
            categoryboolean3=three,
            categorysample4=four,
            categoryboolean5=five,
            categorysupport=support,
        )
    )
    matrix = filter_to_categorised(anno_matrix)
    assert matrix.count_rows() == length


def test_aip_clinvar_default(clinvar_prepared_mt):
    """
    no private annotations applied
    Args:
        clinvar_prepared_mt ():
    """

    mt = annotate_aip_clinvar(
        hl.read_matrix_table(clinvar_prepared_mt), clinvar='absent'
    )
    assert mt.count_rows() == 2
    assert not [x for x in mt.info.clinvar_aip.collect() if x == 1]
    assert not [x for x in mt.info.clinvar_aip_strong.collect() if x == 1]


@pytest.mark.parametrize(
    'rating,stars,rows,regular,strong',
    [
        ('benign', 0, 1, 0, 0),  # with private data, any benign is removed
        ('benign', 1, 1, 0, 0),
        ('other', 7, 2, 0, 0),
        ('pathogenic', 0, 2, 1, 0),
        ('pathogenic', 1, 2, 1, 1),
    ],
)
def test_annotate_aip_clinvar(
    rating, stars, rows, regular, strong, tmp_path, clinvar_prepared_mt
):
    """
    Test intention
    - take a VCF of two variants w/default clinvar annotations
    - create a single variant annotation table with each run
    - apply the parametrized annotations to the table
    """

    # make into a data frame
    table = hl.Table.from_pandas(
        pd.DataFrame(
            [
                {
                    'locus': hl.Locus(contig='chr20', position=63406931),
                    'alleles': ['C', 'CGG'],
                    'rating': rating,
                    'stars': stars,
                    'allele_id': 1,
                }
            ]
        ),
        key=['locus', 'alleles'],
    )
    table_path = str(tmp_path / 'anno.ht')
    table.write(table_path)

    returned_table = annotate_aip_clinvar(
        hl.read_matrix_table(clinvar_prepared_mt), clinvar=table_path
    )
    assert returned_table.count_rows() == rows
    assert (
        len([x for x in returned_table.info.clinvar_aip.collect() if x == 1]) == regular
    )
    assert (
        len([x for x in returned_table.info.clinvar_aip_strong.collect() if x == 1])
        == strong
    )
