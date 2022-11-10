"""
tests for clinvar manual summaries
"""


from copy import deepcopy
from datetime import datetime

from clinvar.conflict_huntr import (
    ACMG_THRESHOLD,
    Consequence,
    Submission,
    county_county,
    check_stars,
    acmg_filter_submissions,
)


CURRENT_TIME = datetime.now()
BASIC_SUB = Submission(CURRENT_TIME, 'submitter', Consequence.UNKNOWN, 'review')
BENIGN_SUB = Submission(CURRENT_TIME, 'submitter', Consequence.BENIGN, 'review')
PATH_SUB = Submission(CURRENT_TIME, 'submitter', Consequence.PATHOGENIC, 'review')
UNCERTAIN_SUB = Submission(CURRENT_TIME, 'submitter', Consequence.UNCERTAIN, 'review')


def test_check_stars_practice():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'practice guideline'
    subs = [BASIC_SUB, sub1]
    assert check_stars(subs) == 4


def test_check_stars_expert():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'reviewed by expert panel'
    subs = [BASIC_SUB, sub1]
    assert check_stars(subs) == 3


def test_check_stars_criteria_prov():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'something, something, criteria provided'
    subs = [BASIC_SUB, sub1]
    assert check_stars(subs) == 1


def test_check_stars_none():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'smithsonian'
    subs = [BASIC_SUB, sub1, BASIC_SUB]
    assert check_stars(subs) == 0


def tests_acmg_filter_neutral():
    """filter submissions against ACMG date threshold"""
    subs = [BASIC_SUB, BASIC_SUB]
    assert acmg_filter_submissions(subs) == subs


def tests_acmg_filter_removes():
    """filter submissions against ACMG date threshold"""
    sub1 = deepcopy(BASIC_SUB)
    sub1.date = datetime(year=1970, month=1, day=1)
    sub2 = deepcopy(BASIC_SUB)
    sub2.date = datetime(year=2000, month=1, day=1)
    subs = [BASIC_SUB, sub1, sub2]
    assert acmg_filter_submissions(subs) == [BASIC_SUB]


def tests_acmg_filter_gte():
    """filter submissions against ACMG date threshold"""
    sub1 = deepcopy(BASIC_SUB)
    sub1.date = ACMG_THRESHOLD
    subs = [BASIC_SUB, sub1]
    assert acmg_filter_submissions(subs) == subs


def test_county_empty():
    """No entries == Uncertain"""
    assert county_county(subs=[]) == Consequence.UNCERTAIN


def test_county_all_path():
    """pathogenic if everything is pathogenic"""
    assert county_county([PATH_SUB] * 10) == Consequence.PATHOGENIC


def test_county_all_benign():
    """benign if everything is benign"""
    assert county_county([BENIGN_SUB] * 10) == Consequence.BENIGN


def test_county_all_uncertain():
    """uncertain if everything is uncertain"""
    assert county_county([UNCERTAIN_SUB] * 10) == Consequence.UNCERTAIN


def test_county_equal_path_benign():
    """equal path and benign == conflicting"""
    subs = ([BENIGN_SUB] * 10) + ([PATH_SUB] * 10)
    assert county_county(subs) == Consequence.CONFLICTING


def test_county_path_majority():
    """path on 60-20 split"""
    subs = ([BENIGN_SUB] * 2) + ([PATH_SUB] * 6) + ([UNCERTAIN_SUB] * 2)
    assert county_county(subs) == Consequence.PATHOGENIC


def test_county_path_almost_majority():
    """conflicting on borderline ratios"""
    subs = ([BENIGN_SUB] * 2) + ([PATH_SUB] * 5) + ([UNCERTAIN_SUB] * 2)
    assert county_county(subs) == Consequence.CONFLICTING


def test_county_benign_majority():
    """benign on 60-20 split"""
    subs = ([BENIGN_SUB] * 6) + ([PATH_SUB] * 2) + ([UNCERTAIN_SUB] * 2)
    assert county_county(subs) == Consequence.BENIGN


def test_county_benign_almost_majority():
    """conflicting on borderline ratios"""
    subs = ([BENIGN_SUB] * 5) + ([PATH_SUB] * 2) + ([UNCERTAIN_SUB] * 2)
    assert county_county(subs) == Consequence.CONFLICTING


def test_county_benign_vs_unknown():
    """benign on slim majority vs VUS"""
    subs = ([BENIGN_SUB] * 5) + ([UNCERTAIN_SUB] * 4)
    assert county_county(subs) == Consequence.BENIGN


def test_county_pathogenic_vs_unknown():
    """pathogenic on slim majority vs VUS"""
    subs = ([PATH_SUB] * 5) + ([UNCERTAIN_SUB] * 4)
    assert county_county(subs) == Consequence.PATHOGENIC


def test_county_uncertain_from_mix():
    """VUS if 50% or more VUS"""
    subs = ([PATH_SUB] * 5) + ([UNCERTAIN_SUB] * 5)
    assert county_county(subs) == Consequence.UNCERTAIN


def test_county_uncertain_from_majority():
    """VUS if 50% or more VUS"""
    subs = ([PATH_SUB] * 5) + ([UNCERTAIN_SUB] * 6)
    assert county_county(subs) == Consequence.UNCERTAIN
