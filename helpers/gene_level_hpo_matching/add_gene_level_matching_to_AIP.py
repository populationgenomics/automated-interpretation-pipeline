"""
Parse a json output from the AIP and add gene level phenotype matching plus export
some results needed for a ROC style analysis.

THIS SCRIPT IS NOT FOR PROD USE. It was only intended as a prototype to test if gene level
matching this way might work for us.


Note: jaccard matching has not been implemented properly here as IC is still being used to
select the high scoring pairs. There is an argument to change this in the semsimian call
but I have not figured out what required argument is 🙃.
"""

import argparse
import json
from collections import defaultdict
from semsimian import Semsimian
import csv


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", type=str, help="AIP results JSON file")
    parser.add_argument("outfile", type=str, help="Modified AIP results JSON file")
    parser.add_argument(
        "genes_to_phenotype_path",
        type=str,
        help="path to genes_to_phenotype.txt. Download from here: https://hpo.jax.org/app/data/annotations#:~:text=download-,GENES,-TO%20PHENOTYPE",
    )
    parser.add_argument(
        "phenio_db_path",
        type=str,
        help="path to phenio.db. Download from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz",
    )
    parser.add_argument(
        "--min_score",
        type=float,
        default=14.3,
        help="Minimum similarity score to consider a match",
    )
    parser.add_argument(
        "--similarity_method",
        type=str,
        choices=["IC", "jaccard"],
        default="IC",
        help="method choice. Note: jaccard is not implemented fully here.",
    )
    parser.add_argument(
        "--seqr_tags",
        type=str,
        help="Optional path to seqr_tags for ROC analysis. This is the tsv exported from seqr tagged variants for a project.",
    )

    return parser.parse_args()


def add_gene_level_hpo_matches(
    aip_json,
    min_score,
    similarity_method,
    phenio_db_path,
    phenotypes_by_gene_symbol,
    seqr_tags,
):
    """
    For each variant in the AIP results, compare HPO terms associated to the gene with
    the HPO terms from the case. If any of the case phenotypes is similar to any of the
    gene phenotypes then consider the variant to have a "gene level" match and record this
    in variant object.

    Also collets some dodgy performance stats.
    """

    if similarity_method == "IC":
        similarity_comparison_field = "ancestor_information_content"
    elif similarity_method == "jaccard":
        similarity_comparison_field = "jaccard_similarity"
    else:
        raise ValueError("Invalid similarity method")

    scores_from_cases_with_positive_variants = []
    hpo_tagged_variant_cnt = 0
    all_causitive_variants = 0
    panel_matched_causitive_variants = 0
    gene_matched_causitive_variants = 0
    any_matched_causitive_variants = 0
    total_panel_matches = 0

    semsimian = Semsimian(
        spo=None,
        predicates=["rdfs:subClassOf"],
        resource_path=phenio_db_path,
    )

    for sgid_id, aip_sample_results in aip_json["results"].items():
        case_phenotypes = {
            hpo.split()[0]
            for hpo in aip_sample_results["metadata"]["phenotypes"]
            if hpo.startswith("HP:")
        }

        sg_has_positive_variant = False
        scores_from_sg = []

        for variant in aip_sample_results["variants"]:
            variant_key = "-".join(
                [
                    variant["var_data"]["coordinates"]["chrom"],
                    str(variant["var_data"]["coordinates"]["pos"]),
                    variant["var_data"]["coordinates"]["ref"],
                    variant["var_data"]["coordinates"]["alt"],
                ]
            )

            # Check if the variant is a known positive in seqr_tags
            # This is only matching on the variant key, not the individual. Rough and ready.
            if variant_key in seqr_tags:
                true_positive = True
                sg_has_positive_variant = True
                all_causitive_variants += 1
            else:
                true_positive = False

            # Pull out all gene symbols associated with this variant. This is not the
            # ideal way to do this, but it saves me from having to ensgid lookup in the prototype.
            if "transcript_consequences" in variant["var_data"]:
                variant_genes = set(
                    cons["symbol"]
                    for cons in variant["var_data"]["transcript_consequences"]
                )
            else:
                variant_genes = set()

            # Get all HPO terms associated with the genes for this variant
            variant_phenotypes = set()
            for symbol in variant_genes:
                variant_phenotypes.update(phenotypes_by_gene_symbol[symbol])

            # Skip any variants that don't have phenotypes
            phenotypes_matches = []
            max_IC_score = 0
            max_jaccard_score = 0
            if case_phenotypes and variant_phenotypes:
                # Compare two sets of phenotypes
                termset_similarity = semsimian.termset_pairwise_similarity(
                    case_phenotypes,
                    variant_phenotypes,
                )

                # Convert subject and object terms to lookups
                subject_termset = {
                    term_dict["id"]: term_dict["label"]
                    for term in termset_similarity["subject_termset"]
                    for term_dict in term.values()
                }
                object_termset = {
                    term_dict["id"]: term_dict["label"]
                    for term in termset_similarity["object_termset"]
                    for term_dict in term.values()
                }

                # Find phenotypes matches that meet the min_score threshold
                # Note: subject_best_matches is currently based on IC only, not jaccard.
                for match in termset_similarity["subject_best_matches"][
                    "similarity"
                ].values():
                    if float(match[similarity_comparison_field]) > min_score:
                        # Record the pairs of above threshold matches. Stuffing this into
                        # a string just so it will display in the current html output
                        # without needing to change the template.
                        phenotypes_matches.append(
                            f'"case_term": {subject_termset[match["subject_id"]]} ({match["subject_id"]}), "gene_term": {object_termset[match["object_id"]]} ({match["object_id"]})'.strip().replace(
                                '"', ""
                            )
                        )

                    # Save the max scores for ROC analysis
                    if float(match["ancestor_information_content"]) > max_IC_score:
                        max_IC_score = float(match["ancestor_information_content"])
                    if float(match["jaccard_similarity"]) > max_jaccard_score:
                        max_jaccard_score = float(match["jaccard_similarity"])

            # Add the phenotype matches to the variant object for output
            variant["panels"]["gene_level"] = phenotypes_matches
            if phenotypes_matches:
                hpo_tagged_variant_cnt += 1

            # Save scores for ROC analysis
            scores_from_sg.append(
                {
                    "true_positive": true_positive,
                    "max_IC_score": max_IC_score,
                    "max_jaccard_score": max_jaccard_score,
                }
            )

            # Do some bodgy counting for the summary stats
            if true_positive:
                if phenotypes_matches or variant["panels"]["matched"]:
                    any_matched_causitive_variants += 1
                if phenotypes_matches:
                    gene_matched_causitive_variants += 1
                if variant["panels"]["matched"]:
                    panel_matched_causitive_variants += 1
            if variant["panels"]["matched"]:
                total_panel_matches += 1

        if sg_has_positive_variant:
            scores_from_cases_with_positive_variants.extend(scores_from_sg)

    print(f"Tagged {hpo_tagged_variant_cnt} variants with gene matched HPO terms")
    print(f"Found {all_causitive_variants} known causative variants in AIP results")
    print(
        f"Matched {any_matched_causitive_variants} ({any_matched_causitive_variants/all_causitive_variants*100:.2f}%) causative variants with either gene or panel level matching"
    )
    print(
        f"Matched {gene_matched_causitive_variants} ({gene_matched_causitive_variants/all_causitive_variants*100:.2f}%) causative variants with gene level matching"
    )
    print(
        f"Matched {panel_matched_causitive_variants} ({panel_matched_causitive_variants/all_causitive_variants*100:.2f}%) causative variants with panel level matching"
    )
    print(
        f"Total variants (known-causative + others) with gene level panel match: {total_panel_matches}"
    )

    return aip_json, scores_from_cases_with_positive_variants


def parse_genes_to_phenotype(genes_to_phenotype_file):
    """Parse the genes to phenotype file and return a dict of gene_symbol -> set of HPO ids"""
    gene_to_phenotype = defaultdict(set)
    with open(genes_to_phenotype_file) as f:
        for line in f:
            ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, disease_id = (
                line.split("\t")
            )
            gene_to_phenotype[gene_symbol].add(hpo_id)
    return gene_to_phenotype


def parse_seqr_tags(seqr_tags_file):
    """
    Parse the seqr tags and return a dict of variant -> set of tags

    `positive_tags` set is used to decide which tags to keep. May need adjustment for different projects.
    """
    seqr_tags = {}
    positive_tags = {
        "AIP training: Possible",
        "Likely Pathogenic",
        "Pathogenic",
        "AIP training: Expected",
    }
    with open(seqr_tags_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variant = "-".join([row["chrom"], row["pos"], row["ref"], row["alt"]])
            tags = row["tags"].split("|")

            # Only save rows with a AIP training, path or known gene for phenotype tag
            tags = set(row["tags"].split("|"))
            if set(tags) & positive_tags:
                seqr_tags[variant] = tags
    print(f"Found {len(seqr_tags)} positive seqr tags")
    return seqr_tags


def main():
    args = parse_arguments()
    print(
        f"Using similarity method: {args.similarity_method} with min score: {args.min_score}"
    )

    # Load gene to phenotype mapping
    phenotypes_by_gene_symbol = parse_genes_to_phenotype(args.genes_to_phenotype_path)

    # Parse seqr tags
    if args.seqr_tags:
        seqr_tags = parse_seqr_tags(args.seqr_tags)
    else:
        seqr_tags = {}

    with open(args.infile, "r") as f:
        input_json = json.load(f)

    output_json, scores_from_cases_with_positive_variants = add_gene_level_hpo_matches(
        input_json,
        args.min_score,
        args.similarity_method,
        args.phenio_db_path,
        phenotypes_by_gene_symbol,
        seqr_tags,
    )

    with open(args.outfile, "w") as f:
        json.dump(output_json, f, indent=4)

    # Save the scores for ROC analysis
    if scores_from_cases_with_positive_variants:
        with open(
            args.infile.rsplit(".json", 1)[0] + ".pos_variant_scores.json", "w"
        ) as f:
            json.dump(scores_from_cases_with_positive_variants, f, indent=4)


if __name__ == "__main__":
    main()
