# Genome-wide Variant Prioritisation

Code Repository for the Rare Disease variant prioritisation pipeline

review [design docs](design_docs) for component descriptions

## Configuration

Configuration settings for this application are provided using the [cpg-utils config functions](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/config.py).
There are two provided configuration templates:

1. The Cohort/Run-specific config [here](reanalysis/reanalysis_cohort.toml)
   - contains settings like which seqr instance to use, and the corresponding project ID

2. The Operational settings config [here](reanalysis/reanalysis_global.toml)
   - contains filter thresholds and parameters

The use of the `cpg-utils` functionality removes the need to manually pass a configuration file around the codebase,
instead generating a single aggregated config file and referencing globally as required.

The main considered use-case is via the CPG Analysis-Runner framework, which handles the merging and storage of config
files, but if required outside of this framework, an issue can be raised and will be addressed.

## Flow

The current implementation consists of the following steps:

### [Annotation](reanalysis/annotation.py)

Summarised here [Annotation stage](design_docs/Annotation.md)

This stage can take a VCF, or an MT  as input (pre-processing step will revert to VCF prior to running annotation)

#### Inputs

- String: path to a VCF
- (Implicit): requires annotation resources set up in expected format (waiting for more documentation)

#### Process

1. We calculate a number of intervals across the genome, so that VCF annotation can be split and parallelised
2. For each interval, we use VEP to annotate the variant data, writing output in JSON format
3. The interval annotation data is parsed using Hail, based on a set VEP schema, creating a Hail Table of annotations
4. Once all annotation jobs are complete, the temporary tables are joined into one large table, containing all the consequence annotations
5. The original VCF is loaded into a Hail MatrixTable, and in a single loop we annotate each variant with population frequencies, variant consequences, and clinvar records

#### Outputs

The resultant MT containing all annotations, written as a `.mt` directory

### [PanelApp Query](reanalysis/query_panelapp.py)

Summarised here [PanelApp stage](design_docs/PanelApp_interaction.md)

#### Inputs

- String: a panelapp panel ID
- Optional(String): path to gene list file
- Optional(String): prior panel version

#### Process

1. Query PanelApp for the latest version of the defined panel
2. Store all Green gene data (ENSG, Symbol, MOI)
3. Optional (if prior data is provided) annotate each current gene entity with 'New' if it has been added since that
   original data

#### Outputs

- JSON: PanelApp data for use downstream

### [Labelling in Hail](reanalysis/hail_filter_and_label.py)

Summarised here [Labelling stage](design_docs/Hail_Filter_and_Label.md)

#### Inputs

- MatrixTable: Annotated Variant Data
- String: Path to PanelApp data (JSON)


#### Process

1. Removal of all variants not on PanelApp 'Green' genes
2. Filtering of the Annotated data to remove common or low quality variants
3. Labelling of individual variants with categories of interest
4. Removal of all un-labelled variants
5. Construction of CSQ strings as annotations

#### Outputs

- VCF: labelled variants

### [MOI Checking](reanalysis/validate_categories.py)

##### Inputs

- VCF: labelled variants
- String: path to PanelApp data (JSON)
- String: path to PED file for cohort

#### Process

1. iterate through the VCF, pulling out the labels & annotations for each variant
2. for each variant in turn:
    1. look up the relevant gene in the PanelApp data to find the MOI to use
    2. for each participant which has this variant, see if the MOI fits with their family structure
    3. If this search succeeds record the variant,  participant, & passing MOI

#### Outputs

- JSON: a list of all confirmation records
