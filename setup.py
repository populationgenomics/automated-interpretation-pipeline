"""
setup.py for the talos package
"""

from setuptools import find_packages, setup

with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


def read_reqs(filename: str) -> list[str]:
    """
    Read requirements from a file, return as a list
    TODO eventually split out the requirements vs. the setup content

    Args:
        filename (str): the requirements file to parse

    Returns:
        list[str]: the requirements
    """
    with open(filename, encoding='utf-8') as filehandler:
        return [line.strip() for line in filehandler if line.strip() and not line.startswith('#')]


setup(
    name='talos',
    description='Centre for Population Genomics Variant Prioritisation',
    long_description=readme,
    version='5.1.3',
    author='Matthew Welland, CPG',
    author_email='matthew.welland@populationgenomics.org.au, cas.simons@populationgenomics.org.au',
    package_data={'talos': ['templates/*.jinja', 'example_config.toml']},
    url='https://github.com/populationgenomics/automated-interpretation-pipeline',
    license='MIT',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=read_reqs('requirements.txt'),
    extras_require={
        'test': read_reqs('requirements-dev.txt'),
        'cpg': read_reqs('requirements-cpg.txt'),
    },
    entry_points={
        'console_scripts': [
            # for use in translating a VEP annotated VCF to a MatrixTable
            'vcf_to_mt = talos.vep_vcf_to_mt:main',
            # CPG internal, scans database for published reports, collects into an index page
            'report_hunter = helpers.report_hunter:run_both',
            # CPG implementation, builds an extended Pedigree format
            'GeneratePED = talos.GeneratePED:main',
            # use the HPO terms to select panels for this analysis
            'GeneratePanelData = talos.GeneratePanelData:cli_main',
            # query PanelApp for those selected panels
            'QueryPanelapp = talos.QueryPanelapp:cli_main',
            # Filter and label a small-variant MatrixTable
            'RunHailFiltering = talos.RunHailFiltering:cli_main',
            # Filter and label a SV MatrixTable
            'RunHailFilteringSV = talos.RunHailFilteringSV:cli_main',
            # Run each of the category-labelled variants through MOI filters
            'ValidateMOI = talos.ValidateMOI:cli_main',
            # CPG internal (?), publish those results as an HTML report
            'CreateTalosHTML = talos.CreateTalosHTML:cli_main',
            # CPG internal (?), generate a file for ingestion by Seqr
            'GenerateSeqrFile = talos.minimise_output_for_seqr:cli_main',
        ],
    },
)
