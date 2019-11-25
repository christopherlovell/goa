# Galaxy Overdensity Analysis (GOA)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1092728.svg)](https://doi.org/10.5281/zenodo.1092728)
[![arXiv](https://img.shields.io/badge/arXiv-1710.02148-red.svg)](https://arxiv.org/abs/1710.02148)

Compute statistics of the galaxy protocluster population from cosmological simulations.

If you use this software, please use the following BibTex citation:

    @article{10.1093/mnras/stx3090,
    author = {Lovell, Christopher C and Thomas, Peter A and Wilkins, Stephen M},
    title = "{Characterising and identifying galaxy protoclusters}",
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume = {474},
    number = {4},
    pages = {4612-4628},
    year = {2017},
    month = {12},
    issn = {0035-8711},
    doi = {10.1093/mnras/stx3090},
    url = {https://doi.org/10.1093/mnras/stx3090},
    eprint = {http://oup.prod.sis.lan/mnras/article-pdf/474/4/4612/23011916/stx3090.pdf},
    }

## Instructions

To download all data used in the analysis you must have access to the [German Astronomical Virtual Observatory](http://www.g-vo.org/), Millennium [database](http://gavo.mpa-garching.mpg.de/MyMillennium/).

You can then run the scripts in the `queries/` folder.

    cd queries
    ./queries.sh
    ./z0_query.sh

To compile the `norm_coods` cythonised module:

    python setup.py build_ext --inplace
