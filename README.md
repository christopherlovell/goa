# Galaxy Overdensity Analysis (GOA)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1092728.svg)](https://doi.org/10.5281/zenodo.1092728)
[![arXiv](https://img.shields.io/badge/arXiv-1710.02148-red.svg)](https://arxiv.org/abs/1710.02148)

Compute statistics of the galaxy protocluster population from cosmological simulations.

If you use this software, please cite the following paper:

    @ARTICLE{2018MNRAS.474.4612L,
       author = {{Lovell}, C.~C. and {Thomas}, P.~A. and {Wilkins}, S.~M.},
       title = "{Characterising and identifying galaxy protoclusters}",
       journal = {\mnras},
       archivePrefix = "arXiv",
       eprint = {1710.02148},
       keywords = {galaxies: clusters: general, galaxies: high-redshift, galaxies: statistics},
       year = 2018,
       month = mar,
       volume = 474,
       pages = {4612-4628},
       doi = {10.1093/mnras/stx3090},
       adsurl = {http://adsabs.harvard.edu/abs/2018MNRAS.474.4612L},
       adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

## Instructions

To download all data used in the analysis you must have access to the [German Astronomical Virtual Observatory](http://www.g-vo.org/), Millennium [database](http://gavo.mpa-garching.mpg.de/MyMillennium/).

You can then run the scripts in the `queries/` folder.

    cd queries
    ./queries.sh
    ./z0_query.sh

To compile the `norm_coods` cythonised module:

    python setup.py build_ext --inplace
