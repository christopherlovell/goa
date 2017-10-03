# Galaxy Overdensity Analysis (GOA)

## Instructions

To download all data used in the analysis you must have access to the [German Astronomical Virtual Observatory](http://www.g-vo.org/), Millennium [database](http://gavo.mpa-garching.mpg.de/MyMillennium/).

You can then run the bash scripts in the `queries/` folder.

    cd queries
    ./queries.sh
    ./z0_query.sh

To compile the `norm_coods` cythonised module:

    python setup.py build_ext --inplace


## Selection

We use four different galaxy selections in the paper:

- $\mathrm{S_{SFR1}}$
- $\mathrm{S_{SFR5}}$
- $\mathrm{S_{MAS9}}$
- $\mathrm{S_{MAS10}}$

