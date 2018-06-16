# primeval

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/carjed/primeval/master)

`primeval` is a small Python module with functions for performing coalescent simulations under various models of human demographic history, using the [msprime](https://github.com/tskit-dev/msprime) reimplementation of Richard Hudson’s seminal `ms` program.

## Quick Start

The easiest way to use `primeval` is with [BinderHub](https://mybinder.org/v2/gh/carjed/primeval/master). This will launch a cloud-based computing environment, with interactive Jupyter notebooks to run and explore the various coalescent models.

## Local Install

To use `primeval` locally, simply clone this repository and ensure that the `msprime` library is installed (we recommend installing with `conda`, but `pip` can also be used--see the [`msprime` documentation](https://msprime.readthedocs.io/en/stable/installation.html) for additional instructions):

```
conda config --add channels conda-forge
conda install msprime

git clone https://github.com/carjed/primeval.git
cd primeval
```

Then, from within the Python console or a script, import the `primeval.py` module:

```python
# load module
from primeval import *

# Print debugging info for the Gutenkunst et al. (2009) model
gutenkunst_model(debug = True)
```

Alternatively, the interactive notebooks in this repository can be accessed locally if you have Jupyter installed on your system.

## Current Models

*primeval* currently includes functions for running coalescent simulations under the following models:

- [Gutenkunst et al., 2009](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695)
- [Fu et al., 2012](https://www.nature.com/articles/nature11690)
- [Chen et al., 2015](https://academic.oup.com/mbe/article/32/11/2996/981518)
- Pull requests for other models are welcome!