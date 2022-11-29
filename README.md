## README

This is a fork of the demo code from Patrick Farrell for multilevel Monte Carlo computations in Python.

For mathematical details, see 

    @article{giles2015,
    author = {Giles, M. B.},
    title = {Multilevel {Monte Carlo} methods},
    journal = {Acta Numerica},
    volume = {24},
    month = {5},
    year = {2015},
    pages = {259--328},
    doi = {10.1017/S096249291500001X},
    }

To run out of the source tree (the easiest way):

    $ export PYTHONPATH=/path/to/pymlmc:$PYTHONPATH

or to install to your user-level Python installation:

    $ python setup.py install --user

or to install to a particular directory:

    $ python setup.py install --prefix=/path/to/directory