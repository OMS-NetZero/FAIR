Introduction
============

The main interaction with FaIR is through the ``FAIR`` class instance::

    from fair import FAIR
    f = FAIR()

``FAIR()`` contains the attributes and methods required to run the model.

Internally, the book-keeping is done with ``numpy`` arrays that have up to five dimensions::

    concentration :: [time, scenario, config, specie, box]

For input and output, we use ``xarray``, which wraps the numpy array with axis labels
and makes model input and output a little more tractable.
