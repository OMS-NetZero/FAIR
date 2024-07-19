Basic example
=============

FaIR v2.2 is object-oriented and designed to be more flexible than its
predecessors. This does mean that setting up a problem is different to
before - gone are the days of 60 keyword arguments to ``fair_scm`` and
we now use classes and functions with fewer arguments that in the long
run should be easier to use. Of course, there is a learning curve, and
will take some getting used to. This tutorial aims to walk through a
simple problem using FaIR 2.2.

The structure of FaIR 2.2 centres around the ``FAIR`` class, which
contains all information about the scenario(s), the forcer(s) we want to
investigate, and any configurations specific to each species and the
response of the climate.

Note
----

The code in this introductory block is explanatory and if you try to
copy and paste it it you’ll get errors. The code in this file is
self-contained below the heading “1. Create FaIR instance” below.
Alternatively, check out the repository from GitHub and run this example
notebook in ``jupyter``. Details
`here <https://docs.fairmodel.net/en/latest/install.html>`__.

Some basics
-----------

There are two main parts to running fair. The first is setting up the
problem definition. This tells fair what you’re including in the run,
how long you are running for, how many scenarios and how many (climate)
ensemble members (known as configs).

Setting up the run
~~~~~~~~~~~~~~~~~~

First make a new instance:

::

   from fair import FAIR
   f = FAIR()

Then, we need to add some information about the time horizon of our
model, forcers we want to run with, their configuration (and the
configuration of the climate), and optionally some model control
options:

::

   from fair.io import read_properties
   f.define_time(2000, 2050, 1)
   f.define_scenarios(['abrupt', 'ramp'])
   f.define_configs(['high', 'central', 'low'])
   species, properties = read_properties()
   f.define_species(species, properties)
   f.ghg_method='Myhre1998'

Finally, we tell fair to generate some empty (all NaN) array variables:
emissions, concentrations, forcing, temperature etc.:

::

   f.allocate()

That’s all you need to set up the run.

Fill in data
~~~~~~~~~~~~

The ``f.allocate()`` step creates ``xarray`` DataArrays, some of which
need to be populated. For example, to fill a constant 40 GtCO2/yr
emissions rate for fossil CO2 emissions into the ‘abrupt’ scenario, we
can do

::

   from fair.interface import fill
   fill(f.emissions, 40, scenario='abrupt', specie='CO2 FFI')
   ...

In this section we also want to tell fair things such as the climate
feedback parameter, carbon cycle sensitivities, aerosol forcing
parameters, and literally hundreds of other things you can vary. There
are convenience functions for reading in emissions and parameter sets
from external files.

Run and analyse model
~~~~~~~~~~~~~~~~~~~~~

Finally, the model is run with

::

   f.run()

Results are stored within the ``FAIR`` instance as ``xarray`` DataArrays
or Dataset, and can be obtained such as

::

   print(fair.temperature)

Multiple ``scenarios`` and ``configs`` can be supplied in a ``FAIR``
instance, and due to internal parallelisation is the fastest way to run
the model (100 ensemble members per second for 1750-2100 on my Mac for
an emissions driven run). The total number of scenarios that will be run
is the product of ``scenarios`` and ``configs``. For example we might
want to run three emissions ``scenarios`` – let’s say SSP1-2.6, SSP2-4.5
and SSP3-7.0 – using climate calibrations (``configs``) from the UKESM,
GFDL, MIROC and NorESM climate models. This would give us a total of 12
ensemble members in total which are run in parallel.

Recommended order for setting up a problem
------------------------------------------

In this tutorial the recommended order in which to define a problem is
set out step by step, and is as follows:

1.  Create the ``FAIR`` instance, inititalised with run control options.
2.  Define the time horizon of the problem with ``FAIR.define_time()``
3.  Define the scenarios to be run (e.g. SSPs, IAM emissions scenarios,
    or anything you want) with ``FAIR.define_scenarios()``.
4.  Define the configurations to be run with ``FAIR.define_configs()``.
    A configuration (``config``) is a set of parameters that describe
    climate response and species response parameters. For example you
    might have a ``config`` with high climate sensitivity and strong
    aerosol forcing, and one with low climate sensitivity and weak
    aerosol forcing.
5.  Define which ``specie``\ s will be included in the problem, and
    their properties including the run mode (e.g. emissions-driven,
    concentration driven) with ``FAIR.define_species()``.
6.  Optionally, modify run control options.
7.  Create input and output ``DataArrays`` with ``FAIR.allocate()``.
8.  Fill in the DataArrays (e.g. emissions), climate configs, and
    species configs, by either working directly with the ``xarray`` API,
    or using FaIR-packaged convenience functions like ``fill``,
    ``initialise``, ``f.fill_from_csv`` and ``f.override_defaults``.
9.  Run: ``FAIR.run()``.
10. Analyse results by accessing the DataArrays that are attributes of
    ``FAIR``.

1. Create FaIR instance
~~~~~~~~~~~~~~~~~~~~~~~

We’ll call our instance ``f``: it’s nice and short and the ``fair`` name
is reserved for the module.

.. code:: ipython3

    from fair import FAIR

.. code:: ipython3

    f = FAIR()

2. Define time horizon
~~~~~~~~~~~~~~~~~~~~~~

There are two different time indicators in FaIR: the ``timebound`` and
the ``timepoint``. ``timebound``\ s, as the name suggests, are at the
edges of each time step; they can be thought of as instantaneous
snapshots. ``timepoint``\ s are what happens between time bounds and are
rates or integral quantities.

The main thing to remember is that only ``emissions`` are defined on
``timepoint``\ s and everything else is defined on ``timebound``\ s, and
when we specify the time horizon in our model, we are defining the
``timebound``\ s of the problem.

Secondly, the number of ``timebound``\ s is one more than the number of
``timepoint``\ s, as the start and end points are included in the
``timebound``\ s.

.. code:: ipython3

    # create time horizon with bounds of 2000 and 2050, at 1-year intervals
    f.define_time(2000, 2050, 1)
    print(f.timebounds)
    print(f.timepoints)

3. Define scenarios
~~~~~~~~~~~~~~~~~~~

The scenarios are a list of strings that label the scenario dimension of
the model, helping you keep track of inputs and outputs.

In this example problem we will create two scenarios: an “abrupt”
scenario (where emissions or concentrations change instantly) and a
“ramp” scenario where they change gradually.

.. code:: ipython3

    # Define two scenarios
    f.define_scenarios(["abrupt", "ramp"])
    f.scenarios

4. Define configs
~~~~~~~~~~~~~~~~~

Similarly to the scenarios, the configs are a labelling tool. Each
config has associated climate- and species-related settings, which we
will come to later.

We’ll use three config sets, crudely corresponding to high, medium and
low climate sensitivity.

.. code:: ipython3

    # Define three scenarios
    f.define_configs(["high", "central", "low"])
    f.configs

5. Define species
~~~~~~~~~~~~~~~~~

This defines the forcers – anthropogenic or natural – that are present
in your scenario. A ``species`` could be something directly emitted like
CO2 from fossil fuels, or it could be a category where forcing has to be
calculate from precursor emissions like aerosol-cloud interactions.

Each ``specie`` is assigned a name that is used to distinguish it from
other species. You can call the species what you like within the model
as long as you are consistent. We also pass a dictionary of
``properties`` that defines how each specie behaves in the model.

In this example we’ll start off running a scenario with CO2 from fossil
fuels and industry, CO2 from AFOLU, CH4, N2O, and Sulfur, and Volcanic
forcing (note you don’t need the full 40 species used in v1.1-1.6, and
some additional default ones are included). From these inputs we also
want to determine forcing from aerosol-radiation and aerosol-cloud
interactions, as well as CO2, CH4 and N2O.

To highlight some of the functionality we’ll run CO2 and Sulfur
emissions-driven, and CH4 and N2O concentration-driven. (This is akin to
an ``esm-ssp585`` kind of run from CMIP6, though with fewer species).
We’ll use totally fake data here - this is not intended to represent a
real-world scenario but just to highlight how FaIR works.

Full simulations may have 50 or more species included and the
``properties`` dictionary can get quite large, so it can be beneficial
to edit it in a CSV and load it in. This is what is done here - we have
taken the default ``species_configs_properties`` file and cut it down to
keep only the species we care about.

Note the label you have given each specie must appear in the first
column of this file.

.. code:: ipython3

    from fair.io import read_properties

.. code:: ipython3

    species, properties = read_properties('data/importing-data/species_configs_properties.csv')
    f.define_species(species, properties)

In total, we have 9 species in this model. We want to run

1. CO2 fossil and industry
2. CO2 AFOLU
3. Sulfur

with specified emissions.

We want to run

4. CH4
5. N2O

with specified concentrations. We want to include an external time
series from

6. Volcanic

We also want to calculate forcing from CO2, so we need to declare the
CO2 as a greenhouse gas in addition to its emitted components:

7. CO2

and we want to calculate forcing from aerosol radiation and aerosol
cloud interactions

8. ERFari
9. ERFaci

Let’s examine them:

.. code:: ipython3

    f.species

.. code:: ipython3

    f.properties

``properties`` is just a dictionary and ``species`` just a list. You can
define them from scratch, though as you can see even with nine species
the dictionary gets quite long so it is easier to read it in. In the
``properties`` dictionary, the keys must match the ``species`` that you
have declared.

The ``properties`` dictionary contains five keys:

-  ``type`` defines the species type such as CO2, an aerosol precursor,
   or volcanic forcing; there’s around 20 pre-defined types in FaIR, and
   the ``type`` defines several hard-coded properties about what the
   specie does. Some can only be defined once per scenario, some can
   have multiple species attached to its type (e.g. ``f-gas``). See the
   cell below for a list.
-  ``input_mode``: how the model should be driven with this ``specie``.
   Valid values are ``emissions``, ``concentration``, ``forcing`` or
   ``calculated`` and not all options are valid for all ``type``\ s
   (e.g. running solar forcing with concentrations). ``calculated``
   means that the emissions/concentration/forcing of this specie depends
   on others, for example aerosol radiative forcing needs precursors to
   be emitted.
-  ``greenhouse_gas``: True if the ``specie`` is a greenhouse gas, which
   means that an associated ``concentration`` can be calculated (along
   with some other species-specific behaviours). Note that CO2 emissions
   from fossil fuels or from AFOLU are not treated as greenhouse gases.
-  ``aerosol_chemistry_from_emissions``: Some routines such as aerosols,
   methane lifetime, or ozone forcing, relate to emissions of
   short-lived climate forcers. If this ``specie`` is one of these, this
   should be set to True.
-  ``aerosol_chemistry_from_concentration``: As above, but if the
   production of ozone, aerosol etc. depends on the concentration of a
   greenhouse gas.

.. code:: ipython3

    # Here's a list of the species types: this cell not necessary for running fair
    
    import pandas as pd
    from fair.structure.species import species_types, valid_input_modes, multiple_allowed
    
    types = pd.DataFrame(
        {
            'type': species_types,
            'valid_input_modes': [valid_input_modes[specie] for specie in species_types],
            'multiple_allowed': [multiple_allowed[specie] for specie in species_types]
        }
    )
    types.set_index('type', inplace=True)
    types

6. Modify run options
~~~~~~~~~~~~~~~~~~~~~

When we initialise the FAIR class, a number of options are given as
defaults.

Let’s say we want to change the greenhouse gas forcing treatment from
Meinshausen et al. 2020 default to Myhre et al. 1998. While this could
have been done when initialising the class, we can also do it by setting
the appropriate attribute.

.. code:: ipython3

    help(f)

.. code:: ipython3

    f.ghg_method

.. code:: ipython3

    f.ghg_method='myhre1998'

.. code:: ipython3

    f.ghg_method

7. Create input and output data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Steps 2–5 above dimensioned our problem; now, we want to actually create
some data to put into it.

First we allocate the data arrays with

.. code:: ipython3

    f.allocate()

This has created our arrays with the correct dimensions as attributes of
the ``FAIR`` class:

.. code:: ipython3

    f.emissions

.. code:: ipython3

    f.temperature

8. Fill in the data
~~~~~~~~~~~~~~~~~~~

The data created is nothing more special than ``xarray`` DataArrays.

8a. fill emissions, concentrations …
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using ``xarray`` methods we can allocate values to the emissions. For
example, to fill CO2 fossil emissions in the abrupt scenario with a
constant emissions rate of 38 GtCO2/yr (about present-day levels), do:

::

   f.emissions.loc[(dict(specie="CO2 FFI", scenario="abrupt"))] = 38

To do this for every specie, especially if you want to create
time-varying scenarios and dimension the scenarios correctly, is very
fiddly and time consuming. Therefore, we have a way to read in scenarios
from CSV files:

.. code:: ipython3

    f.fill_from_csv(
        emissions_file='data/basic_run_example/emissions.csv',
        concentration_file='data/basic_run_example/concentration.csv',
        forcing_file='data/basic_run_example/forcing.csv'
    )

Let’s have a look at what we’re reading in:

.. code:: ipython3

    pd.read_csv('data/basic_run_example/emissions.csv')

.. code:: ipython3

    pd.read_csv('data/basic_run_example/concentration.csv')

.. code:: ipython3

    pd.read_csv('data/basic_run_example/forcing.csv')

The csv reader will also interpolate, so you don’t have to specify every
year in your scenario. We have separate functions for reading in files
from the Reduced Complexity Model Intercomparison Project (the SSPs).

Note also we haven’t filled in every species. CO2 is ``calculated`` (the
sum of ``CO2 FFI`` and ``CO2 AFOLU``), which will correctly determine
emissions, concentrations and forcing. Aerosol-radiation interactions
and Aerosol-cloud interactions are also ``calculated``, from emissions
of Sulfur. If the ``properties``, particularly ``type`` and
``input_mode`` are correctly specified for each specie, ``fair`` knows
what to do with your data.

The other thing that we have to do is define the initial conditions of
our data. If you forget to do this, you might get NaN value errors in
``fair``; this is deliberate, we want the user to think about how they
engage with the model!

Using non-zero initial conditions can be useful for “restart runs”:
switching from concentration-driven to emissions-driven (ZECMIP);
running constant forcing commitments; running interative/adaptive
emissions scenarios; the possibilities are endless.

As CO2 is emissions-driven, we also need to define its starting
concentration.

Again, we have a convenience function that can handle some of the heavy
lifting:

.. code:: ipython3

    # Define initial conditions
    from fair.interface import initialise
    
    initialise(f.concentration, 278.3, specie='CO2')
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)

8b. Fill in ``configs``
^^^^^^^^^^^^^^^^^^^^^^^

This defines how the forcing is calculated, and how the model responds
to a forcing.

There are two ``xarray.Dataset``\ s in ``fair`` that define the
behaviour of the model, which are ``f.climate_configs`` and
``f.species_configs``.

Many, many ``species_configs`` parameters have sensible defaults that
would make little impact, or little sense, to change. We can add these
parameters to the ``species_configs_properties`` file that we read in
earlier that would autopopulate most fields. Let’s do this.

.. code:: ipython3

    f.fill_species_configs('data/importing-data/species_configs_properties.csv')

.. code:: ipython3

    f.species_configs

However, there’s no fun if they are all the same for each config - we
want to use fair to sample an ensemble. Again we can read these in, with
columns following a naming convention. This also fills the
``climate_configs``, which are NaN by default (again we don’t want
people running things without thinking).

.. code:: ipython3

    f.override_defaults('data/basic_run_example/configs_ensemble.csv')

and this is how it looks

.. code:: ipython3

    pd.read_csv('data/basic_run_example/configs_ensemble.csv', index_col=0)

9. run FaIR
~~~~~~~~~~~

.. code:: ipython3

    f.run()

10. plot results
~~~~~~~~~~~~~~~~

.. code:: ipython3

    import matplotlib.pyplot as pl

.. code:: ipython3

    pl.plot(f.timebounds, f.temperature.loc[dict(scenario='ramp', layer=0)], label=f.configs)
    pl.title('Ramp scenario: temperature')
    pl.xlabel('year')
    pl.ylabel('Temperature anomaly (K)')
    pl.legend()

.. code:: ipython3

    pl.plot(f.timebounds, f.concentration.loc[dict(scenario='ramp', specie='CO2')], label=f.configs)
    pl.title('Ramp scenario: CO2')
    pl.xlabel('year')
    pl.ylabel('CO2 (ppm)')
    pl.legend()

.. code:: ipython3

    pl.plot(f.timebounds, f.forcing.loc[dict(scenario='ramp', specie='Aerosol-cloud interactions')], label=f.configs)
    pl.title('Ramp scenario: forcing')
    pl.xlabel('year')
    pl.ylabel('ERF from aerosol-cloud interactions (W m$^{-2}$)')
    pl.legend()

.. code:: ipython3

    pl.plot(f.timebounds, f.forcing_sum.loc[dict(scenario='ramp')], label=f.configs)
    pl.title('Ramp scenario: forcing')
    pl.xlabel('year')
    pl.ylabel('Total ERF (W m$^{-2}$)')
    pl.legend()

.. code:: ipython3

    pl.plot(f.timebounds, f.temperature.loc[dict(scenario='abrupt', layer=0)], label=f.configs)
    pl.title('Abrupt scenario: temperature')
    pl.xlabel('year')
    pl.ylabel('Temperature anomaly (K)')
    pl.legend()

.. code:: ipython3

    pl.plot(f.timebounds, f.forcing_sum.loc[dict(scenario='abrupt')], label=f.configs)
    pl.title('Abrupt scenario: forcing')
    pl.xlabel('year')
    pl.ylabel('Total ERF (W m$^{-2}$)')
    pl.legend()

.. code:: ipython3

    pl.plot(f.timebounds, f.concentration.loc[dict(scenario='abrupt', specie='CO2')], label=f.configs)
    pl.title('Abrupt scenario: CO2')
    pl.xlabel('year')
    pl.ylabel('CO2 (ppm)')
    pl.legend()

.. code:: ipython3

    f.species_configs['g0'].loc[dict(specie='CO2')]

.. code:: ipython3

    f.forcing[-1, :, 1, :]

