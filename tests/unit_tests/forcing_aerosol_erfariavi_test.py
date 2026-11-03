"""Module for aviation ERFari tests."""

import numpy as np

from fair.forcing.aerosol.erfariavi import calculate_erfariavi_forcing

def test_aerosol_erfariavi_linear():

    result = calculate_erfariavi_forcing(
        emissions=(
            np.array([
                3.29469043281061,      
                0.2987951381038286,    
                0.009671277770680551,  
            ])
            * np.ones((1, 1, 1, 3))
        ),
        forcing_scaling=1.0,
        linear_radiative_efficiency=(
            np.array([
                -0.0017640043363624,   
                -0.0121358061533981,   
                 0.1193628593608486,  
            ])
            * np.ones((1, 1, 1, 3))
        ),
        nonlinear_radiative_efficiency=np.zeros((1, 1, 1, 3)),
        linear_emissions_indices=np.array([True, True, True]),
        nonlinear_emissions_indices=np.array([False, False, False]),
    )

    print(result[0, 0, 0, 0])

    np.testing.assert_almost_equal(
        result[0, 0, 0, 0],
        -0.008283576717674062,
    )


def test_aerosol_erfariavi_nonlinear():

    result = calculate_erfariavi_forcing(
        emissions=(
            np.array([
                3.29469043281061,      
                0.2987951381038286,    
                0.009671277770680551,  
            ])
            * np.ones((1, 1, 1, 3))
        ),
        forcing_scaling=1.0,
        linear_radiative_efficiency=np.zeros((1, 1, 1, 3)),
        nonlinear_radiative_efficiency=(
            np.array([
                -0.0039865548220842,
                -0.0039865548220842,
                 0.0,
            ])
            * np.ones((1, 1, 1, 3))
        ),
        linear_emissions_indices=np.array([False, False, False]),
        nonlinear_emissions_indices=np.array([True, True, False]),
    )

    print(result[0, 0, 0, 0])

    np.testing.assert_almost_equal(
        result[0, 0, 0, 0],
        -0.007849027988839437,
    )