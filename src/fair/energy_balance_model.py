"""n-layer energy balance representation of Earth's climate."""

import numpy as np
import scipy.linalg
import scipy.sparse.linalg
import scipy.stats
import xarray as xr

from .constants import DOUBLING_TIME_1PCT
from .earth_params import earth_radius, seconds_per_year


class EnergyBalanceModel:
    """Energy balance model that converts forcing to temperature.

    The energy balance model is converted to an impulse-response formulation
    (hence the IR part of FaIR) to allow efficient evaluation. The benefits of
    this are increased as once derived, the "layers" of the energy balance
    model do not communicate with each other. The model description can be
    found in [Leach2021]_, [Cummins2020]_, [Tsutsui2017]_ and [Geoffroy2013]_.

    Parameters
    ----------
    ocean_heat_capacity : ``np.ndarray``
        Ocean heat capacity of each layer (top first), W m-2 yr K-1
    ocean_heat_transfer : ``np.ndarray``
        Heat exchange coefficient between ocean layers (top first). The
        first element of this array is akin to the climate feedback
        parameter, with the convention that stabilising feedbacks are
        positive (opposite to most climate sensitivity literature).
        W m-2 K-1
    deep_ocean_efficacy : float
        efficacy of deepest ocean layer. See e.g. [Geoffroy2013]_.
    forcing_4co2 : float
        effective radiative forcing from a quadrupling of atmospheric
        CO2 concentrations above pre-industrial.
    stochastic_run : bool
        Activate the stochastic variability component from [Cummins2020]_.
    sigma_eta : float
        Standard deviation of stochastic forcing component from [Cummins2020]_.
    sigma_xi : float
        Standard deviation of stochastic disturbance applied to surface
        layer. See [Cummins2020]_.
    gamma_autocorrelation : float
        Stochastic forcing continuous-time autocorrelation parameter.
        See [Cummins2020]_.
    seed : int or None
        Random seed to use for stochastic variability.
    timestep : float
        Time interval of the model (yr)

    Raises
    ------
    ValueError
        if the shapes of ``ocean_heat_capacity`` and ``ocean_heat_transfer`` differ.
    ValueError
        if there are not at least two layers in the energy balance model.
    """

    def __init__(
        self,
        ocean_heat_capacity,
        ocean_heat_transfer,
        deep_ocean_efficacy=1,
        forcing_4co2=8,
        stochastic_run=False,
        sigma_eta=0.5,
        sigma_xi=0.5,
        gamma_autocorrelation=2,
        seed=None,
        timestep=1,
        n_timesteps=1,
    ):
        """Initialise the EnergyBalanceModel."""
        # adjust ocean heat capacity to be a rate: units W m-2 K-1
        self.ocean_heat_transfer = np.asarray(ocean_heat_transfer)
        self.deep_ocean_efficacy = deep_ocean_efficacy
        self.forcing_4co2 = forcing_4co2
        self.stochastic_run = stochastic_run
        self.sigma_eta = sigma_eta
        self.sigma_xi = sigma_xi
        self.gamma_autocorrelation = gamma_autocorrelation
        self.ocean_heat_capacity = np.asarray(ocean_heat_capacity) / timestep
        self.n_temperature_boxes = len(self.ocean_heat_capacity)
        if len(self.ocean_heat_transfer) != self.n_temperature_boxes:
            raise ValueError(
                "ocean_heat_capacity and ocean_heat_transfer must be arrays of the "
                "same shape."
            )
        if self.n_temperature_boxes < 2:
            raise ValueError(
                "At least two ocean layers must be specified in the energy balance "
                "model."
            )
        self.temperature = np.zeros((1, self.n_temperature_boxes + 1))
        self.n_timesteps = n_timesteps
        self.n_matrix = self.n_temperature_boxes + 1
        self.seed = seed
        self.timestep = timestep

    def _eb_matrix(self):
        """Define the matrix of differential equations.

        Returns
        -------
        eb_matrix_eigenvalues : `np.ndarray`
            1D array of eigenvalues of the energy balance matrix.
        eb_matrix_eigenvectors : `np.ndarray`
            2D array of eigenvectors (an array of 1D eigenvectors) of the
            energy balance matrix.
        """
        # two box model
        # [x  x]
        # [x  x]

        # three box model
        # [x  x  0]
        # [x  x ex]
        # [0  x  x]

        # four box model
        # [x  x  0  0]
        # [x  x  x  0]
        # [0  x  x ex]
        # [0  0  x  x]

        # put the efficacy of deep ocean in the right place
        # making a vector avoids if statements
        n_box = self.n_temperature_boxes
        eb_matrix = np.zeros((n_box, n_box))
        epsilon_array = np.ones(n_box)
        epsilon_array[n_box - 2] = self.deep_ocean_efficacy

        # First row
        eb_matrix[0, :2] = [
            -(
                self.ocean_heat_transfer[0]
                + epsilon_array[0] * self.ocean_heat_transfer[1]
            )
            / self.ocean_heat_capacity[0],
            epsilon_array[0]
            * self.ocean_heat_transfer[1]
            / self.ocean_heat_capacity[0],
        ]

        # Last row
        eb_matrix[-1, -2:] = [
            self.ocean_heat_transfer[-1] / self.ocean_heat_capacity[-1],
            -self.ocean_heat_transfer[-1] / self.ocean_heat_capacity[-1],
        ]

        # Intermediate rows where n>2
        for row in range(1, n_box - 1):
            eb_matrix[row, row - 1 : row + 2] = [
                self.ocean_heat_transfer[row] / self.ocean_heat_capacity[row],
                -(
                    self.ocean_heat_transfer[row]
                    + epsilon_array[row] * self.ocean_heat_transfer[row + 1]
                )
                / self.ocean_heat_capacity[row],
                epsilon_array[row]
                * self.ocean_heat_transfer[row + 1]
                / self.ocean_heat_capacity[row],
            ]

        # Prepend eb_matrix with stochastic terms if this is a stochastic run:
        # Cummins et al. (2020) eqs. 13 and 14
        eb_matrix = np.insert(eb_matrix, 0, np.zeros(n_box), axis=0)
        prepend_col = np.zeros(n_box + 1)
        prepend_col[0] = -self.gamma_autocorrelation
        prepend_col[1] = 1 / self.ocean_heat_capacity[0]
        eb_matrix = np.insert(eb_matrix, 0, prepend_col, axis=1)
        return eb_matrix

    @property
    def eb_matrix_d(self):
        """Return the discretised matrix exponential."""
        _eb_matrix_d = scipy.linalg.expm(self._eb_matrix())
        return _eb_matrix_d

    def _forcing_vector(self):
        forcing_vector = np.zeros(self.n_temperature_boxes + 1)
        forcing_vector[0] = self.gamma_autocorrelation
        return forcing_vector

    @property
    def forcing_vector_d(self):
        """Return the discretised forcing vector."""
        return scipy.linalg.solve(
            self._eb_matrix(),
            (self.eb_matrix_d - np.identity(self.n_temperature_boxes + 1))
            @ self._forcing_vector(),
        )

    @property
    def stochastic_d(self):
        """Return the stochastic matrix."""
        # define stochastic matrix
        _stochastic_d = np.zeros((self.n_timesteps, self.n_temperature_boxes + 1))

        # stochastic stuff
        if self.stochastic_run:
            eb_matrix = self._eb_matrix()
            q_mat = np.zeros((self.n_matrix, self.n_matrix))
            q_mat[0, 0] = self.sigma_eta**2
            q_mat[1, 1] = (self.sigma_xi / self.ocean_heat_capacity[0]) ** 2
            # use Van Loan (1978) to compute the matrix exponential
            h_mat = np.zeros((self.n_matrix * 2, self.n_matrix * 2))
            h_mat[: self.n_matrix, : self.n_matrix] = -eb_matrix
            h_mat[: self.n_matrix, self.n_matrix :] = q_mat
            h_mat[self.n_matrix :, self.n_matrix :] = eb_matrix.T
            g_mat = scipy.sparse.linalg.expm(h_mat)
            q_mat_d = (
                g_mat[self.n_matrix :, self.n_matrix :].T
                @ g_mat[: self.n_matrix, self.n_matrix :]
            )
            q_mat_d = q_mat_d.astype(np.float64)
            _stochastic_d = scipy.stats.multivariate_normal.rvs(
                size=self.n_timesteps,
                mean=np.zeros(self.n_matrix),
                cov=q_mat_d,
                random_state=self.seed,
            )

        return _stochastic_d

    def impulse_response(self):
        """Convert the energy balance to impulse response representation."""
        eb_matrix = self._eb_matrix()

        # calculate the eigenvectors and eigenvalues on the energy balance
        # (determininstic) part of the matrix, these are the timescales of responses
        eb_matrix_eigenvalues, eb_matrix_eigenvectors = scipy.linalg.eig(
            eb_matrix[1:, 1:]
        )
        self.timescales = -self.timestep / (np.real(eb_matrix_eigenvalues))
        self.response_coefficients = (
            self.timescales
            * (
                eb_matrix_eigenvectors[0, :]
                * scipy.linalg.inv(eb_matrix_eigenvectors)[:, 0]
            )
            / (self.ocean_heat_capacity[0] * self.timestep)
        )

    def emergent_parameters(self, forcing_2co2_4co2_ratio=0.5):
        """Calculate emergent parameters from the energy balance parameters.

        Parameters
        ----------
        forcing_2co2_4co2_ratio : float
            ratio of (effective) radiative forcing converting a quadrupling of
            CO2 to a doubling of CO2.
        """
        # requires impulse response step
        if not hasattr(self, "timescales"):
            self.impulse_response()
        self.ecs = (
            self.forcing_4co2
            * forcing_2co2_4co2_ratio
            * np.sum(self.response_coefficients)
        )
        self.tcr = (
            self.forcing_4co2
            * forcing_2co2_4co2_ratio
            * np.sum(
                self.response_coefficients
                * (
                    1
                    - self.timescales
                    / DOUBLING_TIME_1PCT
                    * (1 - np.exp(-DOUBLING_TIME_1PCT / self.timescales))
                )
            )
        )

    def add_forcing(self, forcing, timestep):
        """Add a forcing time series to EnergyBalanceModel.

        Parameters
        ----------
        forcing : np.ndarray
            time series of [effective] radiative forcing
        timestep : float
            Model timestep, in years
        """
        self.forcing = forcing
        self.timestep = timestep
        self.n_timesteps = len(forcing)

    def run(self):
        """Run the EnergyBalanceModel."""
        # internal variables
        n_box = self.n_matrix - 1
        forcing_vector = self._forcing_vector()

        # Calculate the matrix exponential
        eb_matrix = self._eb_matrix()
        eb_matrix_d = scipy.linalg.expm(eb_matrix)

        # Solve for temperature
        forcing_vector_d = scipy.linalg.solve(
            eb_matrix, (eb_matrix_d - np.identity(self.n_matrix)) @ forcing_vector
        )

        solution = np.zeros((self.n_timesteps, self.n_matrix))
        solution[0, :] = self.temperature[0, :]
        for i in range(1, self.n_timesteps):
            solution[i, :] = (
                eb_matrix_d @ solution[i - 1, :]
                + forcing_vector_d * self.forcing[i - 1]
                + self.stochastic_d[i - 1, :]
            )

        self.temperature = solution[:, 1:]
        self.stochastic_forcing = solution[:, 0]
        self.toa_imbalance = (
            self.forcing
            - self.ocean_heat_transfer[0] * self.temperature[:, 0]
            + (1 - self.deep_ocean_efficacy)
            * self.ocean_heat_transfer[n_box - 1]
            * (self.temperature[:, n_box - 2] - self.temperature[:, n_box - 1])
        )
        self.ocean_heat_content_change = np.cumsum(
            self.toa_imbalance
            * self.timestep
            * earth_radius**2
            * 4
            * np.pi
            * seconds_per_year
        )


def multi_ebm(
    configs,
    ocean_heat_capacity,
    ocean_heat_transfer,
    deep_ocean_efficacy,
    stochastic_run,
    sigma_eta,
    sigma_xi,
    gamma_autocorrelation,
    seed,
    use_seed,
    forcing_4co2,
    timestep,
    timebounds,
):
    """Create several instances of the EnergyBalanceModel.

    This allows efficient parallel implementation in FaIR.
    We have to use a for loop in this function as is does not look like the linear
    algebra functions in scipy are naturally parallel.

    Parameters
    ----------
    configs : list
        A named list of climate configurations.
    ocean_heat_capacity : `np.ndarray`
        Ocean heat capacity of each layer (top first), W m-2 yr K-1
    ocean_heat_transfer : `np.ndarray`
        Heat exchange coefficient between ocean layers (top first). The
        first element of this array is akin to the climate feedback
        parameter, with the convention that stabilising feedbacks are
        positive (opposite to most climate sensitivity literature).
        W m-2 K-1
    deep_ocean_efficacy : float
        efficacy of deepest ocean layer. See e.g. [Geoffroy2013]_.
    stochastic_run : bool
        Activate the stochastic variability component from [Cummins2020]_.
    sigma_eta : float
        Standard deviation of stochastic forcing component from [Cummins2020]_.
    sigma_xi : float
        Standard deviation of stochastic disturbance applied to surface
        layer. See [Cummins2020]_.
    gamma_autocorrelation : float
        Stochastic forcing continuous-time autocorrelation parameter.
        See [Cummins2020]_.
    seed : int or None
        Random seed to use for stochastic variability.
    use_seed : bool
        Whether or not to use the random seed.
    forcing_4co2 : float
        effective radiative forcing from a quadrupling of atmospheric
        CO2 concentrations above pre-industrial.
    timestep : float
        Time interval of the model (yr)
    timebounds : np.ndarray
        Vector representing the time snapshots to calculate temperature on.
    """
    n_configs = ocean_heat_capacity.shape[0]
    n_layers = ocean_heat_capacity.shape[1]
    n_timebounds = len(timebounds)
    ebms = xr.Dataset(
        {
            "eb_matrix_d": (
                ["config", "eb_dim0", "eb_dim1"],
                np.ones((n_configs, n_layers + 1, n_layers + 1)) * np.nan,
            ),
            "forcing_vector_d": (
                ["config", "eb_dim0"],
                np.ones((n_configs, n_layers + 1)) * np.nan,
            ),
            "stochastic_d": (
                ["timebounds", "config", "eb_dim0"],
                np.ones((n_timebounds, n_configs, n_layers + 1)) * np.nan,
            ),
            "ecs": (["config"], np.ones(n_configs) * np.nan),
            "tcr": (["config"], np.ones(n_configs) * np.nan),
        },
        coords={
            "timebounds": timebounds,
            "config": configs,
            "eb_dim0": np.arange(-1, n_layers),
            "eb_dim1": np.arange(-1, n_layers),
        },
    )

    for i_conf, config in enumerate(configs):
        ebm = EnergyBalanceModel(
            ocean_heat_capacity=ocean_heat_capacity[i_conf, :],
            ocean_heat_transfer=ocean_heat_transfer[i_conf, :],
            deep_ocean_efficacy=deep_ocean_efficacy[i_conf],
            stochastic_run=stochastic_run[i_conf],
            sigma_eta=sigma_eta[i_conf],
            sigma_xi=sigma_xi[i_conf],
            gamma_autocorrelation=gamma_autocorrelation[i_conf],
            seed=seed.data[i_conf] if use_seed[i_conf] else None,
            forcing_4co2=forcing_4co2[i_conf],
            timestep=timestep,
            n_timesteps=n_timebounds,
        )
        ebms["eb_matrix_d"].loc[dict(config=config)] = ebm.eb_matrix_d
        ebms["forcing_vector_d"].loc[dict(config=config)] = ebm.forcing_vector_d
        ebms["stochastic_d"].loc[dict(config=config)] = ebm.stochastic_d
        ebm.emergent_parameters()
        ebms["ecs"].loc[dict(config=config)] = ebm.ecs
        ebms["tcr"].loc[dict(config=config)] = ebm.tcr

    return ebms


def step_temperature(state_old, eb_matrix_d, forcing_vector_d, stochastic_d, forcing):
    """Advance parallel energy balance models forward one timestep.

    Parameters
    ----------
    state_old : np.ndarray
        stacked arrays of forcing and temperature of layers in previous timestep
    eb_matrix_d : np.ndarray
        stacked discretised energy balance matrices
    forcing_vector_d : np.ndarray
        stacked discretised forcing vectors
    _stochastic_d : np.ndarray
        stacked matrices of stochastic internal variability
    forcing: np.ndarray
        stacked vectors of [effective] radiative forcing

    Returns
    -------
    state_new : np.ndarray
        stacked arrays of forcing and temperature of layers in this timestep
    """
    state_new = (
        (eb_matrix_d[0, ...] @ state_old[0, ..., None])[..., 0]
        + forcing_vector_d[0, ...] * forcing[0, ..., 0, None]
        + stochastic_d[0, ...]
    )

    return state_new


def calculate_toa_imbalance_postrun(
    state,
    forcing,
    ocean_heat_transfer,
    deep_ocean_efficacy,
):
    """Calculate top of atmosphere energy imbalance.

    The calculation is performed after the scenario has been run to avoid
    looping, since no dynamic state changes affect the calculation.

    Parameters
    ----------
    state : np.ndarray
        stacked arrays of forcing and temperature of layers across the run
    forcing : np.ndarray
        stacked arrays of [effective] radiative forcing across the run
    ocean_heat_transfer : np.ndarray
        Heat exchange coefficient between ocean layers (top first). The
        first element of this array is akin to the climate feedback
        parameter, with the convention that stabilising feedbacks are
        positive (opposite to most climate sensitivity literature).
        W m-2 K-1
    deep_ocean_efficacy : np.ndarray
        efficacy of deepest ocean layer.

    Returns
    -------
    toa_imbalance : np.ndarray
        Top of atmsophere energy imbalance.
    """
    toa_imbalance = (
        forcing
        - ocean_heat_transfer[..., 0] * state[..., 1]
        + (1 - deep_ocean_efficacy)
        * ocean_heat_transfer[..., -1]
        * (state[..., -2] - state[..., -1])
    )
    return toa_imbalance
