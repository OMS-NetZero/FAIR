"""Module for greenhouse gas forcing."""

import numpy as np


def etminan2016(
    concentration,
    baseline_concentration,
    forcing_scaling,
    radiative_efficiency,
    co2_indices,
    ch4_indices,
    n2o_indices,
    minor_greenhouse_gas_indices,
    a1=-2.4e-7,
    b1=7.2e-4,
    c1=-2.1e-4,
    d1=5.36,
    a2=-8e-6,
    b2=4.2e-6,
    c2=-4.9e-6,
    d2=0.117,
    a3=-1.3e-6,
    b3=-8.2e-6,
    d3=0.043,
):
    """Greenhouse gas forcing from concentrations.

    This uses the [Etminan2016]_ relationships for CO2, CH4 and N2O
    including band overlaps between these three gases. Forcing from minor
    greenhouse gases are a linear function of their concentration based on their
    radiative efficiency.

    Note that the relationship can be poorly behaved outside of the range of
    valid values in [Etminan2016]_, (180-2000 ppm CO2, 340-3500 ppb CH4, 200-525
    ppb N2O) particularly for CO2 concentrations above 2000 ppm. It is
    recommended to use "meinshausen2020", which is a re-fit of the coefficients
    and extension to be better behaved outside of the valid concentration range.

    Parameters
    ----------
    concentration : np.ndarray
        concentration of greenhouse gases. "CO2", "CH4" and "N2O" must be
        included in units of [ppm, ppb, ppb]. Other GHGs are units of ppt.
    baseline_concentration : np.ndarray
        pre-industrial concentration of the gases (see above).
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : np.ndarray
        radiative efficiency to use for linear-forcing gases, in W m-2 ppb-1
    co2_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CO2.
    ch4_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CH4.
    n2o_indices : np.ndarray of bool
        index along SPECIES AXIS relating to N2O.
    minor_greenhouse_gas_indices : np.ndarray of bool
        indices of other GHGs that are not CO2, CH4 or N2O.
    a1 : float
        fitting parameter (see [Etminan2016]_)
    b1 : float
        fitting parameter (see [Etminan2016]_)
    c1 : float
        fitting parameter (see [Etminan2016]_)
    d1 : float
        fitting parameter (see [Etminan2016]_)
    a2 : float
        fitting parameter (see [Etminan2016]_)
    b2 : float
        fitting parameter (see [Etminan2016]_)
    c2 : float
        fitting parameter (see [Etminan2016]_)
    d2 : float
        fitting parameter (see [Etminan2016]_)
    a3 : float
        fitting parameter (see [Etminan2016]_)
    b3 : float
        fitting parameter (see [Etminan2016]_)
    d3 : float
        fitting parameter (see [Etminan2016]_)

    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W/m2) from greenhouse gases
    """
    erf_out = np.ones_like(concentration) * np.nan

    # easier to deal with smaller arrays
    co2 = concentration[..., co2_indices]
    ch4 = concentration[..., ch4_indices]
    n2o = concentration[..., n2o_indices]
    co2_base = baseline_concentration[..., co2_indices]
    ch4_base = baseline_concentration[..., ch4_indices]
    n2o_base = baseline_concentration[..., n2o_indices]

    # CO2
    erf_out[..., co2_indices] = (
        (
            a1 * (co2 - co2_base) ** 2
            + b1 * np.abs(co2 - co2_base)
            + c1 * 0.5 * (n2o + n2o_base)
            + d1
        )
        * np.log(co2 / co2_base)
        * (forcing_scaling[..., co2_indices])
    )

    # CH4
    erf_out[..., ch4_indices] = (
        (a3 * 0.5 * (ch4 + ch4_base) + b3 * 0.5 * (n2o + n2o_base) + d3)
        * (np.sqrt(ch4) - np.sqrt(ch4_base))
    ) * (forcing_scaling[..., ch4_indices])

    # N2O
    erf_out[..., n2o_indices] = (
        (
            a2 * 0.5 * (co2 + co2_base)
            + b2 * 0.5 * (n2o + n2o_base)
            + c2 * 0.5 * (ch4 + ch4_base)
            + d2
        )
        * (np.sqrt(n2o) - np.sqrt(n2o_base))
    ) * (forcing_scaling[..., n2o_indices])

    # linear for other gases
    # TODO: move to a general linear function
    erf_out[..., minor_greenhouse_gas_indices] = (
        (
            concentration[..., minor_greenhouse_gas_indices]
            - baseline_concentration[..., minor_greenhouse_gas_indices]
        )
        * radiative_efficiency[..., minor_greenhouse_gas_indices]
        * 0.001  # unit handling
    ) * (forcing_scaling[..., minor_greenhouse_gas_indices])

    return erf_out


def meinshausen2020(
    concentration,
    reference_concentration,
    forcing_scaling,
    radiative_efficiency,
    co2_indices,
    ch4_indices,
    n2o_indices,
    minor_greenhouse_gas_indices,
    a1=-2.4785e-07,
    b1=0.00075906,
    c1=-0.0021492,
    d1=5.2488,
    a2=-0.00034197,
    b2=0.00025455,
    c2=-0.00024357,
    d2=0.12173,
    a3=-8.9603e-05,
    b3=-0.00012462,
    d3=0.045194,
):
    """Greenhouse gas forcing from concentrations.

    This uses the [Meinshausen2020]_ relationships for CO2, CH4 and
    N2O including band overlaps between these three gases. This is a rescaled
    [Etminan2016]_ function with improved stability outside the
    range of validity in [Etminan2016]_. Minor greenhouse gases are a linear function of
    their concentration based on their radiative efficiency.

    Parameters
    ----------
    concentration : np.ndarray
        concentration of greenhouse gases. "CO2", "CH4" and "N2O" must be
        included in units of [ppm, ppb, ppb]. Other GHGs are units of ppt.
    reference_concentration : np.ndarray
        pre-industrial concentration of the gases (see above) used as the
        reference to calculate the forcing.
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : np.ndarray
        radiative efficiency to use for linear-forcing gases, in W m-2 ppb-1
    co2_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CO2.
    ch4_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CH4.
    n2o_indices : np.ndarray of bool
        index along SPECIES AXIS relating to N2O.
    minor_greenhouse_gas_indices : np.ndarray of bool
        indices of other GHGs that are not CO2, CH4 or N2O.
    a1 : float
        fitting parameter (see [Meinshausen2020]_)
    b1 : float
        fitting parameter (see [Meinshausen2020]_)
    c1 : float
        fitting parameter (see [Meinshausen2020]_)
    d1 : float
        fitting parameter (see [Meinshausen2020]_)
    a2 : float
        fitting parameter (see [Meinshausen2020]_)
    b2 : float
        fitting parameter (see [Meinshausen2020]_)
    c2 : float
        fitting parameter (see [Meinshausen2020]_)
    d2 : float
        fitting parameter (see [Meinshausen2020]_)
    a3 : float
        fitting parameter (see [Meinshausen2020]_)
    b3 : float
        fitting parameter (see [Meinshausen2020]_)
    d3 : float
        fitting parameter (see [Meinshausen2020]_)

    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W/m2) from greenhouse gases
    """
    erf_out = np.ones_like(concentration) * np.nan

    # easier to deal with smaller arrays
    co2 = concentration[..., co2_indices]
    ch4 = concentration[..., ch4_indices]
    n2o = concentration[..., n2o_indices]
    co2_base = reference_concentration[..., co2_indices]
    ch4_base = reference_concentration[..., ch4_indices]
    n2o_base = reference_concentration[..., n2o_indices]

    # CO2
    ca_max = co2_base - b1 / (2 * a1)
    where_central = np.asarray((co2_base < co2) & (co2 <= ca_max)).nonzero()
    where_low = np.asarray((co2 <= co2_base)).nonzero()
    where_high = np.asarray((co2 > ca_max)).nonzero()
    alpha_p = np.ones_like(co2) * np.nan
    alpha_p[where_central] = (
        d1
        + a1 * (co2[where_central] - co2_base[where_central]) ** 2
        + b1 * (co2[where_central] - co2_base[where_central])
    )
    alpha_p[where_low] = d1
    alpha_p[where_high] = d1 - b1**2 / (4 * a1)
    alpha_n2o = c1 * np.sqrt(n2o)
    erf_out[..., co2_indices] = (
        (alpha_p + alpha_n2o)
        * np.log(co2 / co2_base)
        * (forcing_scaling[..., co2_indices])
    )

    # CH4
    erf_out[..., ch4_indices] = (
        (a3 * np.sqrt(ch4) + b3 * np.sqrt(n2o) + d3)
        * (np.sqrt(ch4) - np.sqrt(ch4_base))
    ) * (forcing_scaling[..., ch4_indices])

    # N2O
    erf_out[..., n2o_indices] = (
        (a2 * np.sqrt(co2) + b2 * np.sqrt(n2o) + c2 * np.sqrt(ch4) + d2)
        * (np.sqrt(n2o) - np.sqrt(n2o_base))
    ) * (forcing_scaling[..., n2o_indices])

    # linear for other gases
    # TODO: move to a general linear function
    erf_out[..., minor_greenhouse_gas_indices] = (
        (
            concentration[..., minor_greenhouse_gas_indices]
            - reference_concentration[..., minor_greenhouse_gas_indices]
        )
        * radiative_efficiency[..., minor_greenhouse_gas_indices]
        * 0.001  # unit handling
    ) * (forcing_scaling[..., minor_greenhouse_gas_indices])

    return erf_out


def myhre1998(
    concentration,
    baseline_concentration,
    forcing_scaling,
    radiative_efficiency,
    co2_indices,
    ch4_indices,
    n2o_indices,
    minor_greenhouse_gas_indices,
    alpha_co2=5.35,
    alpha_ch4=0.036,
    alpha_n2o=0.12,
    alpha_ch4_n2o=0.47,
    a1=2.01e-5,
    exp1=0.75,
    a2=5.32e-15,
    exp2=1.52,
):
    """Greenhouse gas forcing from concentrations.

    Band overlaps are included between CH4 and N2O. This relationship comes from
    [Myhre1998]_, and was used up until the IPCC's Fifth Assessment
    Report. Minor greenhouse gases are a linear function of their concentration
    based on their radiative efficiency.

    Parameters
    ----------
    concentration : np.ndarray
        concentration of greenhouse gases. "CO2", "CH4" and "N2O" must be
        included in units of [ppm, ppb, ppb]. Other GHGs are units of ppt.
    baseline_concentration : np.ndarray
        pre-industrial concentration of the gases (see above).
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : np.ndarray
        radiative efficiency to use for linear-forcing gases, in W m-2 ppb-1
    co2_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CO2.
    ch4_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CH4.
    n2o_indices : np.ndarray of bool
        index along SPECIES AXIS relating to N2O.
    minor_greenhouse_gas_indices : np.ndarray of bool
        indices of other GHGs that are not CO2, CH4 or N2O.
    alpha_co2 : float
        factor relating logarithm of CO2 conentration to radiative forcing.
    alpha_ch4: float
        factor relating square root of CH4 conentration to radiative forcing.
    alpha_n2o : float
        factor relating square root of N2O conentration to radiative forcing.

    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W/m2) from greenhouse gases
    """

    def ch4_n2o_overlap(ch4, n2o, alpha_ch4_n2o, a1, exp1, a2, exp2):
        return alpha_ch4_n2o * np.log(
            1 + a1 * (ch4 * n2o) ** exp1 + a2 * ch4 * (ch4 * n2o) ** exp2
        )

    erf_out = np.ones_like(concentration) * np.nan

    # easier to deal with smaller arrays
    co2 = concentration[..., co2_indices]
    ch4 = concentration[..., ch4_indices]
    n2o = concentration[..., n2o_indices]
    co2_base = baseline_concentration[..., co2_indices]
    ch4_base = baseline_concentration[..., ch4_indices]
    n2o_base = baseline_concentration[..., n2o_indices]

    # CO2
    erf_out[..., co2_indices] = (
        alpha_co2 * np.log(co2 / co2_base) * (forcing_scaling[..., co2_indices])
    )

    # CH4
    erf_out[..., ch4_indices] = (
        alpha_ch4 * (np.sqrt(ch4) - np.sqrt(ch4_base))
        - (
            ch4_n2o_overlap(ch4, n2o_base, alpha_ch4_n2o, a1, exp1, a2, exp2)
            - ch4_n2o_overlap(ch4_base, n2o_base, alpha_ch4_n2o, a1, exp1, a2, exp2)
        )
    ) * forcing_scaling[..., ch4_indices]

    # N2O
    erf_out[..., n2o_indices] = (
        alpha_n2o * (np.sqrt(n2o) - np.sqrt(n2o_base))
        - (
            ch4_n2o_overlap(ch4_base, n2o, alpha_ch4_n2o, a1, exp1, a2, exp2)
            - ch4_n2o_overlap(ch4_base, n2o_base, alpha_ch4_n2o, a1, exp1, a2, exp2)
        )
    ) * forcing_scaling[..., n2o_indices]

    # linear for other gases
    # TODO: move to a general linear function
    erf_out[..., minor_greenhouse_gas_indices] = (
        (
            concentration[..., minor_greenhouse_gas_indices]
            - baseline_concentration[..., minor_greenhouse_gas_indices]
        )
        * radiative_efficiency[..., minor_greenhouse_gas_indices]
        * 0.001  # unit handling
    ) * (forcing_scaling[..., minor_greenhouse_gas_indices])

    return erf_out


def leach2021ghg(
    concentration,
    baseline_concentration,
    forcing_scaling,
    radiative_efficiency,
    co2_indices,
    ch4_indices,
    n2o_indices,
    minor_greenhouse_gas_indices,
    f1_co2=4.57,
    f3_co2=0.086,
    f3_ch4=0.038,
    f3_n2o=0.106,
):
    """Greenhouse gas forcing from concentrations.

    This uses the [Leach2021]_ relationships for CO2, CH4 and N2O
    that do not include band overlaps between these gases, allowing
    single-forcing runs. This is the default treatment in FaIR2.0. This is a
    re-fit of the [Etminan2016]_ formulation with improved
    coefficient fits. Minor greenhouse gases are a linear function of their
    concentration based on their radiative efficiency.

    Parameters
    ----------
    concentration : np.ndarray
        concentration of greenhouse gases. "CO2", "CH4" and "N2O" must be
        included in units of [ppm, ppb, ppb]. Other GHGs are units of ppt.
    baseline_concentration : np.ndarray
        pre-industrial concentration of the gases (see above).
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : np.ndarray
        radiative efficiency to use for linear-forcing gases, in W m-2 ppb-1
    co2_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CO2.
    ch4_indices : np.ndarray of bool
        index along SPECIES_AXIS relating to CH4.
    n2o_indices : np.ndarray of bool
        index along SPECIES AXIS relating to N2O.
    minor_greenhouse_gas_indices : np.ndarray of bool
        indices of other GHGs that are not CO2, CH4 or N2O.
    f1_co2 : float
        factor relating logarithm of CO2 conentration to radiative forcing.
    f3_co2 : float
        factor relating square root of CO2 conentration to radiative forcing.
    f3_ch4 : float
        factor relating square root of CH4 conentration to radiative forcing.
    f3_n2o : float
        factor relating square root of N2O conentration to radiative forcing.

    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W/m2) from greenhouse gases
    """
    erf_out = np.ones_like(concentration) * np.nan

    # easier to deal with smaller arrays
    co2 = concentration[..., co2_indices]
    ch4 = concentration[..., ch4_indices]
    n2o = concentration[..., n2o_indices]
    co2_base = baseline_concentration[..., co2_indices]
    ch4_base = baseline_concentration[..., ch4_indices]
    n2o_base = baseline_concentration[..., n2o_indices]

    # CO2
    erf_out[..., co2_indices] = (
        f1_co2 * np.log(co2 / co2_base) + f3_co2 * (np.sqrt(co2) - np.sqrt(co2_base))
    ) * forcing_scaling[..., co2_indices]

    # CH4
    erf_out[..., ch4_indices] = (
        f3_ch4 * (np.sqrt(ch4) - np.sqrt(ch4_base)) * forcing_scaling[..., ch4_indices]
    )

    # N2O
    erf_out[..., n2o_indices] = (
        f3_n2o * (np.sqrt(n2o) - np.sqrt(n2o_base)) * forcing_scaling[..., n2o_indices]
    )

    # linear for other gases
    # TODO: move to a general linear function
    erf_out[..., minor_greenhouse_gas_indices] = (
        (
            concentration[..., minor_greenhouse_gas_indices]
            - baseline_concentration[..., minor_greenhouse_gas_indices]
        )
        * radiative_efficiency[..., minor_greenhouse_gas_indices]
        * 0.001  # unit handling
    ) * (forcing_scaling[..., minor_greenhouse_gas_indices])

    return erf_out
