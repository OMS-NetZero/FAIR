import os

import pandas as pd

gas_params_filename = os.path.join(
    os.path.dirname(__file__), "default_gas_parameters.csv"
)
gas_params_df = pd.read_csv(gas_params_filename, skiprows=1, index_col="Gas").T


def get_gas_params(gas_list=None):
    """Get gas parameters from the defauls.

    Parameters
    ----------
    gas_list : list[str] or None
        list of gases to calculate. If None, use all defaults

    Returns
    -------
    gas_params_df : :obj:`pandas.DataFrame`
        data frame of requested gases
    """
    res = pd.DataFrame(
        index=[
            "a1",
            "a2",
            "a3",
            "a4",
            "tau1",
            "tau2",
            "tau3",
            "tau4",
            "r0",
            "rC",
            "rT",
            "rA",
            "PI_conc",
            "emis2conc",
            "f1",
            "f2",
            "f3",
            "aer_conc",
        ]
    )
    res.columns.name = "Forcing"

    if gas_list is not None:
        for gas in gas_list:
            forcings = gas_params_df.filter(
                regex=r"^"
                + gas.replace("|", "\|")  # noqa: W605
                + "$"
                + "|"
                + "^"
                + gas.replace("|", "\|")  # noqa: W605
                + "\|"  # noqa: W605
            )
            if forcings.shape[1] != 0:
                for forcing in forcings.columns:
                    res[forcing] = gas_params_df[forcing]
            elif "|" in gas:
                sub_gas = gas[0 : gas.index("|")]
                sub_gas_params = get_gas_params[sub_gas]
                for forcing in sub_gas_params.columns:
                    forcing_string = forcing[forcing.index(sub_gas) + len(sub_gas) :]
                    res[gas + forcing_string] = sub_gas_params[forcing]
    else:
        res = gas_params_df
    return res
