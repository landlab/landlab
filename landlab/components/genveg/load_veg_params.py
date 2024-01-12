"""
Load veg veg_params
May need to change this to output a landlab configured input file
Want to retain formatted Excel file due to type of data required
"""
import numpy as np
import pandas as pd
import pathlib
import yaml
from scipy.optimize import curve_fit


# Create veg_params class type for input data and calculated params
class VegParams:
    """
    Load vegetation parameters from formatted Excel file
    into a Pandas dateframe that is converted to parameter
    dictionaries stored in a yaml for landlab input

    Parameters
    ----------
    fpath: Pathfile object, file path where input file is located. Must be structured Excel or csv file or yaml.
           If blank, model will assume species is corn.
    processes: a list of vegetation processes to initialize parameters for 'plantsize','dispersal','mortality','colonization'.
               If blank, model will run basic growth only.
    outfile: optional input to allow for custom file name
    """

    def __init__(
        self, fpath="None", outfile="veg_params.yml", processes=[], vegparams={}
    ):
        if fpath == "None":
            self.veg_params = {
                "Corn": {
                    "plant_factors": {
                        "species": "Corn",
                        "growth_habit": "forb_herb",
                        "monocot_dicot": "monocot",
                        "angio_gymno": "angiosperm",
                        "annual_perennial": "annual",
                        "leaf_retention": "deciduous",
                        "growth_form": "single_stem",
                        "shape": "erect",
                        "ptype": "C3",
                    },
                    "grow_params": {
                        "growing_season_start": 91,
                        "growing_season_end": 290,
                        "senescence_start": 228,
                        "respiration_coefficient": [0.015, 0.015, 0.03],
                        "glucose_requirement": [1.444, 1.513, 1.463],
                        "k_light_extinct": 0.02,
                        "light_half_sat": 9,
                        "p_max": 0.055,
                        "root_to_leaf_coeffs": [0.031, 0.951, 0],
                        "root_to_stem_coeffs": [-0.107, 1.098, 0.0216],
                        "plant_part_min": [0.01, 0.1, 0.5],
                        "plant_part_max": [6, 25, 30],
                    },
                }
            }
            if "plantsize" in processes:
                self.size_params = {
                    "max_plant_density": 1,
                    "max_n_stems": 3,
                    "max_height_stem": 2.5,
                    "max_mass_stem": 72,
                    "total_cs_area_stems": 0.231,
                }
                if "dispersion" in processes:
                    self.disp_params = {
                        "max_dist_dispersal": 2,
                        "disp_size_rat": 0.5,
                        "disp_cost": 0,
                    }
                else:
                    self.disp_params = {}
                self.veg_params["Corn"]["dispersal_params"] = {**self.dispersal_params}
            else:
                self.size_params = {}
            self.veg_params["Corn"]["size_params"] = {**self.size_params}
            if "colonize" in processes:
                self.col_params = {"col_prob": 0.01, "col_dt": 365}
            else:
                self.col_params = {}
            self.veg_params["Corn"]["col_params"] = {**self.col_params}
            if "mortality" in processes:
                self.mort_params = {
                    "mort_factor_1": "Mortality factor",
                    "mort_factor_1_duration": 365,
                    "mort_factor_1_coeffs": [0, 0],
                }
            else:
                self.mort_params = {}
            self.vegparams["Corn"]["mort_params"] = {**self.mort_params}
        else:
            ispathvalid = fpath.is_file()
            if not ispathvalid:
                raise ValueError("File path is not valid")
            self.fpath = fpath
            # add check for file extension
            exten = pathlib.Path(self.fpath).suffix
            if exten == "yml":
                print(
                    "File already in correct file format. Use Landlab load_params function."
                )
                pass
            else:
                if "xls" in exten:
                    # Read Excel file data into dataframe
                    xlin = pd.ExcelFile(fpath)
                    species = xlin.sheet_names
                    param_dict = {}
                    # Generate nested parameter dictionary for each species
                    for i in species:
                        # Parse data from each sheet, replace NaNs with -9999, and group duplicate index values to a list
                        nested_dict = {}
                        df_in = xlin.parse(i, usecols="B,C,D,E")
                        df_in.fillna(-9999, inplace=True)
                        df_grouped = df_in.groupby(
                            ["Process", "Variable Name", "Descriptor"], dropna=False
                        ).agg(pd.Series.tolist)
                        # Get first level of multiindex to determine parameter groups
                        var_groups = df_grouped.index.levels[0].values
                        for group in var_groups:
                            if group == -9999:
                                continue
                            else:
                                # Create new dictionary level for the parameter group and take a cross-section of the dataframe
                                nested_dict[group] = {}
                                group_vars = df_grouped.xs(group)
                                # Get list of variable names in the parameter group
                                var_names = group_vars.index.get_level_values(0).values
                                for var in var_names:
                                    # Create new dictionary level for each variable
                                    nested_dict[group][var] = {}
                                    var_vals = group_vars.xs(var)
                                    # Get descriptor names for multi-part variables
                                    descriptor_names = var_vals.index.get_level_values(
                                        0
                                    ).values
                                    for descriptor in descriptor_names:
                                        # Format variable values and assign them directly to the variable name or with descriptor labels for multi-part variables
                                        entry = var_vals["Values"].loc[descriptor]
                                        entry = (lambda x: x[0] if len(x) == 1 else x)(
                                            entry
                                        )
                                        # All simple variables have no descriptor
                                        if descriptor == -9999:
                                            nested_dict[group][var] = entry
                                        else:
                                            nested_dict[group][var][descriptor] = entry

                        # Calculate derived parameters
                        nested_dict["grow_params"]["total_min_biomass"] = sum(
                            nested_dict["grow_params"]["plant_part_min"].values()
                        )
                        nested_dict["grow_params"]["total_max_biomass"] = sum(
                            nested_dict["grow_params"]["plant_part_max"].values()
                        )
                        nested_dict["grow_params"]["growth_min_biomass"] = (
                            nested_dict["grow_params"]["total_min_biomass"]
                            - nested_dict["grow_params"]["plant_part_min"][
                                "reproductive"
                            ]
                        )
                        nested_dict["grow_params"]["growth_max_biomass"] = (
                            nested_dict["grow_params"]["total_max_biomass"]
                            - nested_dict["grow_params"]["plant_part_max"][
                                "reproductive"
                            ]
                        )

                        # Replace null values for coefficients with Poorter-derived coefficients if necessary
                        woody_herb = ("herb", "woody")[
                            nested_dict["plant_factors"]["growth_habit"] == "shrub"
                        ]
                        opt_2 = (
                            nested_dict["plant_factors"]["monocot_dicot"],
                            nested_dict["plant_factors"]["angio_gymno"],
                        )[woody_herb == "woody"]
                        options = {
                            "woody": {
                                "angiosperm": {
                                    "root_to_leaf": {
                                        "a": 0.090,
                                        "b1": 0.889,
                                        "b2": -0.0254,
                                    },
                                    "root_to_stem": {
                                        "a": -0.097,
                                        "b1": 1.071,
                                        "b2": 0.0179,
                                    },
                                    "max_nsc_content": {
                                        "root": 0.13607,
                                        "leaf": 0.15300,
                                        "stem": 0.12857,
                                        "reproductive": 0.13607,
                                    },
                                    "nsc_content": {
                                        "root": 0.06750,
                                        "leaf": 0.11815,
                                        "stem": 0.06321,
                                        "reproductive": 0.06750,
                                    },
                                    "min_nsc_content": {
                                        "root": 0.00643,
                                        "leaf": 0.04821,
                                        "stem": 0.01286,
                                        "reproductive": 0.00643,
                                    },
                                    "incremental_nsc": {
                                        "root": [0.5, -0.75, 0, 0.25],
                                        "leaf": [-1, 0.25, 0.75, -1],
                                        "stem": [0.5, -0.75, 0, 0.5],
                                        "reproductive": [
                                            0.625,
                                            -0.5625,
                                            0.0625,
                                            0.3125,
                                        ],
                                    },
                                },
                                "gymnosperm": {
                                    "root_to_leaf": {
                                        "a": 0.243,
                                        "b1": 0.924,
                                        "b2": -0.0282,
                                    },
                                    "root_to_stem": {
                                        "a": -0.070,
                                        "b1": 1.236,
                                        "b2": -0.0186,
                                    },
                                    "max_nsc_content": {
                                        "root": 0.06107,
                                        "leaf": 0.36629,
                                        "stem": 0.07286,
                                        "reproductive": 0.06107,
                                    },
                                    "nsc_content": {
                                        "root": 0.02893,
                                        "leaf": 0.14792,
                                        "stem": 0.02571,
                                        "reproductive": 0.02893,
                                    },
                                    "min_nsc_content": {
                                        "root": 0.01071,
                                        "leaf": 0.05714,
                                        "stem": 0.00750,
                                        "reproductive": 0.01071,
                                    },
                                    "incremental_nsc": {
                                        "root": [-0.25, 0, 0.25, 0.5],
                                        "leaf": [-0.5, 1.5, -1.25, -0.75],
                                        "stem": [0.5, 0.25, 0, -0.5],
                                        "reproductive": [
                                            -0.1875,
                                            0.0625,
                                            0.3125,
                                            0.625,
                                        ],
                                    },
                                },
                            },
                            "herb": {
                                "monocot": {
                                    "root_to_leaf": {"a": 0.031, "b1": 0.951, "b2": 0},
                                    "root_to_stem": {
                                        "a": -0.107,
                                        "b1": 1.098,
                                        "b2": 0.0216,
                                    },
                                    "max_nsc_content": {
                                        "root": 0.36643,
                                        "leaf": 0.36629,
                                        "stem": 0.30964,
                                        "reproductive": 0.36643,
                                    },
                                    "nsc_content": {
                                        "root": 0.21429,
                                        "leaf": 0.17396,
                                        "stem": 0.11143,
                                        "reproductive": 0.21429,
                                    },
                                    "min_nsc_content": {
                                        "root": 0.01071,
                                        "leaf": 0.01548,
                                        "stem": 0.00750,
                                        "reproductive": 0.01071,
                                    },
                                    "incremental_nsc": {
                                        "root": [1.25, -2.5, 0, 2],
                                        "leaf": [1.25, 0, -1, 0.5],
                                        "stem": [0, -0.5, 0, 0.5],
                                        "reproductive": [1.5625, -1.875, 0.0625, 2.5],
                                    },
                                },
                                "dicot": {
                                    "root_to_leaf": {"a": 0.259, "b1": 0.916, "b2": 0},
                                    "root_to_stem": {"a": -0.111, "b1": 1.029, "b2": 0},
                                    "max_nsc_content": {
                                        "root": 0.36643,
                                        "leaf": 0.36629,
                                        "stem": 0.30964,
                                        "reproductive": 0.36643,
                                    },
                                    "nsc_content": {
                                        "root": 0.21429,
                                        "leaf": 0.17396,
                                        "stem": 0.11143,
                                        "reproductive": 0.21429,
                                    },
                                    "min_nsc_content": {
                                        "root": 0.01071,
                                        "leaf": 0.01548,
                                        "stem": 0.00750,
                                        "reproductive": 0.01071,
                                    },
                                    "incremental_nsc": {
                                        "root": [1.25, -2.5, 0, 2],
                                        "leaf": [1.25, 0, -1, 0.5],
                                        "stem": [0, -0.5, 0, 0.5],
                                        "reproductive": [1.5625, -1.875, 0.0625, 2.5],
                                    },
                                },
                            },
                        }
                        df_fill = options[woody_herb][opt_2]

                        replace_vars = [
                            "root_to_leaf",
                            "root_to_stem",
                            "min_nsc_content",
                            "nsc_content",
                            "max_nsc_content",
                            "incremental_nsc",
                        ]
                        for var in replace_vars:
                            if (
                                list(nested_dict["grow_params"][var].values())[0]
                                == -9999
                            ):
                                nested_dict["grow_params"][var] = df_fill[var]

                        # Calculate sigmoid mortality curves
                        sigmoid_coeffs = {}
                        for descriptor in nested_dict["mortality_params"]["response"]:
                            response = np.array(
                                nested_dict["mortality_params"]["response"][descriptor]
                            )
                            predictor = np.array(
                                nested_dict["mortality_params"]["predictor"][descriptor]
                            )
                            response = response[response != -9999]
                            predictor = predictor[predictor != -9999]
                            sigmoid_coeffs[descriptor] = self._build_logistic(
                                predictor, response, fit_method="dogbox"
                            ).tolist()
                        nested_dict["mortality_params"]["coeffs"] = sigmoid_coeffs
                        # Add species nested dictionary to master parameter dictionary
                        param_dict[
                            nested_dict["plant_factors"]["species"]
                        ] = nested_dict
                else:
                    if exten == "csv":
                        # Add Carra's code here and load into dict called x
                        pass
                    else:
                        raise ValueError("File extension not recognized")

            self.vegparams = param_dict

        with open(outfile, "w") as outfile:
            yaml.dump(self.vegparams, outfile, default_flow_style=True)

    # Private method to build logistic mortality function for up to five acute mortality factors
    def _build_logistic(self, xs, ys, fit_method):
        if len(xs) != len(ys):
            msg = "Predictor and response variable arrays must be same length"
            raise ValueError(msg)
        ys[ys <= 0] = 0.0001
        ys[ys >= 1] = 0.9999
        # Keep only unique y values and average x values for repeated ys
        uni_ys, inverses = np.unique(ys, return_inverse=True)
        uni_xs = np.zeros_like(uni_ys)
        for i in range(len(uni_ys + 1)):
            uni_xs[i] = np.mean(xs[inverses == i])
        # Assume no mortality if not enough data provided for sigmoid curve estimate
        if len(uni_xs) <= 1:
            msg = "Not enough points to generate logistic function. Assuming zero mortality."
            print(msg)
            S = [0, 0]
        # Direct solve for coefficients if only two points provided to prevent solver errors
        elif len(uni_xs) == 2:
            # Solve for constant (see function below)
            b = -(
                    np.log((1 - uni_ys[1]) / uni_ys[1])
                    - np.log((1 - uni_ys[0]) / uni_ys[0])
            ) / (uni_xs[1] - uni_xs[0])
            a = ((1 - uni_ys[1]) / uni_ys[1]) / np.exp(-uni_xs[1] * b)
            S = [a, b]
        # Use scipy solver to estimate sigmoid coefficients
        else:
            S = self._get_best_a_b_guess(uni_xs, uni_ys)
        return S

    def _get_best_a_b_guess(self, uni_xs, uni_ys):
        """
        This function finds the best guess for 'a' and 'b' for _cfunc based on the values uni_xs
        and uni_ys values. List of priority for choosing a pair of (x, y) points is the uni_ys
        point closest to 0.5 survival and near sigmoid limit and then finding the best fit best
        fit using min/max points, min point/0.5 survival, 0.5 survival / max point, no intial guess
        inputs:
         - uni_xs: unique xs values
         - uni_ys: unique ys values
        outputs:
         - S: vector of a, b values
        """
        idx_05 = (np.abs(uni_ys - 0.5)).argmin()
        idx_limit = np.gradient(np.gradient(uni_ys, uni_xs), uni_xs).argmax()

        # Assign sigma weights to prioritize points near 0 and 1
        weights = np.ones(len(uni_ys))
        weights[(uni_ys < 0.1) | (uni_ys > 0.9)] = 10
        weights[(uni_ys < 0.02) | (uni_ys > 0.98)] = 500

        if idx_05 != idx_limit:
            x = [uni_xs[min(idx_05, idx_limit)], uni_xs[max(idx_05, idx_limit)]]
            y = [uni_ys[min(idx_05, idx_limit)], uni_ys[max(idx_05, idx_limit)]]
            guess = self._a_b_func(x, y)
            S, pcov = curve_fit(
                self._cfunc,
                uni_xs,
                uni_ys,
                p0=[guess['a'], guess['b']],
                sigma=weights,
                method="dogbox"
            )
        else:
            max_val_indx = np.argmax(uni_ys)
            min_val_indx = np.argmin(uni_ys)

            guess_min = None
            guess_max = None
            guess_min_max = None
            if min_val_indx != idx_05:
                x = [uni_xs[min_val_indx], uni_xs[idx_05]]
                y = [uni_ys[min_val_indx], uni_ys[idx_05]]
                guess_min = self._a_b_func(x, y)
            if max_val_indx != idx_05:
                x = [uni_xs[idx_05], uni_xs[max_val_indx]]
                y = [uni_ys[idx_05], uni_ys[max_val_indx]]
                guess_max = self._a_b_func(x, y)
            if min_val_indx != max_val_indx:
                x = [uni_xs[min_val_indx], uni_xs[max_val_indx]]
                y = [uni_ys[min_val_indx], uni_ys[max_val_indx]]
                guess_min_max = self._a_b_func(x, y)

            S_vals = []
            mse_vals = []
            for guess in [None, guess_min, guess_max, guess_min_max]:
                if guess != None:
                    p0 = [guess['a'], guess['b']]
                else:
                    p0 = None
                S_temp, pcov = curve_fit(self._cfunc, uni_xs, uni_ys, p0=p0, sigma=weights, method="dogbox")
                S_vals.append(S_temp)
                mse_vals.append(self.mse(uni_xs, uni_ys, S_temp))

            lowest_mse_idx = np.argmin(mse_vals)

            S = S_vals[lowest_mse_idx]

        return S

    def mse(self, x, y, coeffs):
        return np.mean((self._cfunc(x, *coeffs) - y) ** 2)

    def _a_b_func(self, x, y):
        """
        a and b guess values based on a pair of x/y points ((x_0, y_0), (x_1, y_1))
        inputs
        - x: vector of x points [x_0, x_1]
        - y: vector of y points [y_0, y_1]
        outputs:
        (dict) of a, b guess
        """
        b_guess = -(np.log((1 - y[1]) / y[1]) - np.log((1 - y[0]) / y[0])) / (x[1] - x[0])
        return {
            'b': b_guess,
            'a': ((1 - y[1]) / y[1]) / np.exp(-x[1] * b_guess)
        }

    def _cfunc(self, x, a, b):
        return 1 / (1 + a * np.exp(-b * x))
