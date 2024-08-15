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
    fpath: Pathfile object, file path where input file is located.
           Must be structured Excel or csv file or yaml.
           If blank, model will assume species is corn.
    processes: a list of vegetation processes to initialize parameters for
               'plantsize','dispersal','mortality','colonization'.
               If blank, model will run basic growth only.
    outfile: optional input to allow for custom file name
    """

    def __init__(self, fpath="None", outfile="veg_params.yml", vegparams={}):
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
                    "duration_params": {
                        "growing_season_start": 91,
                        "growing_season_end": 290,
                        "peak_biomass": 227,
                        "reproduction_start": 180,
                        "reproduction_end": 227,
                        "senescence_start": 228,
                        "max_age": 365,
                    },
                    "photo_params": {
                        "vcmax": 200,
                        "kc": 30000,
                        "ko": 300,
                        "ci": 245,
                        "co": 209,
                        "spec_factor_25": 2600,
                        "stomatal_conductance": 30000,
                    },
                    "grow_params": {
                        "respiration_coefficient": {
                            "root": 0.015,
                            "leaf": 0.015,
                            "stem": 0.03,
                            "reproductive": 0.01,
                        },
                        "glucose_requirement": {
                            "root": 1.444,
                            "leaf": 1.513,
                            "stem": 1.463,
                            "reproductive": 1.414,
                        },
                        "root_to_leaf": [0.031, 0.951, 0],
                        "root_to_stem": [-0.107, 1.098, 0.0216],
                        "plant_part_min": {
                            "root": 0.01,
                            "leaf": 0.1,
                            "stem": 0.5,
                            "reproductive": 0.0,
                        },
                        "plant_part_max": {
                            "root": 6,
                            "leaf": 25,
                            "stem": 30,
                            "reproductive": 200,
                        },
                    },
                    "morph_params": {
                        "max_plant_density": 1,
                        "max_n_stems": 3,
                        "max_height": 2.5,
                        "max_shoot_sys_width": 0.5,
                        "max_root_root_sys_depth": 0.33,
                        "min_height": 0.075,
                        "min_shoot_sys_width": 0.01,
                        "min_root_sys_width": 0.01,
                        "sp_leaf_area": 0.02,
                        "biomass_decay_rate": 0.07,
                        "lai_cr": 4,
                    },
                    "dispersal_params": {
                        "max_dist_dispersal": 2,
                        "min_size_dispersal": 0.5,
                        "unit_cost_dispersal": 1.2,
                    },
                    "col_params": {
                        "prob_colonization": 0.01,
                        "time_to_colonization": 365,
                    },
                    "mortality_params": {
                        "mort_variable_name": {1: "Mortality factor"},
                        "duration": {1: 365},
                        "period": {1: "during growing season"},
                        "predictor": {1: [0, 0]},
                        "response": {1: [0, 0]},
                    },
                }
            }
        else:
            ispathvalid = fpath.is_file()
            if not ispathvalid:
                raise ValueError("File path is not valid")
            self.fpath = fpath
            # add check for file extension
            exten = pathlib.Path(self.fpath).suffix
            if exten == "yml":
                print(
                    "File already in correct file format."
                    "Use Landlab load_params function."
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
                        # Parse data from each sheet, replace NaNs with -9999,
                        # and group duplicate index values to a list
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
                                # Create new dictionary level for the parameter group
                                # and take a cross-section of the dataframe
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
                                        # Format variable values & assign them directly
                                        # to the variable name or with descriptor labels
                                        # for multi-part variables
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

                        # Replace null values for coefficients with
                        # Poorter-derived coefficients if necessary
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
                        descriptor_wout_data = []
                        for descriptor in nested_dict["mortality_params"][
                            "mort_variable_name"
                        ]:
                            var_name = nested_dict["mortality_params"][
                                "mort_variable_name"
                            ][descriptor]
                            if var_name == -9999:
                                descriptor_wout_data.append(descriptor)
                        for var_dict_name in nested_dict["mortality_params"]:
                            for descriptor in descriptor_wout_data:
                                del nested_dict["mortality_params"][var_dict_name][
                                    descriptor
                                ]
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
                        param_dict[nested_dict["plant_factors"]["species"]] = (
                            nested_dict
                        )
                else:
                    if exten == "csv":
                        # Add Carra's code here and load into dict called x
                        pass
                    else:
                        raise ValueError("File extension not recognized")

            self.vegparams = param_dict

        with open(outfile, "w") as outfile:
            yaml.dump(self.vegparams, outfile, default_flow_style=True)

    # Private method to build logistic mortality
    # function for up to five acute mortality factors
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
            msg = (
                "Not enough points to generate logistic function."
                "Assuming zero mortality."
            )
            print(msg)
            S = np.array([0, 0])
        # Direct solve for coefficients if only two
        # points provided to prevent solver errors
        elif len(uni_xs) == 2:
            # Solve for constant (see function below)
            guess = self._a_b_func(uni_xs, uni_ys)
            S = np.array([guess["a"], guess["b"]])
        # Use scipy solver to estimate sigmoid coefficients
        else:
            S = self._get_general_s_curve(uni_xs, uni_ys, fit_method=fit_method)

        return S

    def _get_general_s_curve(self, uni_xs, uni_ys, fit_method):
        """
        This function yeilds a general S curve based on the _cfunc method. It
        utilizes uni_xs, and uni_ys values. This method will utilizes SciPy's curve_fit
        function. It is assumed that using the points closest to 0.5 survival and near
        the sigmoid limit finds the best gerneral s curve. Weights are also used to
        priortize points near 0 and 1.  If there is only one point that is both the 0.5
        survival and near the sigmoid limit the TAM method will be used to create more
        points and then the curve_fit function will be used to obtain and a and b
        estimation
        """

        # find the index for the 0.5 survival and sigmoid limit
        idx_05, idx_limit = self.find_index_values(uni_xs=uni_xs, uni_ys=uni_ys)
        # check to see if they are the same point
        if idx_05 != idx_limit:
            # get x and y coordinates
            x = [uni_xs[min(idx_05, idx_limit)], uni_xs[max(idx_05, idx_limit)]]
            y = [uni_ys[min(idx_05, idx_limit)], uni_ys[max(idx_05, idx_limit)]]
            # get initial a and b guess
            initial_guess = self._a_b_func(x, y)
            # get weights to prioritize points near 0 and 1
            weights = self.get_weights(uni_xs)
            # get the a and be estimation from curve_fit
            S, _ = curve_fit(
                self._cfunc,
                uni_xs,
                uni_ys,
                p0=[initial_guess["a"], initial_guess["b"]],
                sigma=weights,
                method=fit_method,
            )
        # goes to TAM method to make the curve
        else:
            S = self._TAM_method(uni_xs=uni_xs, uni_ys=uni_ys, method=fit_method)

        return S

    def _mse(self, x, y, coeffs):
        # currently this is not used but will be when we produce the graphs and show
        # a warning message for bad fit
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
        b_guess = -(np.log((1 - y[1]) / y[1]) - np.log((1 - y[0]) / y[0])) / (
            x[1] - x[0]
        )
        return {"b": b_guess, "a": ((1 - y[1]) / y[1]) / np.exp(-x[1] * b_guess)}

    def _TAM_method(self, uni_xs, uni_ys, method="dogbox"):
        """
        This will be used when the data List of priority for choosing a pair of (x, y)
        points is the uni_ys point closest to 0.5 survival and near sigmoid limit is
        the same point. This is another method to estimate the build logistic mortality
        function for up to five acute mortality This is a python updated version of the
        logistic function from
        "Ecological Model Development: Toolkit for interActive Modeling (TAM)"
        by Carrillo, Carra C. et al. 2022. Paper suggested "The simplest approach
        is to identify the points where the index value is approximately 1 and where
        the index value is approximately 0" and so we will be using
        the unique point where unique y is closest 0 and 1
        """

        # Get the index value for the min and max data point of uni_xs and uni_ys
        min_val_idx = np.argmin(uni_ys)
        max_val_idx = np.argmax(uni_ys)

        x_1 = uni_xs[min_val_idx]
        x_2 = uni_xs[max_val_idx]
        y_1 = uni_ys[min_val_idx]
        y_2 = uni_ys[max_val_idx]

        # smooth uni_xs values for a fine estimation (note this will create 50
        # points inclusively between x1 and x2.)
        uni_xs_range = np.linspace(x_1, x_2)

        # Equations found in index of referenced paper above
        G = np.log(y_1 / (1 - y_1))
        F = np.log(y_2 / (1 - y_2))
        B = (G - F) / (x_1 - x_2)
        A = G - B * x_1
        Z = np.exp(A + B * uni_xs_range)
        si_x = Z / (1 + Z)

        # Optaion the a, b, parameter from _cfunc
        tam_weights = self.get_weights(si_x)

        # get point values closes to 0.5 survival and near sigmoid limit for TAM curve
        tam_idx_05, tam_idx_limit = self.find_index_values(
            uni_xs=uni_xs_range, uni_ys=si_x
        )
        x = [
            uni_xs_range[min(tam_idx_05, tam_idx_limit)],
            uni_xs_range[max(tam_idx_05, tam_idx_limit)],
        ]
        y = [si_x[min(tam_idx_05, tam_idx_limit)], si_x[max(tam_idx_05, tam_idx_limit)]]

        # a, b guess from TAM curve
        tam_a_b = self._a_b_func(x, y)

        # a, b estimation
        tam_S, _ = curve_fit(
            self._cfunc,
            uni_xs_range,
            si_x,
            p0=[tam_a_b["a"], tam_a_b["b"]],
            sigma=tam_weights,
            method=method,
        )
        return tam_S

    def get_weights(self, y_vals):
        # Assign sigma weights to prioritize points near 0 and 1
        weights = np.ones(len(y_vals))
        weights[(y_vals < 0.1) | (y_vals > 0.9)] = 10
        weights[(y_vals < 0.02) | (y_vals > 0.98)] = 500

        return weights

    def find_index_values(self, uni_xs, uni_ys):
        # Return the index value of the point closest to the 0.5 survival
        # and near sigmoid limit
        idx_05 = (np.abs(uni_ys - 0.5)).argmin()
        idx_limit = np.gradient(np.gradient(uni_ys, uni_xs), uni_xs).argmax()
        return idx_05, idx_limit

    def _cfunc(self, x, a, b):
        return 1 / (1 + a * np.exp(-b * x))
