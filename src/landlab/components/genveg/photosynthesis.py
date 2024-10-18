import numpy as np


class Photosynthesis(object):
    def __init__(
        self,
        latitude,
        _current_day=0,
        _CO2_atmos=245,
        photo_params={
            "vcmax": 200,
            "kc": 300,
            "ko": 300000,
            "ci": 245,
            "co": 210000,
            "spec_factor_25": 2600.0,
            "stomatal_conductance": 1000000000,
        },
        gauss_integration_params=[
            (0.0469101, 0.1184635),
            (0.2307534, 0.2393144),
            (0.5000000, 0.2844444),
            (0.7692465, 0.2393144),
            (0.9530899, 0.1184635),
        ],
    ):
        self.latitude = latitude
        self._solar_declination = 0.0
        self._sunrise = 0.0
        self._sunset = 0.0
        self._sunlit_increment = 0.0
        self.gauss_integration_params = gauss_integration_params
        self.O2_max_coeff = photo_params["ko"]
        self.CO2_max_coeff = photo_params["kc"]
        self.O2_conc = photo_params["co"]
        self.init_CO2_conc = photo_params["ci"]
        self._CO2_atmos = _CO2_atmos
        self.Vc_max_rate = photo_params["vcmax"]
        self.spec_factor_base = photo_params["spec_factor_25"]
        self.stomatal_conductance = photo_params["stomatal_conductance"]
        self.update_solar_variables(_current_day)
        self.assim_limits_by_temp = self.calculate_assimilation_limits()

    def photosynthesize(
        self,
        grid_par_W_per_sqm,
        _min_temperature,
        _max_temperature,
        lai,
        last_biomass,
        _current_day,
    ):
        self.update_solar_variables(_current_day)
        total_canopy_assimilated_CO2 = hourly_gross_assimilation = gphot_CH20 = (
            np.zeros_like(last_biomass["shoot_sys_width"])
        )
        CO2_conc = self.init_CO2_conc
        for day_increment in self.gauss_integration_params:
            (abscissa, weight) = day_increment
            increment_hour = self._sunrise + abscissa * self._sunlit_increment
            leaf_temp_est = self.calculate_hour_temp(
                increment_hour, _min_temperature, _max_temperature
            )
            solar_elevation = self.calculate_solar_elevation(
                increment_hour
            ) * np.ones_like(grid_par_W_per_sqm)
            (
                absorbed_PAR_sunlit,
                absorbed_PAR_shaded,
            ) = self.calculate_absorbed_incremental_PAR(
                increment_hour,
                solar_elevation,
                grid_par_W_per_sqm,
                lai,
                _current_day,
            )
            sunlit_assimilated_CO2 = self.calculate_leaf_assimilation(
                absorbed_PAR_sunlit, CO2_conc, leaf_temp_est
            )
            shaded_assimilated_CO2 = self.calculate_leaf_assimilation(
                absorbed_PAR_shaded, CO2_conc, leaf_temp_est
            )
            (
                sunlit_lai,
                shaded_lai,
            ) = self.calculate_sunlit_shaded_lai_proportion(solar_elevation, lai)

            hourly_gross_assimilation = (
                sunlit_assimilated_CO2 * (sunlit_lai / lai)
            ) + (shaded_assimilated_CO2 * (shaded_lai / lai))
            total_canopy_assimilated_CO2 += (
                hourly_gross_assimilation * weight * 3600 * self._sunlit_increment
            )
            CO2_conc = (
                self._CO2_atmos - hourly_gross_assimilation / self.stomatal_conductance
            )
        gphot_CH20 = (
            total_canopy_assimilated_CO2 * last_biomass["live_leaf_area"] * 30 / 1000000
        )
        gphot_CH20[gphot_CH20 < 0] = 0.0
        return gphot_CH20

    def calculate_hour_temp(
        self,
        increment_hour,
        _min_temperature=25.792034827239835,
        _max_temperature=25.792034827239835,
    ):
        # We are assuming here that air temperature ~canopy temperature
        offset = 1.5
        temp_rise = self._sunrise + offset
        if (increment_hour >= temp_rise) & (increment_hour <= self._sunset):
            tau = np.pi * (increment_hour - temp_rise) / (self._sunset - self._sunrise)
            hour_temp = _min_temperature + (
                _max_temperature - _min_temperature
            ) * np.sin(tau)
        else:
            if increment_hour < temp_rise:
                increment_hour += 24
            tau = np.pi * (self._sunset - temp_rise) / (self._sunset - self._sunrise)
            sunset_temp = _min_temperature + (
                _max_temperature - _min_temperature
            ) * np.sin(tau)
            interp_slope = (_min_temperature - sunset_temp) / (
                temp_rise + 24 - self._sunset
            )
            hour_temp = sunset_temp + interp_slope * (increment_hour - self._sunset)
        return hour_temp

    def calculate_leaf_assimilation(
        self,
        par,
        CO2_conc,
        hour_temp,
    ):

        max_rubisco = self.calculate_rubisco_limits(CO2_conc, hour_temp)
        max_light_limit = self.calculate_light_limits(par, CO2_conc, hour_temp)
        max_sink_limit = self.get_sink_limits(hour_temp)
        min_assim = np.minimum(max_rubisco, max_light_limit)
        min_assim = np.minimum(min_assim, max_sink_limit)
        return min_assim

    def calculate_hourly_direct_light_extinction(self, solar_elevation):
        light_extinction = np.zeros_like(solar_elevation)
        light_extinction[solar_elevation > 0.00001] = 0.5 / np.sin(
            solar_elevation[solar_elevation > 0.00001]
        )
        return light_extinction

    def calculate_hourly_diffuse_light_extinction(self, lai):
        return (1.0 + 0.1174 * lai**0.5) / (1.0 + 0.3732 * lai**0.5)

    def update_solar_variables(self, _current_day):
        self._solar_declination = self.calculate_solar_declination(_current_day)
        self._sunset = self.calculate_solar_sunset(self._solar_declination)
        self._sunrise = self.calculate_solar_sunrise(self._sunset)
        self._sunlit_increment = self._sunset - self._sunrise

    def calculate_solar_declination(self, _current_jday):
        return -0.4093 * np.cos(2 * np.pi * (_current_jday + 10) / 365)

    def calculate_solar_sunset(self, _solar_declination):
        dA = np.sin(_solar_declination) * np.sin(self.latitude)
        dB = np.cos(_solar_declination) * np.cos(self.latitude)
        return 12 * np.arccos(-dA / dB) / np.pi + 12

    def calculate_solar_sunrise(self, sunset):
        return 24 - sunset

    def calculate_solar_elevation(self, increment_hour):
        dA = np.sin(self._solar_declination) * np.sin(self.latitude)
        dB = np.cos(self._solar_declination) * np.cos(self.latitude)
        dHa = np.pi * (increment_hour - 12) / 12
        return np.arcsin(dA + dB * np.cos(dHa))

    def calculate_absorbed_incremental_PAR(
        self, increment_hour, solar_elevation, grid_par_W_per_sqm, lai, current_day
    ):
        P = 0.04  # // reflection coefficient
        A = 0.80  # // scatter coefficient
        S = A**0.5  # // scatter correction
        conv = 4.55  # conversion factor to umol

        hourly_direct_PAR, hourly_diffuse_PAR = self.calculate_incremental_PAR(
            increment_hour, solar_elevation, 2 * grid_par_W_per_sqm, current_day
        )
        hourly_direct_light_extinction_k = (
            self.calculate_hourly_direct_light_extinction(solar_elevation)
        )
        hourly_diffuse_light_extinction_k = (
            self.calculate_hourly_diffuse_light_extinction(lai)
        )
        dI_direct = (1.0 - P) * hourly_direct_PAR
        dIpdr = dI_direct * np.exp(-hourly_direct_light_extinction_k * S * lai)
        dIpdrdr = dI_direct * np.exp(
            -hourly_direct_light_extinction_k * lai
        )  # direct of direct
        dIpdra = (dIpdr - dIpdrdr) / 2  # scatter beams
        dN = hourly_diffuse_light_extinction_k * S * lai
        dI_diffuse = (1.0 - P) * hourly_diffuse_PAR
        dIpdf = dI_diffuse * (1.0 - np.exp(-dN)) / dN  # diffuse
        absorbed_PAR_sunlit = (
            conv
            * A
            * (hourly_direct_light_extinction_k * hourly_direct_PAR + dIpdf + dIpdra)
        )
        absorbed_PAR_shaded = conv * A * (dIpdf + dIpdra)
        return (
            absorbed_PAR_sunlit,
            absorbed_PAR_shaded,
        )

    def calculate_hourly_ET_rad(self, solar_elevation, _current_day):
        E0 = 1 + 0.033 * np.cos(2 * np.pi * (_current_day - 10) / 365)
        ET_rad = 1370 * E0 * np.sin(solar_elevation)
        ET_rad[ET_rad < 0.0] = 0.0
        return ET_rad

    def calculate_incremental_PAR(
        self, increment_hour, solar_elevation, grid_par_W_per_sqm, _current_day
    ):
        dA = np.sin(self._solar_declination) * np.sin(self.latitude)
        dB = np.cos(self._solar_declination) * np.cos(self.latitude)
        dAoB = dA / dB
        dPhi = (np.pi * grid_par_W_per_sqm) / (
            dA * np.arccos(-dAoB) + dB * np.sqrt(1 - dAoB**2)
        )
        dCoefA = -dB * dPhi
        dCoefB = dA * dPhi
        total_incremental_rad = dCoefA * np.cos(np.pi * increment_hour / 12) + dCoefB
        total_ET_rad = self.calculate_hourly_ET_rad(solar_elevation, _current_day)
        R = (
            0.847
            - 1.61 * np.sin(solar_elevation)
            + 1.04 * (np.sin(solar_elevation)) ** 2
        )
        K = (1.47 - R) / 1.66
        transferred_rad = np.ones_like(grid_par_W_per_sqm)
        transferred_rad[total_ET_rad > 0.0] = (
            total_incremental_rad[total_ET_rad > 0.0] / total_ET_rad[total_ET_rad > 0.0]
        )
        condition_list = [
            transferred_rad > K,
            (transferred_rad <= K) & (transferred_rad > 0.35),
            (transferred_rad <= 0.35) & (transferred_rad > 0.22),
            (total_incremental_rad < 0.0) & (np.sin(solar_elevation) < 0),
        ]
        option_list = [
            R * np.ones_like(total_incremental_rad),
            1.47 - 1.66 * transferred_rad,
            1 - 6.4 * (transferred_rad - 0.22) ** 2,
            np.ones_like(total_incremental_rad),
        ]
        diffuse_frac = np.select(condition_list, option_list)
        total_incremental_rad[total_incremental_rad < 0.0] = 0.0
        total_incremental_par = 0.5 * total_incremental_rad
        diffuse_PAR = diffuse_frac * total_incremental_par
        direct_PAR = total_incremental_par - diffuse_PAR
        return (
            direct_PAR,
            diffuse_PAR,
        )

    def calculate_sunlit_shaded_lai_proportion(self, solar_elevation, lai):
        sunlit_lai = np.zeros_like(lai)
        hourly_direct_light_extinction_k = (
            self.calculate_hourly_direct_light_extinction(solar_elevation)
        )
        filter = np.nonzero(hourly_direct_light_extinction_k > 0.0)
        sunlit_lai[filter] = (
            1.0 - np.exp(-hourly_direct_light_extinction_k[filter] * lai[filter])
        ) / hourly_direct_light_extinction_k[filter]
        shaded_lai = lai - sunlit_lai
        return (sunlit_lai, shaded_lai)

    def calculate_light_limits(self, par, CO2_conc, hour_temp):
        quantum_yield = 0.06
        absorption_frac = 0.8
        CO2_comp = self.get_CO2_comp(hour_temp)
        dA = absorption_frac * quantum_yield * par * (CO2_conc - CO2_comp)
        dB = CO2_conc + 2 * CO2_comp
        return dA / dB

    def calculate_rubisco_limits(self, CO2_conc, leaf_temp):
        CO2_comp = self.get_CO2_comp(leaf_temp)
        O2_coeff = self.O2_max_coeff * 1.2 ** ((leaf_temp - 25) / 10)  # 142.86
        CO2_coeff = self.CO2_max_coeff * 2.1 ** ((leaf_temp - 25) / 10)
        Vc_denom = 1 + np.exp(0.128 * (leaf_temp - 40))
        Vc_adj = self.Vc_max_rate * 2.4 ** ((leaf_temp - 25) / 10) / Vc_denom
        dcoeffm = CO2_coeff * (1 + self.O2_conc / O2_coeff)
        dA = Vc_adj * (CO2_conc - CO2_comp)
        dB = CO2_conc + dcoeffm
        rubisco_limits = dA / dB
        return rubisco_limits

    def calculate_assimilation_limits(self):
        # this needs to happen at init then have leaf temp interpolated
        leaf_temp = np.arange(-50, 50, 0.25)
        Vc_denom = 1 + np.exp(0.128 * (leaf_temp - 40))
        Vc_adj = self.Vc_max_rate * 2.4 ** ((leaf_temp - 25) / 10) / Vc_denom
        CO2_comp = (
            0.5
            * self.O2_conc
            / (self.spec_factor_base * 0.57 ** ((leaf_temp - 25) / 10))
        )
        sink_limits = 0.5 * Vc_adj
        dtypes = [
            ("leaf_temp", float),
            ("sink_limits", float),
            ("CO2_comp", float),
        ]

        limit_map = np.column_stack((leaf_temp, sink_limits, CO2_comp))
        limit_map = list(map(tuple, limit_map))
        limit_lookup = np.array(limit_map, dtypes)
        return limit_lookup

    def get_sink_limits(self, hour_temp):
        limits = np.interp(
            hour_temp,
            self.assim_limits_by_temp["leaf_temp"],
            self.assim_limits_by_temp["sink_limits"],
        )
        return limits

    def get_CO2_comp(self, hour_temp):
        comp_pt = np.interp(
            hour_temp,
            self.assim_limits_by_temp["leaf_temp"],
            self.assim_limits_by_temp["CO2_comp"],
        )
        return comp_pt


class C3(Photosynthesis):
    def __init__(
        self,
        latitude,
        photo_params={
            "vcmax": 200,
            "kc": 300,
            "ko": 300000,
            "ci": 245,
            "co": 210000,
            "spec_factor_base": 2600.0,
        },
    ):
        super().__init__(latitude, photo_params=photo_params)


class C4(Photosynthesis):
    def __init__(
        self,
        latitude,
        _CO2_atmos=400,
        photo_params={
            "vcmax": 200,
            "kc": 300,
            "ko": 300000,
            "ci": 245,
            "co": 210000,
            "spec_factor_base": 2600.0,
            "stomatal_conductance": 0.5,
        },
    ):
        super().__init__(latitude, _CO2_atmos=_CO2_atmos, photo_params=photo_params)
        # This will be updated with C4 changes


class Cam(Photosynthesis):
    def __init__(self, latitude, photo_params):
        super().__init__(latitude, photo_params=photo_params)
