import numpy as np


class Photosynthesis(object):
    def __init__(self, latitude, _current_day=0):
        self.latitude = latitude
        self._solar_declination = 0.0
        self._sunrise = 0.0
        self._sunset = 0.0
        self._sunlit_increment = 0.0
        self.gauss_integration_params = [
            (0.0469101, 0.1184635),
            (0.2307534, 0.2393144),
            (0.5000000, 0.2844444),
            (0.7692465, 0.2393144),
            (0.9530899, 0.1184635),
        ]
        self.O2_max_coeff = 300
        self.CO2_max_coeff = 300000
        self.O2_conc = 210000
        self.CO2_conc = 245
        self.Vc_max_rate = 200
        self.spec_factor_base = 2600.0
        self.assim_limits_by_temp = self.calculate_assimilation_limits()

    def photosynthesize(
        self,
        grid_par_W_per_sqm,
        _min_temperature,
        _max_temperature,
        last_biomass,
        lai,
        _current_day,
    ):
        print("I am photosynthesizing during the growing season")
        self.update_solar_variables(_current_day)
        total_canopy_assimilated_CO2 = np.zeros_like(last_biomass["leaf_biomass"])
        total_gross_photosynthesis = np.zeros_like(last_biomass["leaf_biomass"])
        for day_increment in self.gauss_integration_params:
            (abscissa, weight) = day_increment
            increment_hour = self._sunrise + abscissa * self._sunlit_increment
            solar_elevation = self.calculate_solar_elevation(increment_hour)
            (
                absorbed_PAR_sunlit,
                absorbed_PAR_shaded,
            ) = self.calculate_absorbed_incremental_PAR(
                increment_hour, solar_elevation, grid_par_W_per_sqm, lai
            )
            sunlit_assimilated_CO2 = self.calculate_leaf_assimilation(
                increment_hour, absorbed_PAR_sunlit, _min_temperature, _max_temperature
            )
            shaded_assimilated_CO2 = self.calculate_leaf_assimilation(
                increment_hour, absorbed_PAR_shaded, _min_temperature, _max_temperature
            )
            sunlit_LAI, shaded_LAI = self.calculate_sunlit_shaded_LAI_proportion(
                increment_hour, lai
            )
            total_gross_photosynthesis += (sunlit_assimilated_CO2 * sunlit_LAI) + (
                shaded_assimilated_CO2 * shaded_LAI
            )
            total_canopy_assimilated_CO2 += (
                total_gross_photosynthesis * weight * 3600 * self._sunlit_increment
            )
        gphot_CH20 = np.zeros_like(last_biomass["leaf_biomass"])
        filter = np.nonzero(total_canopy_assimilated_CO2 > 0)
        gphot_CH20 = total_canopy_assimilated_CO2[filter] * 30 / 1000000
        return gphot_CH20

    # Ignore leaf assimilation for now
    def calculate_leaf_assimilation(
        self, increment_hour, par, _min_temperature, _max_temperature
    ):
        # We are assuming here that air temperature ~canopy temperature
        offset = 1.5
        temp_rise = self._sunrise + offset
        if (increment_hour >= temp_rise) & (increment_hour < self._sunset):
            tau = (
                np.pi
                * (increment_hour * self._sunrise - offset)
                / (self._sunset - self._sunrise)
            )
            hour_temp = _min_temperature + (
                _max_temperature - _min_temperature
            ) * np.sin(tau)
        else:
            if increment_hour < temp_rise:
                increment_hour += 24
            tau = (
                np.pi
                * (increment_hour * self._sunrise - offset)
                / (self._sunset - self._sunrise)
            )
            sunset_temp = _min_temperature + (
                _max_temperature - _min_temperature
            ) * np.sin(tau)
            interp_slope = (_min_temperature - sunset_temp) / (
                temp_rise + 24 - self._sunset
            )
            hour_temp = sunset_temp + interp_slope * (increment_hour * self._sunset)
        max_rubisco = self.get_rubsico_limits(hour_temp)
        max_light_limit = self.calculate_light_limits(par, hour_temp)
        max_sink_limit = self.get_sink_limits(hour_temp)
        print("Rubsico")
        print(max_rubisco)
        print("Light")
        print(max_light_limit)
        print("Sink")
        print(max_sink_limit)
        min_assim = np.min([max_rubisco, max_light_limit, max_sink_limit])
        return min_assim

    def calculate_hourly_direct_light_extinction(self, solar_elevation):
        if solar_elevation < 0.00000001:
            return 0.0
        else:
            return 0.5 / np.sin(solar_elevation)

    def calculate_hourly_diffuse_light_extinction(self, lai):
        return (1.0 + 0.1174 * lai**2) / (1.0 + 0.3732 * lai**2)

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
        self, increment_hour, solar_elevation, grid_par_W_per_sqm, lai
    ):
        P = 0.04  # // reflection coefficient
        A = 0.80  # // scatter coefficient
        S = A**0.5  # // scatter correction
        conv = 4.55  # conversion factor to umol

        hourly_PAR = self.calculate_incremental_PAR(
            increment_hour, grid_par_W_per_sqm
        )  # Need to find this
        hourly_direct_light_extinction_k = (
            self.calculate_hourly_direct_light_extinction(solar_elevation)
        )
        hourly_diffuse_light_extinction_k = (
            self.calculate_hourly_diffuse_light_extinction(lai)
        )
        dI = (1.0 - P) * hourly_PAR
        dIpdr = dI * np.exp(-hourly_direct_light_extinction_k * S * lai)
        dIpdrdr = dI * np.exp(
            -hourly_direct_light_extinction_k * lai
        )  # direct of direct
        dIpdra = (dIpdr - dIpdrdr) / 2  # scatter beams
        dN = hourly_diffuse_light_extinction_k * S * lai
        dI = (1.0 - P) * hourly_PAR
        dIpdf = dI * (1.0 - np.exp(-dN)) / dN  # diffuse

        absorbed_PAR_sunlit = (
            conv * A * (hourly_direct_light_extinction_k * hourly_PAR + dIpdf + dIpdra)
        )
        absorbed_PAR_shaded = conv * A * (dIpdf + dIpdra)
        return (absorbed_PAR_sunlit, absorbed_PAR_shaded)

    def calculate_incremental_PAR(self, increment_hour, grid_par_W_per_sqm):
        dA = np.sin(self._solar_declination) * np.sin(self.latitude)
        dB = np.cos(self._solar_declination) * np.cos(self.latitude)
        dAoB = dA / dB
        dPhi = (np.pi * grid_par_W_per_sqm / 86400) / (
            dA * np.arccos(-dAoB) + dB * np.sqrt(1 - dAoB * dAoB)
        )
        dCoefA = -dB * dPhi
        dCoefB = dA * dPhi
        total_incremental_PAR = dCoefA * np.cos(np.pi * increment_hour / 12) + dCoefB
        return total_incremental_PAR

    def calculate_sunlit_shaded_LAI_proportion(self, increment_hour, lai):
        sunlit_lai = np.zeros_like(lai)
        hourly_direct_light_extinction_k = (
            self.calculate_hourly_direct_light_extinction(increment_hour)
        )
        if hourly_direct_light_extinction_k > 0.0:
            sunlit_lai = (
                1.0 - np.exp(-hourly_direct_light_extinction_k * lai)
            ) / hourly_direct_light_extinction_k
        shaded_lai = lai - sunlit_lai
        return (sunlit_lai, shaded_lai)

    def calculate_light_limits(self, par, hour_temp):
        quantum_yield = 0.06
        absorption_frac = 0.8
        CO2_comp = self.get_CO2_comp(hour_temp)
        dA = absorption_frac * quantum_yield * par * (self.CO2_conc - CO2_comp)
        dB = self.CO2_conc + 2 * CO2_comp
        return dA / dB

    def calculate_assimilation_limits(self):
        # this needs to happen at init then have leaf temp interpolated
        leaf_temp = np.arange(-50, 50, 5)
        O2_coeff = self.O2_max_coeff * 2.1 ** ((leaf_temp - 25) / 10)
        CO2_coeff = self.CO2_max_coeff * 1.2 ** ((leaf_temp - 25) / 10)
        dcoeffm = O2_coeff * (1 + (self.O2_conc / CO2_coeff))
        Vc_exponent = 1 + np.exp(0.128 * (leaf_temp - 40))
        Vc_adj = (self.Vc_max_rate * 2.4 ** ((leaf_temp - 25) / 10)) / Vc_exponent
        CO2_comp = (
            0.5
            * self.O2_conc
            / (self.spec_factor_base * 0.57 ** ((leaf_temp - 25) / 10))
        )
        dA = Vc_adj * (self.CO2_conc - CO2_comp)
        dB = self.CO2_conc + dcoeffm
        rubisco_limits = dA / dB
        sink_limits = 0.5 * Vc_adj
        dtypes = [
            ("leaf_temp", float),
            ("rubisco_limits", float),
            ("sink_limits", float),
            ("CO2_comp", float),
        ]

        limit_map = np.column_stack((leaf_temp, rubisco_limits, sink_limits, CO2_comp))
        limit_map = list(map(tuple, limit_map))
        limit_lookup = np.array(limit_map, dtypes)
        return limit_lookup

    def calculate_light_limited_assimilation(self):
        pass

    def get_rubsico_limits(self, hour_temp):
        limits = np.interp(
            hour_temp,
            self.assim_limits_by_temp["leaf_temp"],
            self.assim_limits_by_temp["rubisco_limits"],
        )
        return limits

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
            self.assim_limits_by_temp["CO2_comp"],
            self.assim_limits_by_temp["sink_limits"],
        )
        return comp_pt


class C3(Photosynthesis):
    def __init__(self, latitude):
        super().__init__(latitude)


class C4(Photosynthesis):
    def __init__(self, latitude):
        super().__init__(latitude)


class Cam(Photosynthesis):
    def __init__(self, latitude):
        super().__init__(latitude)
