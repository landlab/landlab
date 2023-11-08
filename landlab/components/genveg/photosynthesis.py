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

    def photosynthesize(self, grid_par_W_per_sqm, last_biomass, lai, _current_day):
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
                absorbed_PAR_sunlit
            )
            shaded_assimilated_CO2 = self.calculate_leaf_assimilation(
                absorbed_PAR_shaded
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
    def calculate_leaf_assimilation(self, par):
        return par

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


class C3(Photosynthesis):
    def __init__(self, latitude):
        super().__init__(latitude)


class C4(Photosynthesis):
    def __init__(self, latitude):
        super().__init__(latitude)


class Cam(Photosynthesis):
    def __init__(self, latitude):
        super().__init__(latitude)
