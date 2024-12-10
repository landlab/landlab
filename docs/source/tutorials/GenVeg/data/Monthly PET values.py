try:
            self._PET = self._grid["cell"][
                "surface__potential_evapotranspiration_rate"
            ][:].copy()
        except KeyError:
            msg = (
                "Potential evapotranspiration and/or soil moisture fields not found."
                "GenVeg will use monthly PET values for Atlantic City NJ."
            )
            print(msg)
            # Data from NE Regional Climate Center for Atlantic City NJ
            # https://www.nrcc.cornell.edu/wxstation/pet/pet.html
            _days = np.array([14, 45, 73, 104, 134, 165, 195, 226, 257, 287, 318, 348])
            _days_in_mo = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
            _PET_in = np.array(
                [0.52, 0.76, 1.48, 2.49, 3.75, 4.36, 4.84, 4.18, 2.79, 1.72, 0.87, 0.53]
            )
            _PET_mm_day = _PET_in / _days_in_mo * 25.4
            self._PET = np.interp(
                self._calc_current_jday(), _days, _PET_mm_day, 365
            ) * np.ones_like(self._par)