"""
Grid-based simulation of coseismic shallow landslides

@author: Suryodoy Ghoshal
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage as _nd
from scipy.ndimage import binary_dilation as _binary_dilation
from scipy.ndimage import generate_binary_structure as _generate_binary_structure
from scipy.ndimage import label as _label
from scipy.special import expit as _expit

from landlab.core.model_component import Component


class _RegionProperties:
    """Subset of ``skimage.measure.regionprops`` used by this component."""

    def __init__(self, labels: np.ndarray, label: int, region_slice: tuple[slice, ...]):
        self.label = label
        self._labels = labels
        self._slice = region_slice
        self._coords = None
        self._inertia_tensor = None
        self._inertia_tensor_eigvals = None

    @property
    def bbox(self) -> tuple[int, int, int, int]:
        return (
            self._slice[0].start,
            self._slice[1].start,
            self._slice[0].stop,
            self._slice[1].stop,
        )

    @property
    def coords(self) -> np.ndarray:
        if self._coords is None:
            local_coords = np.argwhere(self._labels[self._slice] == self.label)
            self._coords = local_coords + np.array(
                [self._slice[0].start, self._slice[1].start]
            )
        return self._coords

    @property
    def area(self) -> float:
        return float(len(self.coords))

    @property
    def centroid(self) -> tuple[float, float]:
        return tuple(self.coords.mean(axis=0))

    @property
    def inertia_tensor(self) -> np.ndarray:
        if self._inertia_tensor is None:
            centered = self.coords.astype(float) - self.coords.mean(axis=0)
            variance_row, variance_col = np.mean(centered**2, axis=0)
            covariance = np.mean(centered[:, 0] * centered[:, 1])
            self._inertia_tensor = np.array(
                [[variance_col, -covariance], [-covariance, variance_row]]
            )
        return self._inertia_tensor

    @property
    def inertia_tensor_eigvals(self) -> np.ndarray:
        if self._inertia_tensor_eigvals is None:
            eigvals = np.linalg.eigvalsh(self.inertia_tensor)
            self._inertia_tensor_eigvals = np.clip(eigvals[::-1], 0.0, None)
        return self._inertia_tensor_eigvals

    @property
    def axis_major_length(self) -> float:
        return float(4.0 * np.sqrt(self.inertia_tensor_eigvals[0]))

    @property
    def axis_minor_length(self) -> float:
        return float(4.0 * np.sqrt(self.inertia_tensor_eigvals[-1]))

    @property
    def orientation(self) -> float:
        a, b, _, c = self.inertia_tensor.flat
        if a == c:
            return np.pi / 4.0 if b < 0 else -np.pi / 4.0
        return float(0.5 * np.arctan2(-2.0 * b, c - a))

    @property
    def eccentricity(self) -> float:
        major, minor = self.inertia_tensor_eigvals
        if major == 0.0:
            return 0.0
        return float(np.sqrt(1.0 - minor / major))


def _regionprops(labels: np.ndarray) -> list[_RegionProperties]:
    """Return properties for each nonzero region in a 2-D label array."""
    labels = np.asarray(labels)
    if labels.ndim != 2:
        raise ValueError("Region properties require a 2-D label array")

    regions = []
    for label, region_slice in enumerate(_nd.find_objects(labels), start=1):
        if region_slice is not None:
            regions.append(_RegionProperties(labels, label, region_slice))
    return regions


class ShallowLandslider(Component):
    r"""
    Predict shallow landslide initiation on a Landlab grid.

    This component computes node-wise stability metrics, identifies and
    sub-groups contiguous unstable regions by aspect, optionally splits groups
    using KDE-informed landslide length-width relationships, and labels the
    resulting candidate landslides.

    The grid must contain topographic elevation, soil depth, and horizontal
    and vertical peak-ground-acceleration fields before initialization.
    """

    _name = "ShallowLandslider"
    _unit_agnostic = True
    _aspect_interval = 20
    _gravity = 9.81
    _handle_small = "merge"

    _info = {
        "topographic__elevation": {
            "intent": "in",
            "mapping": "node",
            "dtype": float,
            "units": "m",
            "optional": False,
            "doc": "Land surface topographic elevation",
        },
        "soil__depth": {
            "intent": "in",
            "mapping": "node",
            "dtype": float,
            "units": "m",
            "optional": False,
            "doc": "Depth of soil or weathered bedrock",
        },
        "earthquake__horizontal_pga": {
            "intent": "in",
            "mapping": "node",
            "dtype": float,
            "units": "g",
            "optional": False,
            "doc": "Horizontal PGA (multiples of g).",
        },
        "earthquake__vertical_pga": {
            "intent": "in",
            "mapping": "node",
            "dtype": float,
            "units": "g",
            "optional": False,
            "doc": "Vertical PGA (multiples of g).",
        },
        "landslide__selected_labels": {
            "intent": "out",
            "mapping": "node",
            "dtype": int,
            "units": "-",
            "optional": False,
            "doc": "Labels of candidate landslides (0 outside a candidate).",
        },
    }

    def __init__(
        self,
        grid,
        cohesion_eff: float = 15000.0,
        angle_int_frict: float = 30.0,
        submerged_soil_proportion: float = 0.5,
    ):
        """
        Initialize the ShallowLandslider component.

        Parameters
        ----------
        grid : landlab.ModelGrid
            Landlab grid on which the component operates.
        cohesion_eff : float
            Effective cohesion (Pa). Scalar; array not required.
        angle_int_frict : float
            Angle of internal friction **in degrees**.
        submerged_soil_proportion : float, optional
            Proportion of submerged soil (0-1); used for suction proxy. Default 0.5.
        """
        super().__init__(grid)
        self._initialize_required_output_fields()
        self.cohesion_eff = float(cohesion_eff)
        self.angle_int_frict = float(np.radians(angle_int_frict))
        self.submerged_soil_proportion = float(submerged_soil_proportion)

        # Internals
        self._fos = None
        self._a_transient = None
        self._a_driving = None
        self._a_diff = None
        self._unstable_mask = None
        self._labels = None
        self._aspect_labels = None
        self._split_labels = None
        self._selected_labels = None

        self._pga_h = self.grid.at_node["earthquake__horizontal_pga"]
        self._pga_v = self.grid.at_node["earthquake__vertical_pga"]

        # Cache aspect
        _asp = self.grid.calc_aspect_at_node(
            elevs="topographic__elevation", unit="degrees", ignore_closed_nodes=True
        )
        self._aspect = np.asarray(_asp, dtype=float).copy()
        self._aspect[self.grid.boundary_nodes] = np.nan

        _slope_rad = self.grid.calc_slope_at_node(
            elevs="topographic__elevation"
        )  # Landlab returns float64
        self._slope_rad64 = np.asarray(_slope_rad, dtype=np.float64)
        self._slope_deg32 = np.degrees(self._slope_rad64).astype(np.float32, copy=False)

    def _initialize_required_output_fields(self):
        """Create required output fields for Landlab component metadata checks."""
        for name, meta in self._info.items():
            if meta["intent"] == "out" and not meta["optional"]:
                at = meta["mapping"]
                if name not in self.grid[at]:
                    self.grid.add_zeros(name, at=at, dtype=meta["dtype"])

    @property
    def results(self) -> dict[str, Any]:
        """
        Return the intermediate and final label arrays.

        Returns
        -------
        dict
            Stability calculations and the connected, aspect-split, KDE-split,
            and final candidate labels.
        """
        return {
            "factor_of_safety": self._fos,
            "a_transient": self._a_transient,
            "a_driving": self._a_driving,
            "a_diff": self._a_diff,
            "unstable_mask": self._unstable_mask,
            "labels": self._labels,
            "aspect_labels": self._aspect_labels,
            "split_labels": self._split_labels,
            "selected_labels": self._selected_labels,
        }

    def run_one_step(
        self, kde_input: dict | None = None, random_seed: int | None = None
    ):
        """
        Identify and label candidate landslides.

        Parameters
        ----------
        kde_input : dict, optional
            Configuration for KDE-based width splitting. If omitted, width
            splitting is skipped.
        random_seed : int, optional
            Seed for reproducible probabilistic inventory selection.
        """
        self._compute_stability()
        self._identify_regions()
        self._filter_by_aspect_and_split(kde_input)
        self._select_groups(random_seed)

    def _select_groups(self, random_seed: int | None):
        """Select a probabilistic inventory from the candidate regions."""
        candidates = (
            self._split_labels
            if self._split_labels is not None
            else self._aspect_labels
        ).reshape(self.grid.shape)
        probabilities, _ = self._generate_landslide_probability(
            labeled_2d=candidates,
            slope_deg_1d=self._slope_deg32,
            critical_accel_1d=self._a_transient,
            h_pga_1d=self._pga_h,
            v_pga_1d=self._pga_v,
            random_seed=random_seed,
        )
        selected, _ = self._probabilistic_group_selection(
            labeled_2d=candidates,
            probability_2d=probabilities,
            random_seed=random_seed,
        )
        self._selected_labels = selected.reshape(self.grid.number_of_nodes)
        self.grid.at_node["landslide__selected_labels"] = self._selected_labels

    # ---------------------------------------------------------------------
    # Pipeline steps
    # ---------------------------------------------------------------------
    def _compute_stability(self):
        """Compute node-wise stability and transient acceleration."""
        self._fos = self._factor_of_safety(
            self.grid,
            self.cohesion_eff,
            self.angle_int_frict,
            submerged_soil_proportion=self.submerged_soil_proportion,
        )
        self._a_transient, self._a_driving, self._a_diff = (
            self._critical_transient_acceleration(
                self.grid,
                self.cohesion_eff,
                self.angle_int_frict,
                submerged_soil_proportion=self.submerged_soil_proportion,
                a_h=self._pga_h * self._gravity,
                a_v=self._pga_v * self._gravity,
            )
        )
        self._unstable_mask = self._a_driving > self._a_transient
        self._unstable_mask[self.grid.boundary_nodes] = False

    def _identify_regions(self):
        """Label contiguous unstable regions."""
        labels, _ = self._calculate_regions(
            self._unstable_mask.reshape(self.grid.shape), connect_val=8
        )
        self._labels = labels.reshape(self.grid.number_of_nodes)

    def _filter_by_aspect_and_split(self, kde_input: dict | None):
        """Split region labels by aspect and KDE-informed width."""
        zones = self._create_zones(interval=self._aspect_interval)
        aspect_subgroups, _, _ = self._split_groups_by_aspect(
            groups=self._labels.reshape(self.grid.shape),
            aspect_array=self._aspect.reshape(self.grid.shape),
            zones=zones,
            handle_small=self._handle_small,
        )
        self._aspect_labels = aspect_subgroups.reshape(self.grid.number_of_nodes)
        self._split_labels = None

        if kde_input is not None:
            split_labels, _ = self._recursive_split_wide_regions(
                labeled_2d=self._aspect_labels.reshape(self.grid.shape),
                aspect_2d=self._aspect.reshape(self.grid.shape),
                slopes_2d=self._slope_deg32.reshape(self.grid.shape),
                kde_results=kde_input.get("kde_data"),
                transform_info=kde_input.get("kde_transform"),
                width_threshold=kde_input.get("width_threshold", 1.5),
                max_iterations=kde_input.get("max_iterations", 10),
                min_region_size=kde_input.get("min_region_size", 10),
                convergence_threshold=kde_input.get("convergence_threshold", 0.75),
            )
            self._split_labels = split_labels.reshape(self.grid.number_of_nodes)

    # STABILITY (inlined from stability.py)
    # ---------------------------------------------------------------------

    def _factor_of_safety(
        self,
        grid,
        cohesion_eff: float,
        angle_int_frict_rad: float,
        submerged_soil_proportion: float = 0.5,
        soil_unit_weight: float = 15e3,
        water_unit_weight: float = 9.8e3,
    ) -> np.ndarray:
        """Compute static factor of safety (FoS) at nodes. (Memory-optimized)"""
        # MEMOPT: avoid full copy
        soil_depth = np.asarray(grid["node"]["soil__depth"], dtype=np.float32)

        # MEMOPT: reuse cached slope in radians (float64 for numerical safety)
        slope = self._slope_rad64  # view

        # Guard small values without creating large temporaries
        soil_depth = soil_depth.copy()  # writing below; keep as float32
        np.maximum(soil_depth, 1e-3, out=soil_depth)  # in-place min clamp
        slope_safe = np.where(slope == 0.0, np.nan, slope)  # small temp; unavoidable

        psi = submerged_soil_proportion * water_unit_weight * soil_depth  # float32
        # Upcasts where needed; results are float64 where slope participates
        fos = (cohesion_eff - psi * np.tan(angle_int_frict_rad)) / (
            soil_unit_weight * soil_depth * np.sin(slope_safe)
        ) + np.tan(angle_int_frict_rad) / np.tan(slope_safe)

        return fos

    def _critical_transient_acceleration(
        self,
        grid,
        cohesion_eff: float,
        angle_int_frict: float,
        submerged_soil_proportion: float,
        a_h: np.ndarray | float = 0.0,
        a_v: np.ndarray | float = 0.0,
        soil_unit_weight: float = 15e3,
        water_unit_weight: float = 9.8e3,
        g: float = 9.81,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute critical transient acceleration (a_c), driving acceleration (a_s),
        and their difference (a_s - a_c).

        Parameters
        ----------
        grid : landlab.ModelGrid
        cohesion_eff : float
        angle_int_frict : float
        submerged_soil_proportion : float
        a_h, a_v : float or array-like
            Horizontal/vertical PGA converted to m/s^2 by multiplying `g` externally.
        soil_unit_weight, water_unit_weight : float
        g : float

        Returns
        -------
        tuple of np.ndarray
            (a_c_transient, a_s_driving, a_difference) in node order.
        """
        soil_depth = np.asarray(grid["node"]["soil__depth"], dtype=np.float64).copy()
        soil_depth[soil_depth == 0] += 0.001  # Avoids division by zero

        slope = self._slope_rad64
        if submerged_soil_proportion >= 0:
            psi = submerged_soil_proportion * water_unit_weight * soil_depth
        else:
            psi = -15e3

        # critical transient acceleration (a_c_transient) in 3D
        a_c_transient = (
            np.tan(angle_int_frict)
            * (g * np.cos(slope) - a_v * np.cos(slope) - a_h * np.sin(slope))
            + ((g * cohesion_eff) - (psi * g * np.tan(angle_int_frict)))
            / (soil_unit_weight * soil_depth)
            - g * np.sin(slope)
        )
        # Driving acceleration downslope
        a_s_t = a_h * np.cos(slope) - a_v * np.sin(slope)

        a_c_transient[grid.boundary_nodes] = 0

        a_difference = a_s_t - a_c_transient

        return a_c_transient, a_s_t, a_difference

    # ---------------------------------------------------------------------
    # REGIONS (inlined from regions.py)
    # ---------------------------------------------------------------------
    def _calculate_regions(
        self,
        binary_grid: np.ndarray,
        connect_val: int = 4,
    ) -> tuple[np.ndarray, int]:
        """Label connected components using four- or eight-cell connectivity."""
        if connect_val == 4:
            struct_arr = _generate_binary_structure(2, 1)
        elif connect_val == 8:
            struct_arr = _generate_binary_structure(2, 2)
        else:
            raise ValueError("Only 4- or 8-connectivity possible")
        return _label(binary_grid > 0, structure=struct_arr)

    def _create_zones(self, interval: int = 20) -> dict[str, tuple[float, float]]:
        """
        Create aspect zones partitioning [0, 360) in `interval` degrees.

        Parameters
        ----------
        interval : int
            Zone width in degrees.

        Returns
        -------
        dict
            Mapping zone-name -> (min_angle, max_angle).
        """
        zones = {}
        start = 0
        while start < 360:
            end = (start + interval) % 360
            zone_name = f"{start:03d}-{end:03d}"
            zones[zone_name] = (start, end)
            start += interval
        return zones

    def _zone_aspects(
        self, aspect_array: np.ndarray, zones: dict | None = None
    ) -> np.ndarray:
        """
        Map aspect values (deg) to integer zone labels.

        Parameters
        ----------
        aspect_array : np.ndarray
            2-D aspect values in degrees (0–360).
        zones : dict, optional
            Zone map from `_create_zones`. If None, 20° zones are used.

        Returns
        -------
        np.ndarray
            2-D integer zone indices, -1 where undefined.
        """
        if zones is None:
            zones = self._create_zones(20)
        aspect_array = np.array(aspect_array, dtype=float) % 360
        zone_array = np.full(aspect_array.shape, -1, dtype=int)
        for zone_idx, (_, (min_aspect, max_aspect)) in enumerate(zones.items()):
            if min_aspect > max_aspect:
                mask = (aspect_array >= min_aspect) | (aspect_array <= max_aspect)
            else:
                mask = (aspect_array >= min_aspect) & (aspect_array < max_aspect)
            zone_array[mask] = zone_idx
        return zone_array

    def _split_groups_by_aspect(
        self,
        groups: np.ndarray,
        aspect_array: np.ndarray,
        zones: dict | None = None,
        min_size: int = 2,
        handle_small: str = "merge",
    ) -> tuple[np.ndarray, np.ndarray, dict]:
        """
        Split connected components by aspect zones; optionally merge/remove small parts.

        Parameters
        ----------
        groups : np.ndarray
            2-D integer labels of unstable regions.
        aspect_array : np.ndarray
            2-D aspect values in degrees.
        zones : dict, optional
            Zone definition from `_create_zones`.
        min_size : int
            Minimum pixels for a component to be kept.
        handle_small : {"merge", "remove", "keep"}
            Strategy for small components.
        Returns
        -------
        tuple
            (new_groups, zone_labels, group_info)
        """
        if zones is None:
            zones = self._create_zones(20)
        zone_labels = self._zone_aspects(aspect_array, zones)
        new_groups = np.zeros_like(groups)
        group_info = {}
        next_label = 1
        zone_names = list(zones.keys())
        small_regions = []
        # Use regionprops to get bounding boxes — avoids full-grid scans per group
        region_props_list = _regionprops(groups)
        for region in region_props_list:
            group_id = region.label
            r0, c0, r1, c1 = region.bbox
            # Pad by 1 so dilation later stays in bounds
            r0p = max(0, r0 - 1)
            r1p = min(groups.shape[0], r1 + 1)
            c0p = max(0, c0 - 1)
            c1p = min(groups.shape[1], c1 + 1)

            sub_groups = groups[r0p:r1p, c0p:c1p]
            sub_zones = zone_labels[r0p:r1p, c0p:c1p]
            group_submask = sub_groups == group_id

            for zone_id in np.unique(sub_zones[group_submask]):
                combined_mask = group_submask & (sub_zones == zone_id)
                sub_labels, num_features = _label(combined_mask)
                for label_name in range(1, num_features + 1):
                    component_mask = sub_labels == label_name
                    component_size = int(np.sum(component_mask))
                    local_rows, local_cols = np.where(component_mask)
                    abs_rows = local_rows + r0p
                    abs_cols = local_cols + c0p
                    if component_size < min_size:
                        small_regions.append(
                            {
                                "rows": abs_rows,
                                "cols": abs_cols,
                                "group_id": group_id,
                                "zone_id": zone_id,
                                "size": component_size,
                                "centroid": (
                                    float(np.mean(abs_rows)),
                                    float(np.mean(abs_cols)),
                                ),
                            }
                        )
                    else:
                        new_groups[abs_rows, abs_cols] = next_label
                        group_info[next_label] = (group_id, zone_names[zone_id])
                        next_label += 1

        if small_regions and handle_small == "merge":
            small_regions.sort(key=lambda x: x["size"])
            for region in small_regions:
                rmin = int(region["rows"].min())
                rmax = int(region["rows"].max()) + 1
                cmin = int(region["cols"].min())
                cmax = int(region["cols"].max()) + 1
                # Bbox padded by 1 for dilation
                r0 = max(0, rmin - 1)
                r1 = min(new_groups.shape[0], rmax + 1)
                c0 = max(0, cmin - 1)
                c1 = min(new_groups.shape[1], cmax + 1)

                local = np.zeros((r1 - r0, c1 - c0), dtype=bool)
                local[region["rows"] - r0, region["cols"] - c0] = True
                dilated_local = _binary_dilation(local)
                neighbor_local = dilated_local & ~local

                sub_new_groups = new_groups[r0:r1, c0:c1]
                neighbor_labels = np.unique(sub_new_groups[neighbor_local])
                neighbor_labels = neighbor_labels[neighbor_labels > 0]
                if len(neighbor_labels) > 0:
                    neighbor_counts = [
                        (nl, np.sum(sub_new_groups[neighbor_local] == nl))
                        for nl in neighbor_labels
                    ]
                    best_neighbor = max(neighbor_counts, key=lambda x: x[1])[0]
                    new_groups[region["rows"], region["cols"]] = best_neighbor
                else:
                    new_groups[region["rows"], region["cols"]] = next_label
                    group_info[next_label] = (
                        region["group_id"],
                        zone_names[region["zone_id"]],
                    )
                    next_label += 1
        elif small_regions and handle_small == "remove":
            pass
        return new_groups, zone_labels, group_info

    def _regions__handle_small_regions(
        self,
        labeled_array: np.ndarray,
        min_size: int,
        method: str = "merge",
        grid=None,
    ) -> np.ndarray:
        """
        Handle labeled regions smaller than a minimum size by merging or removing.

        Parameters
        ----------
        labeled_array : np.ndarray
        min_size : int
        method : {"merge", "remove"}
        grid : landlab.ModelGrid, optional
        Returns
        -------
        np.ndarray
            Modified labeled array.
        """
        if method not in ["merge", "remove"]:
            return labeled_array
        modified = labeled_array.copy()
        props = _regionprops(labeled_array)
        smalls = [
            {
                "label": r.label,
                "mask": labeled_array == r.label,
                "centroid": r.centroid,
                "size": r.area,
            }
            for r in props
            if r.area < min_size
        ]
        if not smalls:
            return labeled_array
        if method == "remove":
            for region in smalls:
                modified[region["mask"]] = 0
            return modified
        smalls.sort(key=lambda x: x["size"])
        for region in smalls:
            mask = region["mask"]
            dilated = _binary_dilation(mask)
            neighbor_mask = dilated & ~mask
            neighbor_labels = np.unique(modified[neighbor_mask])
            neighbor_labels = neighbor_labels[
                (neighbor_labels > 0) & (neighbor_labels != region["label"])
            ]
            if len(neighbor_labels) > 0:
                neighbor_counts = [
                    (nl, np.sum(modified[neighbor_mask] == nl))
                    for nl in neighbor_labels
                ]
                best_neighbor = max(neighbor_counts, key=lambda x: x[1])[0]
                modified[mask] = best_neighbor
        return modified

    def _generate_landslide_probability(
        self,
        labeled_2d: np.ndarray,
        slope_deg_1d: np.ndarray | None,
        critical_accel_1d: np.ndarray | None,
        h_pga_1d: np.ndarray,
        v_pga_1d: np.ndarray,
        random_seed: int | None,
        normalize: bool = True,
    ) -> tuple[np.ndarray, dict]:
        """
        Build per-group failure probabilities (optionally normalized).

        Parameters
        ----------
        labeled_2d : np.ndarray
            2-D labels for candidate regions.
        slope_deg_1d : np.ndarray or None
            Node-wise slope angles (degrees).
        critical_accel_1d : np.ndarray or None
            Critical acceleration (m/s^2). If None, default 0.2 used.
        h_pga_1d, v_pga_1d : np.ndarray
            Horizontal/vertical PGA (multiples of g).
        random_seed : int or None
            Seed for reproducibility.
        normalize : bool
            Min-max normalize group probabilities.

        Returns
        -------
        probability_2d : np.ndarray
            2-D per-node probabilities by group.
        metadata : dict
            Details of normalization and group-level metrics.
        """
        h_grid = np.asarray(h_pga_1d).reshape(self.grid.shape)
        v_grid = np.asarray(v_pga_1d).reshape(self.grid.shape)
        if random_seed is not None:
            np.random.seed(random_seed)
        if critical_accel_1d is None:
            crit_grid = np.full(self.grid.shape, 0.2, dtype=np.float32)
        else:
            crit_grid = np.asarray(critical_accel_1d).reshape(self.grid.shape)
        slope_grid = None
        if slope_deg_1d is not None:
            slope_grid = np.asarray(slope_deg_1d).reshape(self.grid.shape)

        unique_labels = np.unique(labeled_2d)[1:]
        if len(unique_labels) == 0:
            return np.zeros_like(labeled_2d, dtype=np.float32), {}

        # ------------------------------------------------------------------
        # Compute all per-group means in one vectorised call each.
        # scipy.ndimage.mean never materialises per-group boolean masks.
        # ------------------------------------------------------------------
        mean_h = np.asarray(_nd.mean(h_grid, labels=labeled_2d, index=unique_labels))
        mean_v = np.asarray(_nd.mean(v_grid, labels=labeled_2d, index=unique_labels))
        mean_crit = np.asarray(
            _nd.mean(crit_grid, labels=labeled_2d, index=unique_labels)
        )
        if slope_grid is not None:
            mean_slope = np.asarray(
                _nd.mean(slope_grid, labels=labeled_2d, index=unique_labels)
            )
        else:
            mean_slope = None

        # Scalar results only — no masks stored anywhere.
        group_probs = {}
        raw_probs = np.zeros(len(unique_labels), dtype=np.float64)

        for i, lab in enumerate(unique_labels):
            local_crit = float(mean_crit[i])
            mh = float(mean_h[i])
            mv = float(mean_v[i])
            resultant = float(np.sqrt(mh**2 + mv**2))
            pga_ratio = resultant / local_crit if local_crit > 0 else np.inf
            base_prob = self._sel___calculate_acceleration_probability(
                pga_ratio, local_crit
            )
            if mean_slope is not None:
                base_prob *= self._sel___calculate_slope_stability_factor(
                    float(mean_slope[i])
                )
            stochastic = float(np.random.lognormal(mean=0, sigma=0.2))
            group_prob = float(np.clip(base_prob * stochastic, 0, 1))
            raw_probs[i] = group_prob
            group_probs[int(lab)] = {
                "probability": group_prob,
                "critical_acceleration": local_crit,
                "resultant_pga": resultant,
                "pga_ratio": pga_ratio,
                "base_probability": base_prob,
                # no "mask" key — that was the memory killer
            }

        # ------------------------------------------------------------------
        # Normalise if requested, then paint prob_2d with a lookup array.
        # lookup[label] = probability; prob_2d = lookup[labeled_2d] is a
        # single vectorised index — no per-group masks, no loop.
        # ------------------------------------------------------------------
        max_label = int(labeled_2d.max())
        prob_lookup = np.zeros(max_label + 1, dtype=np.float32)

        if normalize and len(raw_probs) > 1:
            min_p, max_p = float(raw_probs.min()), float(raw_probs.max())
            if max_p > min_p:
                norm_vals = (raw_probs - min_p) / (max_p - min_p)
                norm_meta = {
                    "performed": True,
                    "min_raw_prob": min_p,
                    "max_raw_prob": max_p,
                }
            else:
                norm_vals = raw_probs.copy()
                norm_meta = {
                    "performed": False,
                    "reason": "All groups have the same probability",
                    "value": min_p,
                }
            for i, lab in enumerate(unique_labels):
                nv = float(norm_vals[i])
                prob_lookup[int(lab)] = nv
                group_probs[int(lab)]["normalized_probability"] = nv
        else:
            norm_meta = {"performed": False, "reason": "Only one group present"}
            for i, lab in enumerate(unique_labels):
                prob_lookup[int(lab)] = raw_probs[i]

        prob_2d = prob_lookup[labeled_2d].astype(np.float32)
        del prob_lookup  # free immediately

        meta = self._sel__create_metadata(
            group_probs, prob_2d, norm_meta, normalized=normalize
        )
        return prob_2d, meta

    def _sel__create_metadata(
        self,
        group_probs: dict,
        prob_arr: np.ndarray,
        norm_meta: dict,
        normalized: bool = False,
    ) -> dict:
        """
        Collect metadata for selection probabilities and normalization.
        """
        meta = {"group_details": [], "normalization": norm_meta}
        for lab, info in group_probs.items():
            entry = {
                "label": lab,
                "critical_acceleration": info["critical_acceleration"],
                "resultant_pga": info["resultant_pga"],
                "pga_ratio": info["pga_ratio"],
                "base_probability": info["base_probability"],
            }
            if normalized and "normalized_probability" in info:
                entry["raw_probability"] = info["probability"]
                entry["final_probability"] = info["normalized_probability"]
            else:
                entry["final_probability"] = info["probability"]
            meta["group_details"].append(entry)
        nz = prob_arr[prob_arr > 0]
        if len(nz) > 0:
            meta["overall_proportion"] = float(np.mean(nz))
            meta["max_proportion"] = float(np.max(nz))
            meta["min_proportion"] = float(np.min(nz))
        else:
            meta["overall_proportion"] = 0.0
            meta["max_proportion"] = 0.0
            meta["min_proportion"] = 0.0
        return meta

    def _sel___calculate_acceleration_probability(
        self, pga_ratio: float, critical_acceleration: float
    ) -> float:
        """
        Combine vulnerability term and PGA ratio via a sigmoid-like mapping.
        """
        if critical_acceleration <= 0:
            return 1.0
        epsilon = 1e-3
        exponent = -5 * (1.0 / max(critical_acceleration, epsilon))
        exponent = np.clip(exponent, -700, 700)
        base_probability = 1 - np.exp(exponent)
        pga_effect = _expit(5 * (pga_ratio - 1))
        return float(np.clip(base_probability * pga_effect, 0, 1))

    def _sel___calculate_slope_stability_factor(self, slope_deg: float) -> float:
        """
        Empirical slope-stability weighting based on slope angle (deg).
        """
        if slope_deg <= 30:
            return 1 + 0.02 * slope_deg
        elif 30 < slope_deg <= 45:
            return 1 + 0.6 * ((slope_deg - 30) / 15)
        else:
            return 2.0

    def _probabilistic_group_selection(
        self,
        labeled_2d: np.ndarray,
        probability_2d: np.ndarray,
        proportion_method: str = "empirical",
        custom_proportion: float | None = None,
        random_seed: int | None = 5000,
        reproducible: bool = True,
    ) -> tuple[np.ndarray, dict]:
        """
        Select groups probabilistically using mean probability and proportion rule.
        """
        unique_labels = np.unique(labeled_2d)
        unique_labels = unique_labels[unique_labels != 0]
        num_groups = len(unique_labels)
        if reproducible and random_seed is not None:
            np.random.seed(random_seed)
        group_means = _nd.mean(probability_2d, labels=labeled_2d, index=unique_labels)
        group_means = np.array(group_means, dtype=float)
        if custom_proportion is not None:
            proportion = custom_proportion
            method = "user_defined"
        else:
            proportion = self._sel___calculate_landslide_proportion(group_means)
            method = "conservative"
        if np.sum(group_means) == 0:
            return np.zeros_like(labeled_2d), {
                "method_used": method,
                "proportion_calculated": 0.0,
                "num_groups_total": num_groups,
                "num_groups_selected": 0,
                "selected_labels": [],
                "group_probabilities": {},
                "selection_probabilities": {},
            }
        normalized = group_means / np.sum(group_means)
        num_to_select = max(1, int(np.ceil(num_groups * proportion)))
        num_to_select = min(num_to_select, num_groups)
        selected = np.random.choice(
            unique_labels, num_to_select, replace=False, p=normalized
        )
        max_label = int(labeled_2d.max())
        sel_lookup = np.zeros(max_label + 1, dtype=bool)
        sel_lookup[selected] = True
        selected_groups = np.where(sel_lookup[labeled_2d], labeled_2d, 0)
        del sel_lookup
        metadata = {
            "method_used": method,
            "proportion_calculated": float(proportion),
            "num_groups_total": num_groups,
            "num_groups_selected": num_to_select,
            "selected_labels": selected.tolist(),
            "group_probabilities": dict(
                zip(unique_labels.tolist(), group_means.tolist())
            ),
            "selection_probabilities": dict(
                zip(unique_labels.tolist(), normalized.tolist())
            ),
        }
        return selected_groups, metadata

    def _sel___calculate_landslide_proportion(self, group_probs: np.ndarray) -> float:
        """Compute a conservative proportion of landslide groups to select."""
        valid = group_probs[group_probs > 0]
        if len(valid) == 0:
            return 0.0
        q90 = float(np.percentile(valid, 90))
        q95 = float(np.percentile(valid, 95))
        mean_probability = float(np.mean(valid))
        threshold = 0.7 * q90 + 0.3 * q95
        proportion = np.sum(valid >= threshold) / len(valid)
        if mean_probability > 0.6:
            proportion *= 1.5
        elif mean_probability > 0.4:
            proportion *= 1.2
        return float(np.clip(proportion, 0.05, 0.50))

    # ---------------------------------------------------------------------
    # SPLIT (KDE-based) (inlined from split.py)
    # ---------------------------------------------------------------------
    def _recursive_split_wide_regions(
        self,
        labeled_2d: np.ndarray,
        aspect_2d: np.ndarray,
        slopes_2d: np.ndarray,
        kde_results,
        transform_info,
        width_threshold: float = 1.5,
        max_iterations: int = 10,
        min_region_size: int = 10,
        convergence_threshold: float = 0.95,
    ) -> tuple[np.ndarray, list]:
        """
        Recursively split regions whose actual width exceeds KDE-expected width.

        Parameters
        ----------
        labeled_2d : np.ndarray
        aspect_2d : np.ndarray
        slopes_2d : np.ndarray
        kde_results : dict
            Must contain key 'overall' with a gaussian_kde instance.
        transform_info : dict
            Contains transformation flags (e.g., log_x/log_y).
        width_threshold : float
            Split when actual/expected width > threshold.
        max_iterations : int
            Maximum recursive iterations.
        min_region_size : int
            Minimum pixels to consider splitting.
        convergence_threshold : float
            Stop if fraction of conforming regions exceeds threshold.
        Returns
        -------
        (final_labels, all_split_info)
        """
        elevation_grid = self.grid.at_node["topographic__elevation"]
        current_labels = labeled_2d.copy()
        all_split_info = []

        for iteration in range(max_iterations):
            unique_labels = np.unique(current_labels)
            unique_labels = unique_labels[unique_labels != 0]

            props = self._split__calculate_region_dimensions(
                current_labels,
                elevation_grid,
                aspect_2d,
                slopes_2d,
                self.grid,
                unique_labels,
            )
            region_df = pd.DataFrame(
                {
                    "label": props["label"],
                    "length_m": props["slope_direction_length_new"],
                    "width_m": props["perpendicular_width_new"],
                    "area": props["area"],
                    "direction_method": props["direction_method"],
                }
            )
            region_df = region_df[region_df["area"] >= min_region_size]

            if len(region_df) == 0:
                break

            new_labels, split_info = self._split__split_wide_regions_single_iteration(
                self.grid,
                current_labels,
                region_df,
                kde_results,
                transform_info,
                width_threshold,
            )
            num_splits = len(split_info)
            total_regions = len(region_df)
            conforming = total_regions - num_splits
            conformance_rate = (
                (conforming / total_regions) if total_regions > 0 else 1.0
            )

            for split in split_info:
                split["iteration"] = iteration + 1
            all_split_info.extend(split_info)

            if num_splits == 0:
                break
            elif conformance_rate >= convergence_threshold:
                break

            current_labels = new_labels

        return current_labels, all_split_info

    def _split__calculate_region_dimensions(
        self,
        labeled_2d: np.ndarray,
        elevation_1d: np.ndarray,
        aspect_2d: np.ndarray,
        slopes_2d: np.ndarray,
        grid,
        unique_labels=None,
    ) -> dict:
        """
        Compute length/width for each region using elevation or aspect direction.
        """
        elevation_grid = elevation_1d.reshape(grid.shape)
        slopes_grid = slopes_2d.reshape(grid.shape)
        aspect_grid = aspect_2d.reshape(grid.shape)
        if unique_labels is None:
            unique_labels = np.unique(labeled_2d)
            unique_labels = unique_labels[unique_labels != 0]
        regions = _regionprops(labeled_2d)
        labels, areas, min_elevations, max_elevations = [], [], [], []
        region_coords_map = {}  # cache coords to avoid re-scanning full grid below
        for region in regions:
            labels.append(region.label)
            areas.append(region.area)
            rc = region.coords  # shape (n_pixels, 2)
            region_coords_map[region.label] = rc
            region_elev = elevation_grid[rc[:, 0], rc[:, 1]]
            min_elevations.append(np.min(region_elev))
            max_elevations.append(np.max(region_elev))
        labels = np.array(labels)
        areas = np.array(areas)
        min_elevations = np.array(min_elevations)
        max_elevations = np.array(max_elevations)
        result_props = {
            "label": labels,
            "area": areas,
            "min_elevation": min_elevations,
            "max_elevation": max_elevations,
            "local_relief": max_elevations - min_elevations,
            "slope_direction_length_new": np.zeros(len(labels)),
            "perpendicular_width_new": np.zeros(len(labels)),
            "direction_method": [""] * len(labels),
        }
        relief_threshold = 2.0
        for label_num in unique_labels:
            if label_num == 0:
                continue
            prop_idx = np.where(result_props["label"] == label_num)[0]
            if len(prop_idx) == 0:
                continue
            prop_idx = prop_idx[0]
            elevation_relief = result_props["local_relief"][prop_idx]
            rc = region_coords_map.get(label_num)
            if rc is None or len(rc) == 0:
                continue
            row_coords = rc[:, 0]
            col_coords = rc[:, 1]
            x_coords = col_coords * grid.dx
            y_coords = row_coords * grid.dy
            coords = np.column_stack([x_coords, y_coords])
            gradient_direction = None
            method_used = None
            if elevation_relief > relief_threshold:
                region_elev = elevation_grid[row_coords, col_coords]
                max_elev = result_props["max_elevation"][prop_idx]
                min_elev = result_props["min_elevation"][prop_idx]
                tol = 0.1 * elevation_relief
                high_idx = region_elev >= (max_elev - tol)
                low_idx = region_elev <= (min_elev + tol)
                high_c = coords[high_idx]
                low_c = coords[low_idx]
                if len(high_c) > 0 and len(low_c) > 0:
                    downslope_vec = np.mean(low_c, axis=0) - np.mean(high_c, axis=0)
                    grad_len = np.linalg.norm(downslope_vec)
                    if grad_len > 0:
                        gradient_direction = downslope_vec / grad_len
                        method_used = "elevation_gradient"
            if gradient_direction is None:
                region_aspects = aspect_grid[row_coords, col_coords]
                region_slopes = slopes_grid[row_coords, col_coords]
                slope_threshold = max(np.percentile(region_slopes, 25), 1.0)
                valid_mask = region_slopes > slope_threshold
                if np.sum(valid_mask) > 0:
                    valid_aspects = region_aspects[valid_mask]
                    valid_slopes = region_slopes[valid_mask]
                    downslope_aspects_rad = (valid_aspects + 180) * np.pi / 180
                    cos_a = np.cos(downslope_aspects_rad) * valid_slopes
                    sin_a = np.sin(downslope_aspects_rad) * valid_slopes
                    mean_cos = np.sum(cos_a) / np.sum(valid_slopes)
                    mean_sin = np.sum(sin_a) / np.sum(valid_slopes)
                    mean_dir_rad = np.arctan2(mean_sin, mean_cos)
                    cart_angle = np.pi / 2 - mean_dir_rad
                    gradient_direction = np.array(
                        [np.cos(cart_angle), np.sin(cart_angle)]
                    )
                    method_used = "aspect_weighted"
                else:
                    valid_aspects = region_aspects[~np.isnan(region_aspects)]
                    if len(valid_aspects) > 0:
                        downslope_aspects_rad = (valid_aspects + 180) * np.pi / 180
                        cos_a = np.cos(downslope_aspects_rad)
                        sin_a = np.sin(downslope_aspects_rad)
                        mean_dir_rad = np.arctan2(np.mean(sin_a), np.mean(cos_a))
                        cart_angle = np.pi / 2 - mean_dir_rad
                        gradient_direction = np.array(
                            [np.cos(cart_angle), np.sin(cart_angle)]
                        )
                        method_used = "aspect_simple"
                    else:
                        r = _regionprops(labeled_2d == label_num)[0]
                        orientation = r.orientation
                        gradient_direction = np.array(
                            [np.cos(orientation), np.sin(orientation)]
                        )
                        method_used = "region_orientation"
            perp_dir = np.array([-gradient_direction[1], gradient_direction[0]])
            centroid = np.mean(coords, axis=0)
            centered = coords - centroid
            grad_proj = np.dot(centered, gradient_direction)
            perp_proj = np.dot(centered, perp_dir)
            result_props["slope_direction_length_new"][prop_idx] = np.max(
                grad_proj
            ) - np.min(grad_proj)
            result_props["perpendicular_width_new"][prop_idx] = np.max(
                perp_proj
            ) - np.min(perp_proj)
            result_props["direction_method"][prop_idx] = method_used
        return result_props

    def _split__split_wide_regions_single_iteration(
        self,
        grid,
        labeled_2d: np.ndarray,
        region_df: pd.DataFrame,
        kde_results,
        transform_info,
        width_threshold: float = 1.5,
        label_col: str = "label",
        length_col: str = "length_m",
        width_col: str = "width_m",
    ) -> tuple[np.ndarray, list]:
        """
        Perform one iteration of splitting using KDE-expected widths and main axis direction.
        """
        new_labels = labeled_2d.copy()
        next_label = np.max(labeled_2d) + 1 if labeled_2d.size > 0 else 1
        split_info = []
        kde = kde_results["overall"]
        log_x = transform_info.get("log_x", False)
        log_y = transform_info.get("log_y", False)
        regions_to_split = []
        for _, row in region_df.iterrows():
            label_id = row[label_col]
            if label_id == 0:
                continue
            length = row[length_col]
            actual_width = row[width_col]
            if length <= 0 or actual_width <= 0:
                continue
            length_t = np.log(length) if log_x else length
            num_samples = 200
            max_attempts = 500
            candidates = kde.resample(max_attempts)  # (2, 500), one call
            tol = 0.05 * (1 + abs(length_t))
            keep = np.abs(candidates[0] - length_t) < tol
            samples = candidates[1, keep][:num_samples]
            if len(samples) < 10:
                samples = kde.resample(num_samples)[1, :]
            expected_widths = np.exp(samples) if log_y else np.asarray(samples)
            expected_width = float(np.median(expected_widths))
            width_ratio = actual_width / expected_width
            if width_ratio > width_threshold:
                main_direction = self._split__extract_direction_from_region_data(
                    row, labeled_2d, grid
                )
                regions_to_split.append(
                    {
                        "label": label_id,
                        "actual_width": float(actual_width),
                        "expected_width": float(expected_width),
                        "ratio": float(width_ratio),
                        "main_direction": main_direction,
                        "direction_method": row.get("direction_method", "unknown"),
                    }
                )
        for region in regions_to_split:
            label_id = region["label"]
            mask = labeled_2d == label_id
            rrows, rcols = np.where(mask)
            if len(rrows) == 0:
                continue
            x_coords = rcols * grid.dx
            y_coords = rrows * grid.dy
            coords = np.column_stack([x_coords, y_coords])
            if region["main_direction"] is not None:
                split_dir = np.array(
                    [-region["main_direction"][1], region["main_direction"][0]]
                )
            else:
                split_dir = self._split__get_pca_split_direction(coords)
            centroid = np.mean(coords, axis=0)
            centered = coords - centroid
            projections = np.dot(centered, split_dir)
            split_mask = projections >= 0
            part1_rows = rrows[split_mask]
            part1_cols = rcols[split_mask]
            part2_rows = rrows[~split_mask]
            part2_cols = rcols[~split_mask]
            if len(part1_rows) == 0 or len(part2_rows) == 0:
                continue
            part1_mask = np.zeros_like(mask, dtype=bool)
            part2_mask = np.zeros_like(mask, dtype=bool)
            part1_mask[part1_rows, part1_cols] = True
            part2_mask[part2_rows, part2_cols] = True
            new_labels[part1_mask] = next_label
            part1_label = next_label
            next_label += 1
            new_labels[part2_mask] = next_label
            part2_label = next_label
            next_label += 1
            region["part1_label"] = int(part1_label)
            region["part2_label"] = int(part2_label)
            region["original_size"] = int(len(rrows))
            region["part1_size"] = int(len(part1_rows))
            region["part2_size"] = int(len(part2_rows))
            region["split_direction"] = (
                split_dir.tolist() if region["main_direction"] is not None else None
            )
            split_info.append(region)
        return new_labels, split_info

    def _split__extract_direction_from_region_data(
        self, region_row: pd.Series, labeled_2d: np.ndarray, grid
    ) -> np.ndarray | None:
        """
        Reconstruct main axis direction using PCA of region coordinates.
        """
        label_id = region_row["label"]
        mask = labeled_2d == label_id
        rrows, rcols = np.where(mask)
        if len(rrows) == 0:
            return None
        x_coords = rcols * grid.dx
        y_coords = rrows * grid.dy
        coords = np.column_stack([x_coords, y_coords])
        return self._split__get_region_major_axis(coords)

    def _split__get_region_major_axis(self, coords: np.ndarray) -> np.ndarray:
        """
        Compute region major-axis unit vector via PCA.
        """
        if len(coords) < 2:
            return np.array([1, 0])
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid
        cov = np.cov(centered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        major_idx = np.argmax(eigvals)
        major_axis = eigvecs[:, major_idx]
        if major_axis[0] < 0:
            major_axis = -major_axis
        return major_axis

    def _split__get_pca_split_direction(self, coords: np.ndarray) -> np.ndarray:
        """
        Return direction perpendicular to the major axis (PCA fallback).
        """
        major_axis = self._split__get_region_major_axis(coords)
        return np.array([-major_axis[1], major_axis[0]])
