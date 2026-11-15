"""Scenario-driven Bekker-Wong calculations; this module does not invent soil inputs."""
from __future__ import annotations

import math

import numpy as np


def particle_size_summary(csv_frame):
    """Validate a mass-fraction particle-size table and return measured PSD quantiles."""
    required = {"particle_size_mm", "mass_fraction"}
    missing = required - set(csv_frame.columns)
    if missing:
        raise ValueError("PSD CSV requires columns: particle_size_mm, mass_fraction")
    sizes = np.asarray(csv_frame["particle_size_mm"], dtype=float)
    fractions = np.asarray(csv_frame["mass_fraction"], dtype=float)
    valid = np.isfinite(sizes) & np.isfinite(fractions) & (sizes > 0) & (fractions >= 0)
    sizes, fractions = sizes[valid], fractions[valid]
    if not sizes.size or fractions.sum() <= 0:
        raise ValueError("PSD requires positive sizes and positive total mass fraction")
    order = np.argsort(sizes)
    sizes, fractions = sizes[order], fractions[order] / fractions.sum()
    cumulative = np.cumsum(fractions)
    quantiles = {}
    for percent in (10, 50, 90):
        index = min(int(np.searchsorted(cumulative, percent / 100, side="left")), len(sizes) - 1)
        quantiles[f"D{percent}_mm"] = float(sizes[index])
    return {"particle_count": int(len(sizes)), "total_fraction": float(fractions.sum()),
            "quantiles": quantiles, "sizes_mm": sizes.tolist(),
            "mass_fractions": fractions.tolist()}


def bekker_wong_scenario(kc, kphi, exponent, wheel_width_m, contact_length_m,
                          wheel_load_n, wheel_radius_m, cohesion_pa,
                          internal_friction_deg, shear_modulus_m, max_shear_displacement_m,
                          curve_points=100):
    """Compute a parameterized pressure-sinkage curve and simplified traction/risk outputs.

    Bekker: p = (kc / b + kphi) z^n.
    Wong-Janosi-Hanamoto shear: tau = (c + p tan(phi)) (1 - exp(-j/K)).
    The risk index is a transparent screening heuristic, not a validated vehicle model.
    """
    positive = {"kc": kc, "kphi": kphi, "exponent": exponent,
                "wheel_width_m": wheel_width_m, "contact_length_m": contact_length_m,
                "wheel_load_n": wheel_load_n, "wheel_radius_m": wheel_radius_m,
                "shear_modulus_m": shear_modulus_m,
                "max_shear_displacement_m": max_shear_displacement_m}
    invalid = [name for name, value in positive.items() if not np.isfinite(value) or value <= 0]
    if invalid:
        raise ValueError("Positive measured/scenario values required for: " + ", ".join(invalid))
    if not np.isfinite(cohesion_pa) or not np.isfinite(internal_friction_deg) or cohesion_pa < 0 or not 0 <= internal_friction_deg < 90:
        raise ValueError("cohesion must be nonnegative and friction angle between 0 and 90 degrees")
    if curve_points < 2:
        raise ValueError("curve_points must be at least 2")

    pressure_coefficient = kc / wheel_width_m + kphi
    contact_area = wheel_width_m * contact_length_m
    normal_pressure = wheel_load_n / contact_area
    sinkage_m = (normal_pressure / pressure_coefficient) ** (1 / exponent)

    pressures = np.linspace(0, normal_pressure * 2, curve_points)
    sinkages = np.where(pressures > 0,
                        (pressures / pressure_coefficient) ** (1 / exponent), 0.0)
    shear_limit = cohesion_pa + normal_pressure * math.tan(math.radians(internal_friction_deg))
    shear_factor = 1 - math.exp(-max_shear_displacement_m / shear_modulus_m)
    gross_traction_n = shear_limit * shear_factor * contact_area
    compaction_resistance_n = (wheel_width_m * pressure_coefficient *
                               sinkage_m ** (exponent + 1) / (exponent + 1))
    drawbar_pull_n = gross_traction_n - compaction_resistance_n
    efficiency = float(np.clip(drawbar_pull_n / gross_traction_n, 0, 1)) if gross_traction_n > 0 else 0.0
    sinkage_fraction = float(np.clip(sinkage_m / wheel_radius_m, 0, 1))
    risk_index = 100 * (0.7 * sinkage_fraction + 0.3 * (1 - efficiency))
    return {
        "normal_pressure_pa": float(normal_pressure), "sinkage_m": float(sinkage_m),
        "pressure_pa": pressures.tolist(), "sinkage_curve_m": sinkages.tolist(),
        "gross_traction_n": float(gross_traction_n),
        "compaction_resistance_n": float(compaction_resistance_n),
        "drawbar_pull_n": float(drawbar_pull_n),
        "drawbar_pull_efficiency": efficiency,
        "rover_mobility_risk_index": float(np.clip(risk_index, 0, 100)),
        "risk_method": "70% sinkage/radius + 30% (1 − drawbar efficiency); heuristic",
    }
