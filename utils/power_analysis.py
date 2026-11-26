"""
Post-hoc power analysis for small-sample hypothesis tests.

With n=3 targets vs n=8 controls, Mann-Whitney U has very limited
statistical power. This module provides honest power estimates
via simulation, so the paper can transparently report what effect
sizes the study design could realistically detect.
"""

import numpy as np
from scipy import stats


def mannwhitney_power(n1, n2, effect_size_d, alpha=0.05, n_sim=10000, seed=42):
    """Estimate power of one-sided Mann-Whitney U via simulation.

    Draws n_sim pairs of normal samples with the given Cohen's d,
    runs MWU on each, and reports the fraction of significant results.
    """
    rng = np.random.RandomState(seed)
    sig_count = 0
    for _ in range(n_sim):
        x = rng.normal(effect_size_d, 1, n1)
        y = rng.normal(0, 1, n2)
        try:
            _, p = stats.mannwhitneyu(x, y, alternative="greater")
            if p < alpha:
                sig_count += 1
        except ValueError:
            pass
    return sig_count / n_sim


def fisher_power(p1, p2, n1, n2, alpha=0.05, n_sim=10000, seed=42):
    """Estimate power of one-sided Fisher's exact test via simulation."""
    rng = np.random.RandomState(seed)
    sig_count = 0
    for _ in range(n_sim):
        a = rng.binomial(n1, p1)
        c = rng.binomial(n2, p2)
        table = [[a, n1 - a], [c, n2 - c]]
        try:
            _, p = stats.fisher_exact(table, alternative="greater")
            if p < alpha:
                sig_count += 1
        except ValueError:
            pass
    return sig_count / n_sim


def minimum_detectable_effect(n1, n2, alpha=0.05, target_power=0.80,
                               d_range=None, n_sim=5000, seed=42):
    """Find minimum Cohen's d detectable at target_power with MWU."""
    if d_range is None:
        d_range = np.arange(0.5, 5.1, 0.25)

    for d in d_range:
        power = mannwhitney_power(n1, n2, d, alpha, n_sim, seed)
        if power >= target_power:
            return round(float(d), 2), round(power, 3)

    return float(d_range[-1]), round(
        mannwhitney_power(n1, n2, float(d_range[-1]), alpha, n_sim, seed), 3
    )


def full_power_report(n1, n2, observed_d, alpha=0.05, n_sim=10000, seed=42):
    """Complete power analysis report for a MWU test."""
    achieved_power = mannwhitney_power(n1, n2, abs(observed_d), alpha, n_sim, seed)

    d_grid = np.arange(0.5, 5.1, 0.25)
    power_curve = {}
    for d in d_grid:
        power_curve[round(float(d), 2)] = round(
            mannwhitney_power(n1, n2, d, alpha, 5000, seed), 3
        )

    mde_d, mde_power = minimum_detectable_effect(n1, n2, alpha, 0.80, d_grid, 5000, seed)

    return {
        "n_group1": n1,
        "n_group2": n2,
        "observed_effect_size_d": round(observed_d, 3),
        "achieved_power": round(achieved_power, 3),
        "alpha": alpha,
        "min_detectable_effect_d_at_80pct": mde_d,
        "interpretation": _interpret(n1, n2, achieved_power, mde_d, observed_d),
        "power_curve": power_curve,
    }


def _interpret(n1, n2, achieved_power, mde_d, observed_d):
    if achieved_power < 0.30:
        qualifier = "severely underpowered"
    elif achieved_power < 0.50:
        qualifier = "underpowered"
    elif achieved_power < 0.80:
        qualifier = "moderately powered"
    else:
        qualifier = "adequately powered"

    return (
        f"With n={n1} vs n={n2}, this test is {qualifier} "
        f"(achieved power={achieved_power:.1%}). "
        f"The minimum detectable effect at 80% power is Cohen's d={mde_d}, "
        f"which corresponds to a {'very ' if mde_d > 2.0 else ''}"
        f"large effect. The observed d={observed_d:.2f} "
        f"{'falls below' if abs(observed_d) < mde_d else 'exceeds'} this threshold. "
        f"Non-significant results should be interpreted as inconclusive, "
        f"not as evidence of no difference."
    )
