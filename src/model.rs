/// Empirical CDF: fraction of `sorted_vals` ≤ x. Matches R's `ecdf`.
pub(crate) fn ecdf(sorted_vals: &[f64], x: f64) -> f64 {
    let n = sorted_vals.len();
    if n == 0 {
        return 0.0;
    }
    let pos = sorted_vals.partition_point(|&v| v <= x);
    pos as f64 / n as f64
}

/// `pctremain(x)` from SEACR's R script:
/// `(len_exp - ecdf(expvec)(x)*len_exp) / (len_both - ecdf(both)(x)*len_both)`
///
/// Returns `None` when the denominator is zero (≡ R's `NA`).
fn pctremain(x: f64, sorted_exp: &[f64], sorted_both: &[f64]) -> Option<f64> {
    let len_exp = sorted_exp.len() as f64;
    let len_both = sorted_both.len() as f64;
    let denom = len_both - ecdf(sorted_both, x) * len_both;
    if denom == 0.0 {
        None
    } else {
        let numer = len_exp - ecdf(sorted_exp, x) * len_exp;
        Some(numer / denom)
    }
}

/// Bandwidth for Gaussian KDE: R's `bw.nrd0`.
///
/// `0.9 * lo * n^(-1/5)` where `lo = min(sd, IQR/1.34)`, falling back to `sd` when `IQR == 0`.
pub(crate) fn bw_nrd0(v: &[f64]) -> f64 {
    let n = v.len();
    assert!(n >= 2, "density requires at least 2 observations");
    let mean = v.iter().sum::<f64>() / n as f64;
    let variance = v.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
    let sd = variance.sqrt();
    let iqr = {
        // R's IQR: quantile(type 7) at 0.75 - quantile(type 7) at 0.25
        let mut s = v.to_vec();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        quantile_type7(&s, 0.75) - quantile_type7(&s, 0.25)
    };
    let lo = if iqr == 0.0 { sd } else { sd.min(iqr / 1.34) };
    0.9 * lo * (n as f64).powf(-0.2)
}

/// Mode of R's `density(v)` — the x at `max(y)` on a 512-point Gaussian KDE.
///
/// Uses direct kernel evaluation over a 512-point grid spanning
/// `[min(v) − 3·bw, max(v) + 3·bw]`; same mode index as R's FFT-based `density()`.
pub fn density_mode(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return v.first().copied().unwrap_or(1.0);
    }
    let bw = bw_nrd0(v);
    if bw == 0.0 {
        return v[0];
    }
    let x_lo = v.iter().cloned().fold(f64::INFINITY, f64::min) - 3.0 * bw;
    let x_hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 3.0 * bw;
    let n_grid = 512usize;
    let dx = (x_hi - x_lo) / (n_grid as f64 - 1.0);
    let inv_bw = 1.0 / bw;
    let n_obs = v.len() as f64;
    let mut best_y = f64::NEG_INFINITY;
    let mut best_x = x_lo;
    for i in 0..n_grid {
        let xi = x_lo + i as f64 * dx;
        // Gaussian kernel: (1/n) * sum_j kernel((xi - v_j) / bw) / bw
        let y = v
            .iter()
            .map(|&vj| gauss_kernel((xi - vj) * inv_bw))
            .sum::<f64>()
            / (n_obs * bw);
        if y > best_y {
            best_y = y;
            best_x = xi;
        }
    }
    best_x
}

#[inline]
fn gauss_kernel(u: f64) -> f64 {
    (-0.5 * u * u).exp() * std::f64::consts::FRAC_1_SQRT_2 / std::f64::consts::SQRT_2
}

/// `(1 − fraction)` quantile of a sorted slice, R type-7 definition.
#[must_use]
pub fn quantile_type7(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sorted[0];
    }
    let h = (n as f64 - 1.0) * p;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
}

/// R's `quantile(x, p, type=7)` on an already-sorted vector.
fn quantile_of_sorted(sorted: &[f64], p: f64) -> f64 {
    quantile_type7(sorted, p)
}

/// Compute thresholds x0, z0, d0 from block AUC and num-interval vectors.
///
/// Faithful port of the pctremain/z0/spurious-correction logic in `SEACR_1.3.R`.
/// Returns `(x0, z0, d0)`.
pub(crate) fn compute_thresholds(
    expvec: &[f64],
    ctrlvec: &[f64],
    expmax: &[f64],
    ctrlmax: &[f64],
) -> (f64, f64, f64) {
    let mut sorted_exp = expvec.to_vec();
    sorted_exp.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut sorted_ctrl = ctrlvec.to_vec();
    sorted_ctrl.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut both: Vec<f64> = expvec.iter().chain(ctrlvec.iter()).copied().collect();
    both.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let x_unique: Vec<f64> = {
        let mut u = both.clone();
        u.dedup();
        u
    };

    let pr = |x: f64| pctremain(x, &sorted_exp, &both);

    // x0: x maximising pctremain among pctremain < 1
    let x0 = {
        let cands: Vec<f64> = x_unique
            .iter()
            .copied()
            .filter(|&x| pr(x).is_some_and(|v| v < 1.0))
            .collect();
        let best = cands
            .iter()
            .copied()
            .max_by(|&a, &b| pr(a).unwrap().partial_cmp(&pr(b).unwrap()).unwrap());
        best.unwrap_or_else(|| *x_unique.last().unwrap_or(&1.0))
    };

    let z0 = compute_z0(x0, &x_unique, pr);
    let (x0, z0) = spurious_correction(x0, z0, &x_unique, pr);
    let d0 = compute_d0(expmax, ctrlmax);

    (x0, z0, d0)
}

/// Compute z0 from x0 using the midpoint-of-curve logic in SEACR_1.3.R.
fn compute_z0(x0: f64, x_unique: &[f64], pr: impl Fn(f64) -> Option<f64>) -> f64 {
    let z: Vec<f64> = x_unique.iter().copied().filter(|&x| x <= x0).collect();
    if z.is_empty() {
        return x0;
    }
    let pr_z_min = z
        .iter()
        .filter_map(|&x| pr(x))
        .fold(f64::INFINITY, f64::min);
    let pr_x0 = pr(x0).unwrap_or(0.0);
    let midpoint = (pr_x0 + pr_z_min) / 2.0;
    let z2 = *z
        .iter()
        .min_by(|&&a, &&b| {
            let da = (pr(a).unwrap_or(0.0) - midpoint).abs();
            let db = (pr(b).unwrap_or(0.0) - midpoint).abs();
            da.partial_cmp(&db).unwrap()
        })
        .unwrap_or(&x0);

    if x0 == z2 {
        // Added 7/15/19 to avoid omitting z when x0==z2
        return x0;
    }
    let z_filtered: Vec<f64> = z.iter().copied().filter(|&x| x > z2).collect();
    if z_filtered.is_empty() {
        return x0;
    }
    let z_max = z_filtered.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let z_min = z_filtered.iter().copied().fold(f64::INFINITY, f64::min);
    let target = z_max - 0.5 * (z_max - z_min);
    *z_filtered
        .iter()
        .min_by(|&&a, &&b| (a - target).abs().partial_cmp(&(b - target).abs()).unwrap())
        .unwrap_or(&x0)
}

/// Spurious-threshold correction: if the sub-curve's max pctremain is ≥ 95%
/// of the full curve's max pctremain, replace x0/z0 with the sub-curve's.
fn spurious_correction(
    x0: f64,
    z0_initial: f64,
    x_unique: &[f64],
    pr: impl Fn(f64) -> Option<f64>,
) -> (f64, f64) {
    if x_unique.len() < 2 {
        return (x0, z0_initial);
    }
    let pairs: Vec<(f64, f64, f64)> = x_unique
        .windows(2)
        .filter_map(|w| {
            let thresh = w[0];
            let pr_a = pr(w[0]);
            let pr_b = pr(w[1]);
            let diff = match (pr_a, pr_b) {
                (Some(a), Some(b)) => (b - a).abs(),
                (Some(a), None) => a.abs(),
                (None, Some(b)) => b.abs(),
                (None, None) => return None,
            };
            let pct = pr(thresh)?;
            Some((thresh, pct, diff))
        })
        .collect::<Vec<_>>()
        .into_iter()
        .filter(|&(_, _, diff)| !diff.is_nan())
        .collect::<Vec<_>>();

    let diffs: Vec<f64> = pairs.iter().map(|&(_, _, d)| d).collect();
    if diffs.is_empty() {
        return (x0, z0_initial);
    }
    let mut test3 = 0.99f64;
    let mut nines = 2usize;
    loop {
        let q = quantile_of_sorted(
            &{
                let mut s = diffs.clone();
                s.sort_by(|a, b| a.partial_cmp(b).unwrap());
                s
            },
            test3,
        );
        if q > 0.0 {
            break;
        }
        nines += 1;
        test3 = 1.0 - 10f64.powi(-(nines as i32));
    }
    let q_thresh = quantile_of_sorted(
        &{
            let mut s = diffs.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            s
        },
        test3,
    );

    let a: Vec<f64> = pairs
        .iter()
        .filter(|&&(_, _, d): &&(f64, f64, f64)| d != 0.0 && d < q_thresh)
        .map(|&(t, _, _)| t)
        .collect();
    if a.is_empty() {
        return (x0, z0_initial);
    }
    let mut sorted_a = a.clone();
    sorted_a.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let a_cands: Vec<f64> = a
        .iter()
        .copied()
        .filter(|&x| pr(x).is_some_and(|v| v < 1.0))
        .collect();
    if a_cands.is_empty() {
        return (x0, z0_initial);
    }
    let a0 = *a_cands
        .iter()
        .max_by(|&&p, &&q| pr(p).unwrap().partial_cmp(&pr(q).unwrap()).unwrap())
        .unwrap();

    let b0 = compute_z0(a0, &sorted_a, &pr);

    let max_pr_a = a_cands
        .iter()
        .filter_map(|&x| pr(x))
        .fold(f64::NEG_INFINITY, f64::max);
    let max_pr_x = x_unique
        .iter()
        .filter_map(|&x| pr(x))
        .filter(|&v| v < 1.0)
        .fold(f64::NEG_INFINITY, f64::max);

    if max_pr_x > 0.0 && max_pr_a / max_pr_x > 0.95 {
        (a0, b0)
    } else {
        (x0, z0_initial)
    }
}

/// d0: min num-intervals value where `pctremain2(x) > 1`.
///
/// `pctremain2(x) = 1 - (ecdf(expmax)(x) - ecdf(ctrlmax)(x))`
pub(crate) fn compute_d0(expmax: &[f64], ctrlmax: &[f64]) -> f64 {
    if expmax.is_empty() && ctrlmax.is_empty() {
        return 1.0;
    }
    let mut sorted_exp = expmax.to_vec();
    sorted_exp.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut sorted_ctrl = ctrlmax.to_vec();
    sorted_ctrl.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let both2: Vec<f64> = {
        let mut v: Vec<f64> = expmax.iter().chain(ctrlmax.iter()).copied().collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v.dedup();
        v
    };

    let pctremain2 = |x: f64| 1.0 - (ecdf(&sorted_exp, x) - ecdf(&sorted_ctrl, x));

    let candidates: Vec<f64> = both2
        .iter()
        .copied()
        .filter(|&x| pctremain2(x) > 1.0)
        .collect();
    if candidates.is_empty() {
        1.0
    } else {
        candidates.iter().cloned().fold(f64::INFINITY, f64::min)
    }
}

/// Knee-point detection for the norm-mode cutoff value.
///
/// Implements the `dist2d` + `expframe`/`ctrlframe` logic from SEACR_1.3.R.
/// Returns the AUC cutoff value below which density mode is computed.
pub(crate) fn knee_value(vec: &[f64]) -> f64 {
    let n = vec.len();
    let ninety_pct_idx = (0.9 * n as f64) as usize;
    if ninety_pct_idx == 0 {
        return *vec
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&1.0);
    }
    let mut sorted = vec.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p90 = sorted[ninety_pct_idx.min(n - 1)];
    let max_val = *sorted.last().unwrap();

    let mut desc = sorted.clone();
    desc.reverse();

    // dist2d(c(count,quant), 0, 1) = count - quant per row (R derivation in SEACR_1.3.R)
    let mut frame: Vec<(f64, f64, f64, f64)> = desc
        .iter()
        .enumerate()
        .map(|(i, &val)| {
            let count = 1.0 - (i as f64) / (n as f64 - 1.0).max(1.0);
            let quant = if max_val == 0.0 { 0.0 } else { val / max_val };
            let diff = (count - quant).abs();
            let _dist = count - quant;
            (count, quant, val, diff)
        })
        .collect();

    let max_diff = frame
        .iter()
        .map(|&(_, _, _, d)| d)
        .fold(f64::NEG_INFINITY, f64::max);
    let threshold_diff = 0.9 * max_diff;
    frame.retain(|&(_, _, _, d)| d > threshold_diff);

    if frame.is_empty() {
        return p90.max(*sorted.last().unwrap());
    }

    let knee_val = frame
        .iter()
        .max_by(|&(ca, qa, _, _), &(cb, qb, _, _)| {
            let da = ca - qa;
            let db = cb - qb;
            da.partial_cmp(&db).unwrap()
        })
        .map(|&(_, _, val, _)| val)
        .unwrap_or(max_val);

    if knee_val > p90 { knee_val } else { p90 }
}

/// For numeric-threshold mode: min value where `1 - ecdf(v)(x) <= frac`.
pub(crate) fn numeric_threshold_min(sorted: &[f64], frac: f64) -> f64 {
    let n = sorted.len() as f64;
    sorted
        .iter()
        .enumerate()
        .find(|&(i, _)| {
            let percentile = 1.0 - (i as f64 + 1.0) / n;
            percentile <= frac
        })
        .map(|(_, &v)| v)
        .unwrap_or(*sorted.last().unwrap_or(&0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_type7_endpoints() {
        let s = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(quantile_type7(&s, 0.0), 1.0);
        assert_eq!(quantile_type7(&s, 1.0), 5.0);
        assert_eq!(quantile_type7(&s, 0.5), 3.0);
    }

    #[test]
    fn density_mode_matches_r() {
        // ctrl subvec from golden fixture: ctrlvec[ctrlvec <= 4.70167]
        let ctrl_sub = vec![
            0.733204, 4.70167, 0.824854, 1.97965, 3.79433, 2.6212, 0.641553, 4.15178, 3.46439,
            2.56621, 1.28311, 4.21592, 3.97763, 4.41755, 1.01732, 0.659884, 0.733204, 2.62121,
            0.797359, 4.61919, 2.23627, 1.09981, 2.50206, 0.339107, 3.92264, 0.403262, 3.776,
            1.77802, 0.632388, 0.559068,
        ];
        // R: density(ctrl_sub)$x[which.max(density(ctrl_sub)$y)] ≈ 0.8467413
        let mode = density_mode(&ctrl_sub);
        assert!(
            (mode - 0.8467413).abs() < 0.005,
            "ctrl density mode {mode} != expected ~0.8467"
        );

        let exp_sub = vec![
            152.116, 130.671, 32.0244, 12.5096, 16.4411, 6.50497, 29.4511, 143.467, 74.6285,
            30.3089, 51.3964, 24.5902, 26.9492, 71.1972, 32.4533, 11.4373, 13.7248, 61.9045,
            7.93464, 17.013, 62.9768, 13.8677, 17.2275, 8.72095, 3.43119, 18.4427, 43.4617,
            53.1121, 2.64488,
        ];
        // R: mode ≈ 17.18486
        let exp_mode = density_mode(&exp_sub);
        assert!(
            (exp_mode - 17.18486).abs() < 0.1,
            "exp density mode {exp_mode} != expected ~17.18"
        );
    }

    #[test]
    fn bw_nrd0_matches_r() {
        // R: bw.nrd0(ctrl_sub) ≈ 0.6912042
        let ctrl_sub = vec![
            0.733204, 4.70167, 0.824854, 1.97965, 3.79433, 2.6212, 0.641553, 4.15178, 3.46439,
            2.56621, 1.28311, 4.21592, 3.97763, 4.41755, 1.01732, 0.659884, 0.733204, 2.62121,
            0.797359, 4.61919, 2.23627, 1.09981, 2.50206, 0.339107, 3.92264, 0.403262, 3.776,
            1.77802, 0.632388, 0.559068,
        ];
        let bw = bw_nrd0(&ctrl_sub);
        assert!(
            (bw - 0.6912042).abs() < 0.0001,
            "bw {bw} != expected ~0.6912"
        );
    }
}
