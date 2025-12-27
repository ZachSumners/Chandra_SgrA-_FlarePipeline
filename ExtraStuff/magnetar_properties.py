#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exp_decay(t, A, tau, C):
    """Model: y = A * exp(-t/tau) + C"""
    return A * np.exp(-t / tau) + C

def fit_exp_decay(t, y, yerr=None, allow_offset=True, plot=True):
    """
    Fit y(t) = A*exp(-t/tau) + C using non-linear least squares.

    Parameters
    ----------
    t, y : array-like
        Data (t must be increasing, positive not required). NaNs are ignored.
    yerr : array-like or None
        1σ uncertainties on y; enables proper weighting & uncertainty scaling.
    allow_offset : bool
        If False, fixes C=0 (pure exponential); uses a 2-parameter fit.
    plot : bool
        If True, makes a data+fit plot and residuals (linear axis).

    Returns
    -------
    result : dict
        Keys: A, tau, C, A_err, tau_err, C_err, half_life, half_life_err,
              rchi2, dof, r2, cov, y_fit
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    mask = np.isfinite(t) & np.isfinite(y)
    if yerr is not None:
        yerr = np.asarray(yerr, float)
        mask &= np.isfinite(yerr) & (yerr > 0)
    t, y = t[mask], y[mask]
    if yerr is not None:
        yerr = yerr[mask]

    # Sort by time (curve_fit doesn't need it, but plotting/diagnostics like it)
    order = np.argsort(t)
    t, y = t[order], y[order]
    if yerr is not None:
        yerr = yerr[order]

    # Initial guesses
    C0 = np.min(y) if allow_offset else 0.0
    # Rough amplitude guess from first point minus baseline
    A0 = max(y[0] - C0, (np.max(y) - np.min(y)) * 0.5, 1e-12)
    # Rough tau guess from decay from first to last; fall back to span/3
    span = max(t[-1] - t[0], 1.0)
    drop = max((y[0] - C0) / max(y[-1] - C0, 1e-12), 1.0001)
    tau0 = span / np.log(drop) if np.isfinite(drop) and drop > 1 else span / 3.0

    if allow_offset:
        p0 = [A0, tau0, C0]
        bounds = ([0.0, 0.0, -np.inf], [np.inf, np.inf, np.inf])  # tau>0, A>=0
        model = exp_decay
    else:
        # C fixed to 0: wrap model to 2 parameters
        def model(t, A, tau): return exp_decay(t, A, tau, 0.0)
        p0 = [A0, tau0]
        bounds = ([0.0, 0.0], [np.inf, np.inf])

    # Fit
    popt, pcov = curve_fit(
        model, t, y,
        sigma=yerr if yerr is not None else None,
        p0=p0, bounds=bounds,
        absolute_sigma=True if yerr is not None else False,
        maxfev=20000
    )
    perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.full(len(popt), np.nan)

    # Unpack with/without offset
    if allow_offset:
        A, tau, C = popt
        A_err, tau_err, C_err = perr
    else:
        A, tau = popt
        A_err, tau_err = perr
        C, C_err = 0.0, 0.0

    y_fit = model(t, *popt)

    # Diagnostics
    resid = y - y_fit
    if yerr is not None:
        chi2 = np.sum((resid / yerr)**2)
        dof = max(len(y) - len(popt), 1)
        rchi2 = chi2 / dof
    else:
        dof = max(len(y) - len(popt), 1)
        # Unweighted pseudo-chi2 on unit variance
        rchi2 = np.sum(resid**2) / dof

    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # Half-life and its uncertainty (propagate from tau)
    ln2 = np.log(2.0)
    half_life = tau * ln2
    half_life_err = tau_err * ln2

    if plot:
        # Data + fit
        plt.figure(figsize=(6, 4))
        if yerr is not None:
            plt.errorbar(t, y, yerr=yerr, fmt='o', linestyle='-', capsize=3, label='data2')
        else:
            plt.plot(t, y, 'o', label='data')
        # smooth fit line
        tt = np.linspace(t.min(), t.max(), 400)
        #plt.plot(tt, model(tt, *popt), '-', label=f'fit: A={A:.3g}, τ={tau:.3g}, C={C:.3g}')
        plt.xlabel('Time Since Outburst (April 24, 2013)'); plt.ylabel('Magnetar Flux')
        plt.grid(alpha=0.3); plt.legend(); plt.tight_layout()
        plt.show()

        # Residuals
        plt.figure(figsize=(6, 2.6))
        if yerr is not None:
            plt.errorbar(t, resid, yerr=yerr, fmt='o', capsize=3)
        else:
            plt.plot(t, resid, 'o')
        plt.axhline(0, ls='--')
        plt.xlabel('t'); plt.ylabel('resid')
        plt.grid(alpha=0.3); plt.tight_layout()
        plt.show()

    return {
        "A": A, "tau": tau, "C": C,
        "A_err": A_err, "tau_err": tau_err, "C_err": C_err,
        "half_life": half_life, "half_life_err": half_life_err,
        "rchi2": rchi2, "dof": dof, "r2": r2,
        "cov": pcov, "t": t, "y": y, "y_fit": y_fit
    }

# --- Example usage ---
if __name__ == "__main__":
    np.random.seed(0)
    t = [18.448391203703068, 41.370462962964666, 69.29023148147826, 94.06527777777956, 109.96049768518424, 129.42741898148233, 143.0039699074041, 149.29366898148146, 163.72410879629751, 176.65090277777927, 187.60149305555387, 303.4822569444441, 324.4291550925918, 345.1034027777787, 369.1184259259244, 391.0184606481489, 405.13021990740526, 436.8724421296283, 448.9523263888914, 465.1512847222257, 493.2035416666695, 544.3460879629638, 660.0261226851871, 750.370694444442, 845.4446643518531, 910.2505208333314, 1263.7947916666672, 1269.4474074074096]
    y = [56.715988837210475, 48.632377339332336, 42.56481226832645, 38.15293945634061, 35.619427472745826, 28.93840031675893, 29.482864605779433, 26.818011662728082, 27.182520251765688, 24.617005663285386, 24.50350925949328, 16.262988843606692, 15.162079400821783, 13.823950886694346, 13.20832975680683, 11.72846902009956, 10.63128997335998, 9.243089081864332, 9.142629206991487, 8.474639942209624, 8.083714882987667, 7.225173517238182, 5.390889351954448, 4.146325894851396, 3.281384020169172, 2.624851741238316, 1.3438986720828472, 1.3145863029517355]

    res = fit_exp_decay(t, y, allow_offset=True, plot=True)
    print(f"A = {res['A']:.4g} ± {res['A_err']:.2g}")
    print(f"tau = {res['tau']:.4g} ± {res['tau_err']:.2g}")
    print(f"C = {res['C']:.4g} ± {res['C_err']:.2g}")
    print(f"t1/2 = {res['half_life']:.4g} ± {res['half_life_err']:.2g}")
    print(f"reduced χ² = {res['rchi2']:.3f},  R² = {res['r2']:.4f},  dof = {res['dof']}")