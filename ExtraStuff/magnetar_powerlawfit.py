#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --------------------------------------------------------------------
# Define power-law model
def powerlaw(x, A, alpha):
    """Power law: y = A * x**(-alpha)"""
    return A * np.power(x, -alpha)

# --------------------------------------------------------------------
def fit_powerlaw(x, y, yerr=None, plot=True, logspace=True):
    """
    Fit y = A * x^-alpha using nonlinear least squares.
    Parameters
    ----------
    x, y : array-like
        Data to fit (must be positive).
    yerr : array-like or None
        Optional uncertainties on y.
    plot : bool
        Whether to make diagnostic plot.
    logspace : bool
        Whether to fit in log space (linear regression on log(y) vs log(x))
        if yerr not provided.
    Returns
    -------
    A, alpha : float
        Best-fit parameters.
    perr : ndarray
        1-sigma uncertainties on parameters.
    """
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)
    mask = (x > 0) & (y > 0)
    x, y = x[mask], y[mask]

    if logspace and yerr is None:
        # Linear fit in log-log space
        lx, ly = np.log10(x), np.log10(y)
        coeffs, cov = np.polyfit(lx, ly, 1, cov=True)
        slope, intercept = coeffs
        perr = np.sqrt(np.diag(cov))
        alpha = -slope
        A = 10**intercept
        if plot:
            plt.figure(figsize=(5,4))
            plt.loglog(x, y, 'o', label='data')
            plt.loglog(x, A*x**(-alpha), '-', label=f'A={A:.3g}, α={alpha:.3f}')
            plt.xlabel('x'); plt.ylabel('y'); plt.legend(); plt.grid(True, which='both', alpha=0.3)
            plt.tight_layout(); plt.show()
        return A, alpha, perr
    else:
        # Nonlinear fit (accounts for yerr)
        p0 = [np.median(y), 1.0]
        popt, pcov = curve_fit(powerlaw, x, y, sigma=yerr, p0=p0, absolute_sigma=True if yerr is not None else False)
        perr = np.sqrt(np.diag(pcov))
        A, alpha = popt

        if plot:
            plt.figure(figsize=(5,4))
            plt.loglog(x, y, 'o', label='data')
            plt.loglog(x, powerlaw(x, *popt), '-', label=f'A={A:.3g}, α={alpha:.3f}')
            plt.xlabel('x'); plt.ylabel('y'); plt.legend(); plt.grid(True, which='both', alpha=0.3)
            plt.tight_layout(); plt.show()
        return A, alpha, perr

# --------------------------------------------------------------------
if __name__ == "__main__":
    # Example usage
    np.random.seed(0)
    x = [56.715988837210475, 48.632377339332336, 42.56481226832645, 38.15293945634061, 35.619427472745826, 28.93840031675893, 29.482864605779433, 26.818011662728082, 27.182520251765688, 24.617005663285386, 24.50350925949328, 16.262988843606692, 15.162079400821783, 13.823950886694346, 13.20832975680683, 11.72846902009956, 10.63128997335998, 9.243089081864332, 9.142629206991487, 8.474639942209624, 8.083714882987667, 7.225173517238182, 5.390889351954448, 4.146325894851396, 3.281384020169172, 2.624851741238316, 1.3438986720828472, 1.3145863029517355]
    y = [18.448391203703068, 41.370462962964666, 69.29023148147826, 94.06527777777956, 109.96049768518424, 129.42741898148233, 143.0039699074041, 149.29366898148146, 163.72410879629751, 176.65090277777927, 187.60149305555387, 303.4822569444441, 324.4291550925918, 345.1034027777787, 369.1184259259244, 391.0184606481489, 405.13021990740526, 436.8724421296283, 448.9523263888914, 465.1512847222257, 493.2035416666695, 544.3460879629638, 660.0261226851871, 750.370694444442, 845.4446643518531, 910.2505208333314, 1263.7947916666672, 1269.4474074074096]

    A, alpha, perr = fit_powerlaw(x, y, plot=True)
    print(f"Best fit: A = {A:.3f} ± {perr[1]:.3f}, α = {alpha:.3f} ± {perr[0]:.3f}")