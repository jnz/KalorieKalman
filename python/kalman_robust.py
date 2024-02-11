import numpy as np
from scipy.stats import chi2
from numpy.linalg import inv, cholesky

def kalman_robust(x, P, dl, R, H, chi2_threshold=None):
    """
    Robust Kalman Filter based on Mahalanobis distance as outlier judging criterion.
    """
    if chi2_threshold is None:
        # Precomputed chi2 table for up to m=20 observations
        # at significance level alpha = 0.05 (5%)
        m = len(dl)
        chi2_threshold = chi2.ppf(0.95, m)

    HPHT = H @ P @ H.T
    P_predicted = HPHT + R
    mahalanobis_dist_squared = dl.T @ inv(P_predicted) @ dl

    if mahalanobis_dist_squared > chi2_threshold:
        f = mahalanobis_dist_squared / chi2_threshold
        Rf = (f - 1) * HPHT + f * R
        outlier = True
    else:
        Rf = R
        outlier = False

    x, P = kalman_takasu(x, P, dl, Rf, H)
    return x, P, outlier

def kalman_takasu(x, P, dl, R, H):
    """
    Kalman Filter equations.
    """
    D = P @ H.T
    S = H @ D + R
    U = cholesky(S)
    E = np.linalg.solve(U, D.T).T  # Equivalent to D/U in MATLAB
    K = np.linalg.solve(U.T, E)  # Equivalent to E/U' in MATLAB
    dx = K @ dl
    x = x + dx
    P = P - E @ E.T
    P = 0.5 * (P + P.T)  # Ensure symmetry
    return x, P

