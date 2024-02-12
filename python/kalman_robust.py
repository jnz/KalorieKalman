import numpy as np
from scipy.stats import chi2
from numpy.linalg import inv, cholesky

def kalman_robust(x, P, dl, R, H, chi2_threshold=5.99):
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

    x, P = kalman(x, P, dl, Rf, H)
    return x, P, outlier

def kalman(x, P, dl, R, H):
    """
    Kalman Filter equations.

    x nx1 A priori state vector (size=n) at epoch k
    P nxn Covariance matrix of state vector x at epoch k
    dl mx1 Measurement difference vector (size=m) at epoch k
           dl = measurement - predicted measurement
           dl = y - H*x
    R mxm Covariance matrix of measurement vector y
    H mxn Observation matrix so that y = H*x

    Return value:
    x nx1 A posteriori state vector at epoch k (corrected by measurements)
    P nxn A posteriori covariance of x at epoch k
    """

    S = H.dot(P).dot(H.T) + R  # System uncertainty
    K = P.dot(H.T).dot(np.linalg.inv(S))  # Kalman gain
    x = x + K.dot(dl)
    I = np.eye(P.shape[0])
    P = (I - K.dot(H)).dot(P)

    return x, P
