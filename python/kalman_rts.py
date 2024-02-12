import numpy as np

def kalman_rts(x_apriori, x_aposteriori, P_apriori, P_aposteriori, phi):
    """
    Rauch-Tung-Striebel Smoother.
    Source:
    - H. Rauch, F. Tung, and C. Striebel. Maximum likelihood estimates of
      lineardynamic systems. AIAA Journal, 3(8):1445 – 1450, 1965
    - P. D. Groves, Principles of GNSS, Inertial, and Multisensor
      Integrated Navigation Systems. Norwood, MA: Artech House, 2008.

    n: Number of elements in column state vector
    N: Number of total epochs

    Parameters:
    - x_apriori: A priori state vector (n by 1 by N) at epoch 1..N
    - x_aposteriori: A posteriori state vector (n by 1 by N) at epoch 1..N
    - P_apriori: A priori covariance matrix of state vector x at epoch 1..N
    - P_aposteriori: A posteriori covariance matrix of state vector x at epoch 1..N
    - phi: State transition matrix at epoch 1..N
    """
    N = x_apriori.shape[1]  # Number of total epochs
    n = x_apriori.shape[0]  # Number of elements in column state vector

    # Initialize smoothed state and covariance matrices
    x_smooth = np.zeros((n, N))
    P_smooth = np.zeros((n, n, N))

    # Set the last smoothed state to the last filtered state
    x_smooth[:, N-1] = x_aposteriori[:, N-1]
    P_smooth[:, :, N-1] = P_aposteriori[:, :, N-1]

    # Iterate backwards through time
    for k in range(N-2, -1, -1):
        K = P_aposteriori[:, :, k] @ phi[:, :, k].T @ np.linalg.inv(P_apriori[:, :, k+1])
        x_smooth_diff = x_smooth[:, k+1] - x_apriori[:, k+1]
        x_smooth[:, k] = x_aposteriori[:, k] + K @ x_smooth_diff
        P_smooth_diff = P_smooth[:, :, k+1] - P_apriori[:, :, k+1]
        P_smooth[:, :, k] = P_aposteriori[:, :, k] + K @ P_smooth_diff @ K.T

    return x_smooth, P_smooth
