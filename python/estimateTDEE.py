import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import glob
from kalman_robust import kalman_robust
from readhealthkitcsv import readhealthkitcsv

def estimate_tdee():
    # Constants
    kcal_per_kg_fat = 7700
    tdee_begin_stddev = 600
    std_weight_measurement = 1.15
    tdee_begin = 2000
    kcal_tracking_precision = 100
    pred_noise_kg = 0.05
    pred_noise_tdee = 10.0

    # Loading CSV files
    bodymass_csv_file = glob.glob('*BodyMass*.csv')[0]
    dietenergy_csv_file = glob.glob('*DietaryEnergyConsumed*.csv')[0]

    bodymass = readhealthkitcsv(bodymass_csv_file)
    energyconsumed = readhealthkitcsv(dietenergy_csv_file)

    # Initialize state and covariance matrices
    initial_weight = bodymass['value'].iloc[0]
    filter_state = np.array([[initial_weight], [tdee_begin]])
    filter_state_P = np.diag([std_weight_measurement**2, tdee_begin_stddev**2])
    R = std_weight_measurement**2

    # Sort by startDate if not already sorted
    bodymass.sort_values(by='startDate', inplace=True)
    energyconsumed.sort_values(by='startDate', inplace=True)

    epochs = len(bodymass) - 1
    history_time_in_days = np.zeros(epochs)
    history_weight_raw = np.zeros(epochs)
    history_kcals = np.zeros(epochs)
    history_dt_in_days = np.zeros(epochs)

    # Kalman Filter and RTS Smoother placeholders
    n = 2  # State vector dimension
    rts_history_filter_state_apriori = np.zeros((n, epochs))
    rts_history_filter_state_aposteriori = np.zeros((n, epochs))
    rts_history_filter_state_P_apriori = np.zeros((n, n, epochs))
    rts_history_filter_state_P_aposteriori = np.zeros((n, n, epochs))
    rts_history_phi = np.zeros((n, n, epochs))

    for i in range(epochs):
        # Time delta in days
        if i < epochs - 1:
            dt = (bodymass['startDate'].iloc[i+1] - bodymass['startDate'].iloc[i]).days
        else:
            dt = 1  # Assuming last measurement has at least one day until the "end" of the data

        history_time_in_days[i] = (bodymass['startDate'].iloc[i] - bodymass['startDate'].iloc[0]).days
        history_weight_raw[i] = bodymass['value'].iloc[i]
        history_dt_in_days[i] = dt

        # Energy consumed calculation
        start_date = bodymass['startDate'].iloc[i]
        if i < epochs - 1:
            end_date = bodymass['startDate'].iloc[i+1]
        else:
            end_date = bodymass['startDate'].iloc[i] + pd.Timedelta(days=1)

        kcal_in_interval = energyconsumed[(energyconsumed['startDate'] >= start_date) & (energyconsumed['startDate'] < end_date)]['value'].sum()
        history_kcals[i] = kcal_in_interval

        # Kalman filter update
        H = np.array([[1, 0]])
        dl = np.array([[bodymass['value'].iloc[i] - filter_state[0, 0]]])
        filter_state, filter_state_P, _ = kalman_robust(filter_state, filter_state_P, dl, R, H, chi2_threshold=5.99)  # Example chi2_threshold

        # State transition and prediction
        phi = np.array([[1, -dt / kcal_per_kg_fat], [0, 1]])
        u = kcal_in_interval
        G = np.array([[1 / kcal_per_kg_fat], [0]])
        std_kcal = kcal_tracking_precision * dt
        Qu = std_kcal**2
        Qnoise = np.diag([pred_noise_kg**2 * dt**2, pred_noise_tdee**2 * dt**2])

        filter_state = phi @ filter_state + G * u
        filter_state_P = phi @ filter_state_P @ phi.T + G * Qu * G.T + Qnoise


if __name__ == '__main__':
    estimate_tdee()

