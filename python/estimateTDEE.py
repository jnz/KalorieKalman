import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import glob
from kalman_robust import kalman_robust
from readhealthkitcsv import readhealthkitcsv
from kalman_rts import kalman_rts

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

    bodymass_entries_before = len(bodymass)
    bodymass.drop_duplicates(subset=['startDate'], inplace=True)
    bodymass_entries_after = len(bodymass)
    print("Removed %i duplicated body mass entries." % (bodymass_entries_before - bodymass_entries_after))

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
        dt_date = (bodymass['startDate'].iloc[i+1] - bodymass['startDate'].iloc[i])
        dt = dt_date.total_seconds()/(3600.0*24)

        history_time_in_days[i] = (bodymass['startDate'].iloc[i] - bodymass['startDate'].iloc[0]).days
        history_weight_raw[i] = bodymass['value'].iloc[i]
        history_dt_in_days[i] = dt

        # Energy consumed calculation
        start_date = bodymass['startDate'].iloc[i]
        end_date = bodymass['startDate'].iloc[i+1]

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

        rts_history_phi[:,:,i] = phi

    # smooth data backwards with the Rauch Tung Striebel smoother
    filter_state_rts, filter_state_rts_P = kalman_rts(
        rts_history_filter_state_apriori,
        rts_history_filter_state_aposteriori,
        rts_history_filter_state_P_apriori,
        rts_history_filter_state_P_aposteriori,
        rts_history_phi)

    # Calculate the standard deviation for TDEE
    cov_tdee = filter_state_rts_P[1, 1, :]  # Extract the covariance of TDEE over time
    std_tdee = np.sqrt(cov_tdee)  # Standard deviation is the square root of the variance

    # Plotting
    plt.figure()
    plt.grid(True)
    plt.xlabel('Time (days)')
    plt.title(f'TDEE estimation: {filter_state_rts[1, -1]:.0f} kcal')

    # Plot mass
    plt.plot(history_time_in_days, filter_state_rts[0, :], 'b', label='Body Mass (kg)')
    plt.ylabel('Mass (kg)')
    plt.gca().set_prop_cycle(None)  # Reset the color cycle to reuse the first color for yyaxis right

    # Create secondary y-axis for TDEE
    ax2 = plt.gca().twinx()
    ax2.plot(history_time_in_days, filter_state_rts[1, :], 'r-', label='TDEE (kcal)')
    ax2.plot(history_time_in_days, filter_state_rts[1, :] + std_tdee, 'k--', label='1 sigma band TDEE')
    ax2.plot(history_time_in_days, filter_state_rts[1, :] - std_tdee, 'k--')
    ax2.set_ylabel('TDEE (kcal)')

    # Adjust the y-limit for TDEE
    ax2.set_ylim(0, np.max(filter_state_rts[1, :] + std_tdee) * 1.2)

    # Because we are using twin axes, we need to manually construct the legend to include all lines
    lines, labels = plt.gca().get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='best')

    plt.show()

if __name__ == '__main__':
    estimate_tdee()

