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
    questionable_measurements = 0; # keep track of prediction vs. measurement mismatches

    # Sort by startDate if not already sorted
    bodymass.sort_values(by='startDate', inplace=True)
    energyconsumed.sort_values(by='startDate', inplace=True)

    epochs = len(bodymass) - 1
    history_time_in_days = np.zeros(epochs)
    history_weight_raw = np.zeros(epochs)
    history_kcals = np.zeros(epochs)
    history_dt_in_days = np.zeros(epochs)
    # Kalman Filter and RTS Smoother placeholders
    n = filter_state.shape[0]  # State vector dimension
    rts_history_filter_state_apriori = np.zeros((n, epochs))
    rts_history_filter_state_aposteriori = np.zeros((n, epochs))
    rts_history_filter_state_P_apriori = np.zeros((n, n, epochs))
    rts_history_filter_state_P_aposteriori = np.zeros((n, n, epochs))
    rts_history_phi = np.zeros((n, n, epochs))

    for i in range(epochs):
        # bookkeeping for plots
        current_body_mass = bodymass['value'].iloc[i] # mass at epoch i
        history_time_in_days[i] = (bodymass['startDate'].iloc[i] - bodymass['startDate'].iloc[0]).days
        history_weight_raw[i] = current_body_mass

        # Kalman Filter Fusion / Correction with current body mass measurement
        # --------------------------------------------------------------------
        rts_history_filter_state_apriori[:, i] = filter_state.flatten()
        rts_history_filter_state_P_apriori[:,:, i] = filter_state_P

        H = np.array([[1, 0]]) # Kalman Filter Design matrix
        # dl is the delta between the measurement and the expected body mass
        dl = np.array([[current_body_mass - filter_state[0, 0]]])
        filter_state, filter_state_P, outlier = kalman_robust(filter_state, filter_state_P, dl, R, H, chi2_threshold=5.99)
        if outlier:
            questionable_measurements += 1

        rts_history_filter_state_aposteriori[:, i] = filter_state.flatten()
        rts_history_filter_state_P_aposteriori[:,:, i] = filter_state_P

        # Predict to next epoch
        # ---------------------
        # Energy consumed calculation
        start_date = bodymass['startDate'].iloc[i]
        end_date = bodymass['startDate'].iloc[i+1]
        dt_date = (end_date - start_date)
        dt = dt_date.total_seconds()/(3600.0*24) # prediction step in days (float)
        assert(dt >= 0.0)
        history_dt_in_days[i] = dt # bookkeeping for plots

        # sum up all calories consumed in the prediction horizon
        kcal_in_interval = energyconsumed[(energyconsumed['startDate'] >= start_date) & (energyconsumed['startDate'] < end_date)]['value'].sum()
        history_kcals[i] = kcal_in_interval # bookkeeping for plots

        # try to  highlight periods without proper food tracking
        prediction_penalty_kg = 0
        if dt >= 1.0: # if for a longer period of time,
            if kcal_in_interval < 0.5*filter_state[1, 0]*dt: # the tracked calories is way below the TDEE, print a note to the console
                print("Note: Low food energy consumption on day %i" % round(history_time_in_days[i]))
                prediction_penalty_kg = 0.25; # assume that the predicted body mass is a bit more inaccurate

        # Kalman state transition and prediction for next epoch
        phi = np.array([[1, -dt / kcal_per_kg_fat], [0, 1]])
        u = kcal_in_interval
        G = np.array([[1 / kcal_per_kg_fat], [0]])
        std_kcal = kcal_tracking_precision * dt
        Qu = std_kcal**2
        Qnoise = np.diag([((pred_noise_kg + prediction_penalty_kg)*dt)**2, (pred_noise_tdee*dt)**2])

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
    plt.gca().set_prop_cycle(None)  # Reset the color cycle to reuse the first color for y-axis right
    ax1 = plt.gca()

    # Create secondary y-axis for TDEE
    ax2 = plt.gca().twinx()
    ax2.plot(history_time_in_days, filter_state_rts[1, :], 'r-', label='TDEE (kcal)')
    # ax2.plot(history_time_in_days, filter_state_rts[1, :] + std_tdee, 'k--', label='1 sigma band TDEE')
    # ax2.plot(history_time_in_days, filter_state_rts[1, :] - std_tdee, 'k--')
    ax2.fill_between(history_time_in_days,
                     filter_state_rts[1, :] - std_tdee,
                     filter_state_rts[1, :] + std_tdee,
                     color='red', alpha=0.15, label='1 σ error band')
    ax2.set_ylabel('TDEE (kcal)')
    # Adjust the y-limit for TDEE
    ax2.set_ylim(np.min(filter_state_rts[1, :] - std_tdee) * 0.9, np.max(filter_state_rts[1, :] + std_tdee) * 1.1)
    # Because we are using twin axes, we need to manually construct the legend to include all lines

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='lower right')

    plt.show()

if __name__ == '__main__':
    estimate_tdee()

