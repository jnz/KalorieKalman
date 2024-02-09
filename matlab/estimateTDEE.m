function [] = estimateTDEE()

close all;

% <const to tune>
kcal_per_kg_fat         = 7700;   % kcal in a kg of fat (constant) Source: Mayo Clinic: "Counting calories: Get back to weight-loss basics." Published: Feb. 12, 2021
TDEE_begin_stddev       = 600;    % kcal std. dev. of initial TDEE est.
std_weight_measurement  = 1.15;   % kg std. deviation of weight measurement. Normal daily weight fluctuation can range from around 0.5 - 1.8 kg due to factors such as fluid balance, sodium intake, carbohydrate intake, bowel movements, and hormonal changes,
TDEE_begin              = 2000;   % kcal TDEE at the start. As the filter runs backwards as well, the initial value is not that critical.
kcal_tracking_precision = 100;    % kcal/day tracking error
pred_noise_kg           = 0.05;   % kg/day body mass prediction noise
pred_noise_tdee         = 10.0;   % kcal/day prediction noise (this is supposed to be the basal metabolic rate (BMR) variation from day to day, which is assumed to be somewhat stable but subject to metabolic adaptions over a longer time).
% </const to tune>

% <input data>
% List all files in the current directory that have 'BodyMass' in their name and a .csv extension
files = dir('*BodyMass*.csv');
if length(files) < 1
    fprintf('No body mass .csv file found')
    return
end
bodymass_csv_file = files(1).name;

files = dir('*DietaryEnergyConsumed*.csv');
if length(files) < 1
    fprintf('No dietary energy consumed .csv file found')
    return
end
dietenergyconsumed_csv_file = files(1).name;
% </input data>

fprintf('Loading body mass from %s... ', bodymass_csv_file);
bodymass = readhealthkitcsv(bodymass_csv_file);
fprintf('%i entries.\n', size(bodymass, 1));
fprintf('Loading energy consumed from %s... ', dietenergyconsumed_csv_file);
energyconsumed = readhealthkitcsv(dietenergyconsumed_csv_file);
fprintf('%i entries.\n', size(energyconsumed, 1));

% Filter out duplicates from bodymass:
epochs_pre = size(bodymass,1) - 1;
[~, ia] = unique(bodymass.startDate);
bodymass = bodymass(ia, :);
epochs = size(bodymass,1) - 1;
if epochs ~= epochs_pre
    fprintf('Filtered %i body mass measurements down to %i unique measurements.\n', epochs_pre+1, epochs);
end

if epochs < 5
    fprint('Not enough measurements');
    return
end

% state vector:
% x =  [ body mass in kg
%        TDEE (total daily energy expenditure in kcal ]
filter_state = [bodymass.value(1); TDEE_begin];
% covariance of state vector:
filter_state_P = diag([std_weight_measurement^2, TDEE_begin_stddev^2]);
R = std_weight_measurement^2; % covariance of weight measurement
questionable_measurements = 0; % keep track of prediction vs. measurement mismatches

history_time_in_days = zeros(epochs, 1);
history_weight_raw = zeros(epochs, 1);
history_kcals = zeros(epochs, 1);
history_dt_in_days = zeros(epochs, 1);
fprintf('Analyzing data between %s and %s.\n', ...
         datestr(bodymass.startDate(1)), ...
         datestr(bodymass.startDate(end)) ); %#ok<DATST>

% for RTS smoother
n = length(filter_state);
rts_history_filter_state_apriori = zeros(n, epochs);
rts_history_filter_state_aposteriori = zeros(n, epochs);
rts_history_filter_state_P_apriori = zeros(n, n, epochs);
rts_history_filter_state_P_aposteriori = zeros(n, n, epochs);
rts_history_phi = zeros(n, n, epochs);

for i=1:epochs

    rts_history_filter_state_apriori(:, i) = filter_state;
    rts_history_filter_state_P_apriori(:,:, i) = filter_state_P;

    % Fusion with latest weight/mass measurement
    % ------------------------------------------
    H = [1 0];
    dl = bodymass.value(i) - H*filter_state; % weight measured - predicted weight
    [filter_state, filter_state_P, outlier] = ...
        kalman_robust(filter_state, ...
                      filter_state_P, ...
                      dl, ...
                      R, H);

    if outlier == true
        questionable_measurements = questionable_measurements + 1;
    end

    rts_history_filter_state_aposteriori(:, i) = filter_state;
    rts_history_filter_state_P_aposteriori(:,:, i) = filter_state_P;

    % bookkeeping for plots
    history_time_in_days(i) = ...
        days(bodymass.startDate(i)-bodymass.startDate(1));
    history_weight_raw(i) = bodymass.value(i);

    % prediction to next epoch (d2)
    % -----------------------------
    d1 = bodymass.startDate(i);
    d2 = bodymass.startDate(i+1);
    d = d2 - d1; % duration
    dt = days(d); % time delta in days
    if dt < 0
        error('Data entries are not sorted by date.')
    end
    history_dt_in_days(i) = dt;

    % energy consumed from food in current time period
    index_energy_entries_in_time_interval = ...
        (energyconsumed.startDate >= d1 & energyconsumed.startDate < d2);
    kcal_in_interval = energyconsumed(index_energy_entries_in_time_interval, {'value'});
    sum_kcal = sum(kcal_in_interval.value);
    history_kcals(i) = sum_kcal;

    % try to  highlight periods without proper food tracking
    prediction_penalty_kg = 0;
    if (dt >= 1) % if for a longer period of time,
        if (sum_kcal < 0.5*filter_state(2)*dt) % the tracked calories is way below the TDEE, print a note to the console
            fprintf('Note: Low food energy consumption on day %i, i.e. between %s and %s. KCAL=%.0f (~ %.0f KCAL/day).\n', round(history_time_in_days(i)), datestr(d1), datestr(d2), sum_kcal, sum_kcal/dt); %#ok<DATST>
            prediction_penalty_kg = 0.25; % assume that the predicted body mass is a bit more inaccurate
        end
    end

    phi = [1   -1/kcal_per_kg_fat * dt; 0 1 ];
    u = sum_kcal; % 1x1 vector

    std_kcal = kcal_tracking_precision*dt;
    Qu = std_kcal^2; % 1x1 matrix
    G = [1/kcal_per_kg_fat; 0];
    Qnoise = diag([((pred_noise_kg + prediction_penalty_kg)*dt)^2 (pred_noise_tdee*dt)^2]);

    filter_state = phi*filter_state + G*u;
    filter_state_P = phi*filter_state_P*phi' + G*Qu*G' + Qnoise;

    % Save state for RTS
    rts_history_phi(:,:,i) = phi;
end

% smooth data backwards with RTS
[filter_state_rts, filter_state_rts_P] = kalman_rts( ...
        rts_history_filter_state_apriori, ...
        rts_history_filter_state_aposteriori, ...
        rts_history_filter_state_P_apriori, ...
        rts_history_filter_state_P_aposteriori, ...
        rts_history_phi);

fprintf('Estimated TDEE at the beginning of input data: %.0f kcal (+- %.1f)\n', filter_state_rts(2, 1), sqrt(filter_state_rts_P(2,2,1)));
fprintf('Estimated TDEE at the end of input data: %.0f kcal (+- %.1f)\n', filter_state_rts(2, end), sqrt(filter_state_rts_P(2,2,end)));
fprintf('Number of epochs with a mismatch between calories and body mass measurements: %i\n', questionable_measurements);

% TDEE estimation (forward filtered only, no RTS)
% figure(1);
% hold on;
% cov_tdee = reshape(rts_history_filter_state_P_aposteriori(2,2,:), [], 1)';
% std_tdee = sqrt(cov_tdee);
% plot(history_time_in_days, rts_history_filter_state_aposteriori(2, :) + std_tdee, 'r--');
% plot(history_time_in_days, rts_history_filter_state_aposteriori(2, :) - std_tdee, 'r--');
% plot(history_time_in_days, rts_history_filter_state_aposteriori(2, :), 'b');
% xlabel('Time (days)');
% ylabel('TDEE in kcal');
% title('TDEE estimation (forward)');
% legend('1 sigma band', 'TDEE');
% ylim([0 max(rts_history_filter_state_aposteriori(2, :) + std_tdee)]);
% hold off;

% figure;
% hold on;
% cov_tdee = reshape(filter_state_rts_P(2,2,:), [], 1)';
% std_tdee = sqrt(cov_tdee);
% xlabel('Time (days)');
% title('TDEE estimation (smooth)');
% plot(history_time_in_days, filter_state_rts(2, :) + std_tdee, 'k--');
% plot(history_time_in_days, filter_state_rts(2, :), 'b');
% plot(history_time_in_days, filter_state_rts(2, :) - std_tdee, 'k--');
% ylabel('TDEE (kcal)');
% ylim([0 max(filter_state_rts(2, :) + std_tdee)*1.2]);
% legend('1 sigma band TDEE', 'TDEE','Location', 'Best');
% hold off;

% TDEE estimation over time
% -------------------------
figure;
clf;
hold on;
grid on;
cov_tdee = reshape(filter_state_rts_P(2,2,:), [], 1)';
std_tdee = sqrt(cov_tdee);
xlabel('Time (days)');
title(sprintf('TDEE estimation: %.0f kcal', filter_state_rts(2, end)));
yyaxis left
% smooth_window_size = ceil(length(filter_state_rts(1, :))*0.2);
% weight_filtered_kg = movmean(filter_state_rts(1, :), smooth_window_size);
% plot(history_time_in_days, weight_filtered_kg, 'b');
plot(history_time_in_days, filter_state_rts(1, :), 'b');
ylabel('Mass (kg)');
yyaxis right
plot(history_time_in_days, filter_state_rts(2, :), 'r-');
plot(history_time_in_days, filter_state_rts(2, :) + std_tdee, 'k--');
plot(history_time_in_days, filter_state_rts(2, :) - std_tdee, 'k--');
ylabel('TDEE (kcal)');
ylim([0 max(filter_state_rts(2, :) + std_tdee)*1.2]);
legend('Body Mass (kg)', 'TDEE (kcal)', '1 sigma band TDEE', 'Location', 'Best');
hold off;

% Body mass over time
% -------------------
figure;
clf;
hold on;
grid on;
plot(history_time_in_days, filter_state_rts(1, :), 'b');
plot(history_time_in_days, history_weight_raw, 'x', 'Color', [0.8 0.8 0.8]);
p = polyfit(history_time_in_days, filter_state_rts(1, :), 1);
slope_kg_per_day = p(1); % weight loss trend kg/day
y_fit = polyval(p, history_time_in_days);
plot(history_time_in_days, y_fit, 'g--');
xlabel('Time (days)');
ylabel('Body mass (kg)');
title(sprintf('Body mass over time. Overall trend: %.2f kg/day', slope_kg_per_day));
legend('Body Mass (filtered)', 'Body Mass (measured)', 'Overall trend', 'Location', 'Best');
hold off;
fprintf('Average weight change per day: %.0f g, %.1f kg/week\n', slope_kg_per_day*1000, slope_kg_per_day*7);

% % Energy balance over time
% % ------------------------
% figure;
% kcal_per_day = history_kcals ./ history_dt_in_days;
% useful_epochs_kcal_per_day = kcal_per_day > median(kcal_per_day)*0.3;
% kcal_per_day = kcal_per_day(useful_epochs_kcal_per_day);
% kcal_balance = kcal_per_day - filter_state_rts(2, useful_epochs_kcal_per_day)';
% clf;
% hold on;
% grid on;
% smooth_window_size = ceil(length(filter_state_rts(1, :))*0.3);
% energy_surplus_filtered_kcal = movmean(kcal_balance, smooth_window_size);
% plot(history_time_in_days(useful_epochs_kcal_per_day), energy_surplus_filtered_kcal, 'r-');
% plot(history_time_in_days(useful_epochs_kcal_per_day), kcal_balance, 'x', 'Color', [0.8 0.8 0.8]);
% % Curve fitting: trend of energy balance
% p = polyfit(history_time_in_days(useful_epochs_kcal_per_day), kcal_balance, 1);
% y_fit = polyval(p, history_time_in_days(useful_epochs_kcal_per_day));
% plot(history_time_in_days(useful_epochs_kcal_per_day), y_fit, 'g--');
% 
% xlabel('Time (days)');
% ylabel('Energy balance (kcal)');
% legend('Energy balance (filtered)', 'Energy balance (data points)', sprintf('Trend %.0f kcal/day', p(1)));
% title(sprintf('Energy balance. Avg: %.0f kcal/day', mean(energy_surplus_filtered_kcal)));
% hold off;

% Energy Balance (alternative way via change in body mass)
% --------------------------------------------------------
time_period_days = ceil(history_time_in_days(end)); % total number of days
time_days = 1:time_period_days; % get integer days
% one body mass in kg per day (interpolate):
body_mass_interp = interp1(history_time_in_days, filter_state_rts(1, :), time_days, 'spline', 'extrap' );
% change in body mass per day:
delta_body_mass_per_day = [0, body_mass_interp(2:end)-body_mass_interp(1:(end-1))];
% calculate energy surplus based on body mass change:
energy_surplus_kcal = kcal_per_kg_fat * delta_body_mass_per_day;

smooth_window_size = ceil(length(filter_state_rts(1, :))*0.35);
energy_surplus_filtered_kcal = movmean(energy_surplus_kcal, smooth_window_size);
figure;
clf;
hold on;
grid on;
plot(time_days, energy_surplus_filtered_kcal, 'r-');
plot(time_days, energy_surplus_kcal, 'x', 'Color', [0.8 0.8 0.8]);
% Curve fitting: trend of energy balance
p = polyfit(time_days, energy_surplus_kcal, 1);
y_fit = polyval(p, time_days);
plot(time_days, y_fit, 'g--');
xlabel('Time (days)');
ylabel('Energy balance (kcal)');
legend('Energy balance (filtered)', 'Energy balance (data points)', sprintf('Trend %.0f kcal/day', p(1)));
title(sprintf('Energy balance. Avg: %.0f kcal/day', mean(energy_surplus_filtered_kcal)));
hold off;

end

