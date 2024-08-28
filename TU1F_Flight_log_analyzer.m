% TU-1.f Flight Log Visualizer
% Author: Jaewook Shim
% Date: 2024-08-19

%% HOUSEKEEPING
clear all; close all; clc;
format longG

%% UNIVERSAL CONSTANTS

g0 = 9.80665; % m/s2, gravitational acceleration
R = 287.05; % J/kg-K, gas constant of air

%% READ PARSED LOG FILE
if isfile("Data.mat")
    load("Data.mat")
else
    Data1 = readtable("2024-08-11 10-30-01 copy.xlsx"); % First log segment during flight
    Data2 = readtable("2024-08-11 10-39-01 copy.xlsx"); % Second log segment during flight
    Data3 = readtable("2024-08-11 10-40-46 copy.xlsx"); % Third log segment during flight
    save("Data.mat", "Data1", "Data2", "Data3");
end

%% Extract Data from Data1

time_ms = Data1.Time_ms_; % ms, timestamp
time_s = time_ms * 0.001; % s, timestamp

% GPS measurement
lat_deg = Data1.Lat_deg_; % deg, geodetic latitude from GPS
lon_deg = Data1.Lon_deg_; % deg, longitude from GPS
alt_msl_m = Data1.AltMSL_m_; % m, MSL altitude from GPS
geoid_m = Data1.Geoid_m_; % m, geoid separation from GPS

% launch site location configuration
ls_lat_deg = Data1.Lat_Launch_Site_deg_(end); % deg, launch site latitude
ls_lon_deg = Data1.Lon_Launch_Site_deg_(end); % deg, launch site longitude
ls_alt_msl_m = Data1.AltMSL_Launch_Site_m_(end); % m, launch site MSL
ls_geoid_m = Data1.Geoid_Launch_Site_m_(end); % m, launch site geoid separation

% linear acceleration in body frame
accX_ms2_body = Data1.AccX_m_s2_; % m/s2, accX in body frame
accY_ms2_body = Data1.AccY_m_s2_; % m/s2, accY in body frame
accZ_ms2_body = Data1.AccZ_m_s2_; % m/s2, accZ in body frame

% angular rate in body frame
gyroX_rads_body = Data1.GyroX_rad_s_; % rad/s, gyroX in body frame
gyroY_rads_body = Data1.GyroY_rad_s_; % rad/s, gyroY in body frame
gyroZ_rads_body = Data1.GyroZ_rad_s_; % rad/s, gyroZ in body frame

% Attitude quaternion from ENU to body
q1 = Data1.qX;
q2 = Data1.qY;
q3 = Data1.qZ;
q4 = Data1.qW;

% Barometer measurement
p_baro_hPa = Data1.p_baro_hPa_; % hPa, ambient pressure from barometer
T_baro_degC = Data1.T_baro_degC_; % degC, ambient(more technically, PCB) temperature from barometer

% Launch Site Configuration for Navigation
p_ls_hPa = Data1.P_Launch_Site_hPa_(end); % hPa, launch site ambient pressure
T_ls_degC = Data1.T_Launch_Site_degC_(end); % degC, launch site ambient(more technically, PCB) temperature
T_ls_K = T_ls_degC + 273.15; % K

% Raw Barometer Altitude, without LPF
h_baro_data1 = (R * T_ls_K / g0) * log(p_ls_hPa ./ p_baro_hPa); % m

% ADC Data
batV = Data1.V_bat_V_; % V, battery voltage per cell
P_chamber_barg = Data1.p_chamber_barg_; % barg, motor chamber pressure

% Navigation Data
posE = Data1.posE_m_; % m, E position
posN = Data1.posN_m_; % m, N position
posU = Data1.posU_m_; % m, U position, i.e., Above Ground Level altitude
velE = Data1.velE_m_s_; % m/s, E velocity
velN = Data1.velN_m_s_; % m/s, N velocity
velU = Data1.velU_m_s_; % m/s, U velocity
max_alt_m = Data1.max_alt_m_; % m, maximum altitude recorded so far

%% Extract Data from Data 2
time_ms_data2 = Data2.Time_ms_;
time_s_data2 = time_ms_data2 * 0.001;

% GPS measurement
lat_deg_data2 = Data2.Lat_deg_; % deg, geodetic latitude from GPS
lon_deg_data2 = Data2.Lon_deg_; % deg, longitude from GPS
alt_msl_m_data2 = Data2.AltMSL_m_; % m, MSL altitude from GPS
geoid_m_data2 = Data2.Geoid_m_; % m, geoid separation from GPS

% linear acceleration in body frame
accX_ms2_body_data2 = Data2.AccX_m_s2_; % m/s2, accX in body frame
accY_ms2_body_data2 = Data2.AccY_m_s2_; % m/s2, accY in body frame
accZ_ms2_body_data2 = Data2.AccZ_m_s2_; % m/s2, accZ in body frame

% angular rate in body frame
gyroX_rads_body_data2 = Data2.GyroX_rad_s_; % rad/s, gyroX in body frame
gyroY_rads_body_data2 = Data2.GyroY_rad_s_; % rad/s, gyroY in body frame
gyroZ_rads_body_data2 = Data2.GyroZ_rad_s_; % rad/s, gyroZ in body frame

q1_data2 = Data2.qX;
q2_data2 = Data2.qY;
q3_data2 = Data2.qZ;
q4_data2 = Data2.qW;

p_baro_hPa_data2 = Data2.p_baro_hPa_; % hPa, ambient pressure from barometer
h_baro_data2 = (R * T_ls_K / g0) * log(p_ls_hPa ./ p_baro_hPa_data2); % m

%% Extract Data from Data 3
time_ms_data3 = Data3.Time_ms_;
time_s_data3 = time_ms_data3 * 0.001;

% GPS measurement
lat_deg_data3 = Data3.Lat_deg_; % deg, geodetic latitude from GPS
lon_deg_data3 = Data3.Lon_deg_; % deg, longitude from GPS
alt_msl_m_data3 = Data3.AltMSL_m_; % m, MSL altitude from GPS
geoid_m_data3 = Data3.Geoid_m_; % m, geoid separation from GPS

% linear acceleration in body frame
accX_ms2_body_data3 = Data3.AccX_m_s2_; % m/s2, accX in body frame
accY_ms2_body_data3 = Data3.AccY_m_s2_; % m/s2, accY in body frame
accZ_ms2_body_data3 = Data3.AccZ_m_s2_; % m/s2, accZ in body frame

% angular rate in body frame
gyroX_rads_body_data3 = Data3.GyroX_rad_s_; % rad/s, gyroX in body frame
gyroY_rads_body_data3 = Data3.GyroY_rad_s_; % rad/s, gyroY in body frame
gyroZ_rads_body_data3 = Data3.GyroZ_rad_s_; % rad/s, gyroZ in body frame

q1_data3 = Data3.qX;
q2_data3 = Data3.qY;
q3_data3 = Data3.qZ;
q4_data3 = Data3.qW;

p_baro_hPa_data3 = Data3.p_baro_hPa_; % hPa, ambient pressure from barometer
h_baro_data3 = (R * T_ls_K / g0) * log(p_ls_hPa ./ p_baro_hPa_data3); % m

%% Plot
figure;
hold on
plot(time_s, q1, 'DisplayName', 'q1')
plot(time_s, q2, 'DisplayName', 'q2')
plot(time_s, q3, 'DisplayName', 'q3')
plot(time_s, q4, 'DisplayName', 'q4')
xlabel("time [sec]")
ylabel("quaternion [-]")
legend show
grid on
xlim([520, 545])

% Construct new time stamp (concatenate three chuncks of data)
t_offset_data2 = 9 * 60; % sec, timestamp offset from data 1
t_offset_data3 = 60 + 45 + t_offset_data2; % sec, timestamp offset from data 1

time_s_data2_shifted = time_s_data2 + t_offset_data2;
time_s_data3_shifted = time_s_data3 + t_offset_data3;

% time_merged_s = [time_s; 
%                 time_s_data2 + t_offset_data2; 
%                 time_s_data3 + t_offset_data3]; % s, merged timestamp for data 1, 2, and 3

% Plot barometer altitude from data 1
figure;
hold on
plot(time_s, posU, 'DisplayName', 'LPF, Data1', 'LineWidth', 1.5);
plot(time_s, h_baro_data1, 'DisplayName', 'Raw, Data1', 'LineWidth', 1.5);

% Plot baromter altitude from data 2
plot(time_s_data2_shifted, h_baro_data2, 'DisplayName', 'Raw, Data2', 'LineWidth', 1.5);

% Plot barometer altitude from data 3
plot(time_s_data3_shifted, h_baro_data3, 'DisplayName', 'Raw, Data3', 'LineWidth', 1.5);

xlabel("time [sec]")
ylabel("altitude [m]")
grid on
title("Altitude from BMP280, Above Ground Level (== posU)")
legend show
xlim([520 650])

%% Animation
quat_offset = [0; 0; sind(-45); cosd(-45)];

quat_history_data1 = [q1, q2, q3, q4];
quat_history_data1 = transpose(quat_history_data1);
% fix orientation offset
for i=1:1:length(quat_history_data1)
    quat_history_data1(:, i) = quaternion_product(quat_history_data1(:, i), quat_offset);
end

quat_history_data2 = [q1_data2, q2_data2, q3_data2, q4_data2];
quat_history_data2 = transpose(quat_history_data2);
% fix orientation offset
for i=1:1:length(quat_history_data2)
    quat_history_data2(:, i) = quaternion_product(quat_history_data2(:, i), quat_offset);
end

quat_history_data3 = [q1_data3, q2_data3, q3_data3, q4_data3];
quat_history_data3 = transpose(quat_history_data3);
% fix orientation offset
for i=1:1:length(quat_history_data3)
    quat_history_data3(:, i) = quaternion_product(quat_history_data3(:, i), quat_offset);
end

% during flight
dummy_time = (time_s(end):0.01:time_s_data2_shifted(1))'; % missing time between data1 and data2
N = length(dummy_time);
dummy_quat = zeros(4, N);

flight_merged_time_s = [time_s; dummy_time; time_s_data2_shifted];
flight_merged_quat_history = [quat_history_data1, dummy_quat, quat_history_data2];
% animate_with_cylinder(flight_merged_time_s(42000:end), flight_merged_quat_history(:, 42000:end));

% % full version including recovery
% full_merged_time_s = [time_s; time_s_data2_shifted; time_s_data3_shifted];
% full_merged_quat_history = [quat_history_data1, quat_history_data2, quat_history_data3];
% animate_with_cylinder(full_merged_time_s(42000:end), full_merged_quat_history(:, 42000:end));

%% %%%%%%%%%%%%%%% Applying Kalman Filter (for data1) %%%%%%%%%%%%%%%

% LPF gyroscope data
Ts = 0.01; % sec
f_cutoff = 5; % Hz
tau = 1 / 2/pi / f_cutoff; % time constant
alpha = Ts / (tau + Ts);

N = length(time_s);

filtered_gyro_x = zeros(N, 1);
filtered_gyro_y = zeros(N, 1);
filtered_gyro_z = zeros(N, 1);
filtered_gyro_x(1) = gyroX_rads_body(1); % rad/s
filtered_gyro_y(1) = gyroY_rads_body(1); % rad/s
filtered_gyro_z(1) = gyroZ_rads_body(1); % rad/s

for i = 2:1:N
    filtered_gyro_x(i) = alpha * gyroX_rads_body(i) + (1 - alpha) * filtered_gyro_x(i-1);
    filtered_gyro_y(i) = alpha * gyroY_rads_body(i) + (1 - alpha) * filtered_gyro_y(i-1);
    filtered_gyro_z(i) = alpha * gyroZ_rads_body(i) + (1 - alpha) * filtered_gyro_z(i-1);
end

dt = diff(time_s);
angular_acc_x = diff(gyroX_rads_body) ./ dt;
angular_acc_y = diff(gyroY_rads_body) ./ dt;
angular_acc_z = diff(gyroZ_rads_body) ./ dt;
filtered_angular_acc_x = diff(filtered_gyro_x) ./ dt;
filtered_angular_acc_y = diff(filtered_gyro_y) ./ dt;
filtered_angular_acc_z = diff(filtered_gyro_z) ./ dt;

figure;
subplot(3, 1, 1);
hold on
plot(time_s, gyroX_rads_body, 'DisplayName', 'gyro_x')
plot(time_s, filtered_gyro_x, 'DisplayName', 'filtered-gyro_x')
legend show
title("Angular Rate: Original vs. Filtered")

subplot(3, 1, 2);
hold on
plot(time_s, gyroY_rads_body, 'DisplayName', 'gyro_y')
plot(time_s, filtered_gyro_y, 'DisplayName', 'filtered-gyro_y')
legend show

subplot(3, 1, 3);
hold on
plot(time_s, gyroZ_rads_body, 'DisplayName', 'gyro_z')
plot(time_s, filtered_gyro_z, 'DisplayName', 'filtered-gyro_z')
legend show

% Plot angular acceleration obtained by finite difference
figure; 
subplot(3, 1, 1);
hold on
plot(time_s(1:end-1), angular_acc_x, 'DisplayName', 'angular-acc_x')
plot(time_s(1:end-1), filtered_angular_acc_x, 'DisplayName', 'filtered-angular-acc_x')
legend show
title("Angular Acceleration: Original vs. Filtered")

subplot(3, 1, 2);
hold on
plot(time_s(1:end-1), angular_acc_y, 'DisplayName', 'angular-acc_-')
plot(time_s(1:end-1), filtered_angular_acc_y, 'DisplayName', 'filtered-angular-acc_y')
legend show

subplot(3, 1, 3);
hold on
plot(time_s(1:end-1), angular_acc_z, 'DisplayName', 'angular-acc_z')
plot(time_s(1:end-1), filtered_angular_acc_z, 'DisplayName', 'filtered-angular-acc_z')
legend show

r_IMU = [0.236287; 0; 0]; % m

angular_acc_x = [0; angular_acc_x];
angular_acc_y = [0; angular_acc_y];
angular_acc_z = [0; angular_acc_z];

acc_ms2_imu = [accX_ms2_body'; accY_ms2_body'; accZ_ms2_body'];
angular_rate = [gyroX_rads_body'; gyroY_rads_body'; gyroZ_rads_body'];
angular_acc = [angular_acc_x'; angular_acc_y'; angular_acc_z'];

acc_ms2_cg = zeros(3, N);

for i = 1:1:N
    acc_ms2_cg(:, i) = acc_ms2_imu(:, i) - cross(angular_acc(:, i), r_IMU) - cross(angular_rate(:, i), cross(angular_rate(:, i), r_IMU));
end

figure;

subplot(3, 1, 1)
title("ACC felt by IMU and ACC of CG")
hold on
plot(time_s, accX_ms2_body, 'DisplayName', 'IMU ACC X')
plot(time_s, acc_ms2_cg(1, :), 'DisplayName', 'CG ACC X')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 2)
title("ACC felt by IMU and ACC of CG")
hold on
plot(time_s, accY_ms2_body, 'DisplayName', 'IMU ACC Y')
plot(time_s, acc_ms2_cg(2, :), 'DisplayName', 'CG ACC Y')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 3)
title("ACC felt by IMU and ACC of CG")
hold on
plot(time_s, accZ_ms2_body, 'DisplayName', 'IMU ACC Z')
plot(time_s, acc_ms2_cg(3, :), 'DisplayName', 'CG ACC Z')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

% Convert Acc from body to ENU frame
acc_cg_in_ENU = zeros(3, N);
for i=1:1:N
    q = quat_history_data1(:, i); % quaternion
    DCM_from_ENU_to_body = quat_2_dcm(q);
    DCM_from_body_to_ENU = transpose(DCM_from_ENU_to_body);
    acc_cg_in_ENU(:, i) = DCM_from_body_to_ENU * acc_ms2_cg(:, i);
end

figure;

subplot(3, 1, 1)
title("ACC in ENU frame")
hold on
plot(time_s, acc_cg_in_ENU(1, :), 'DisplayName', 'CG ACC E')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 2)
hold on
plot(time_s, acc_cg_in_ENU(2, :), 'DisplayName', 'CG ACC N')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 3)
hold on
plot(time_s, acc_cg_in_ENU(3, :), 'DisplayName', 'CG ACC U')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

% Apply Kalman Filter

u = acc_cg_in_ENU(:, 42272:end);
z = h_baro_data1(42272:end); 
dt = dt(42272:end);

N = length(u);

x = zeros(6, N);

x(:, 1) = [0; 0; 0; 0; 0; 0]; % initialize state vector

P = zeros(6, 6, N);
P(:, :, 1) = 10 * eye(6); % initialize state covariance

for i = 2:1:N
    inc_t = dt(i-1);
    x_prev = x(:, i-1);
    u_prev = u(:, i-1);
    P_prev = P(:, :, i-1);
    z_curr = z(i-1);

    [x_curr, P_curr] = KalmanFilter(inc_t, x_prev, u_prev, P_prev, z_curr);
    x(:, i) = x_curr;
    P(:, :, i) = P_curr;
end

figure
hold on
plot(time_s(42272:end), x(1, :), 'DisplayName', 'pos E (m)')
plot(time_s(42272:end), x(2, :), 'DisplayName', 'pos N (m)')
plot(time_s(42272:end), x(3, :), 'DisplayName', 'pos U (m)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
title("Estimated Position in ENU frame")
legend show

figure
hold on
plot(time_s(42272:end), x(4, :), 'DisplayName', 'vel E (m/s)')
plot(time_s(42272:end), x(5, :), 'DisplayName', 'vel N (m/s)')
plot(time_s(42272:end), x(6, :), 'DisplayName', 'vel U (m/s)')
xlabel("time [sec]")
ylabel("velocity [m/s]")
grid on
title("Estimated Velocity in ENU frame")
legend show

figure
plot(time_s(42272:end), vecnorm(x(4:6, :)), 'DisplayName', 'speed (m/s)')
xlabel("time [sec]")
ylabel("speed [m/s]")
grid on
title("Estimated Speed")
legend show

figure
plot3(x(1, :), x(2, :), x(3, :), 'DisplayName', 'Flight Path', 'LineWidth', 1.5)
xlabel("pos E [m]")
ylabel("pos N [m]")
zlabel("pos U [m]")
grid on
title("Flight Path")
legend show
axis equal

%% Compare with GPS lat lon result
ls_pos_ECEF = lla2ecef([ls_lat_deg, ls_lon_deg, ls_alt_msl_m + ls_geoid_m])'; % m

sLam = sind(ls_lon_deg);
cLam = cosd(ls_lon_deg);
sPhi = sind(ls_lat_deg);
cPhi = cosd(ls_lat_deg);

DCM_ECEF_to_ENU = [-sLam, cLam, 0;
                   -sPhi * cLam, -sPhi * sLam, cPhi;
                   cPhi * cLam, cPhi * sLam, sPhi];

N = length(lat_deg);
r_ECEF = zeros(3, N); % m, ECEF position of the vehicle
r_rel_ECEF = zeros(3, N); % m, relative position of the vehicle wrt origin, in ECEF
pos_ENU_from_gps = zeros(3, N); % m, ENU position directly obtained from GPS lla

for i=1:1:N
    r_ECEF(:, i) = lla2ecef([lat_deg(i), lon_deg(i), alt_msl_m(i) + geoid_m(i)])'; % m
    r_rel_ECEF(:, i) = r_ECEF(:, i) - ls_pos_ECEF;
    pos_ENU_from_gps(:, i) = DCM_ECEF_to_ENU * r_rel_ECEF(:, i);
end

figure
hold on
plot(time_s(42272:end), pos_ENU_from_gps(1, 42272:end), 'DisplayName', 'pos E (m)')
plot(time_s(42272:end), pos_ENU_from_gps(2, 42272:end), 'DisplayName', 'pos N (m)')
plot(time_s(42272:end), pos_ENU_from_gps(3, 42272:end), 'DisplayName', 'pos U (m)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
title("GPS measured position in ENU frame")
legend show

% Compare GPS measurement and Kalman Filtered result
figure;

subplot(3, 1, 1)
title("Comparison of GPS measurement and KF results")
hold on
plot(time_s(42272:end), x(1, :), 'DisplayName', 'posE (KF)')
plot(time_s(42272:end), pos_ENU_from_gps(1, 42272:end), 'DisplayName', 'posE (GPS)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
legend show
legend('Location','northeastoutside')

subplot(3, 1, 2)
hold on
plot(time_s(42272:end), x(2, :), 'DisplayName', 'posN (KF)')
plot(time_s(42272:end), pos_ENU_from_gps(2, 42272:end), 'DisplayName', 'posN (GPS)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
legend show
legend('Location','northeastoutside')

subplot(3, 1, 3)
hold on
plot(time_s(42272:end), x(3, :), 'DisplayName', 'posU (KF)')
plot(time_s(42272:end), pos_ENU_from_gps(3, 42272:end), 'DisplayName', 'posU (GPS)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
legend show
legend('Location','northeastoutside')

%% %%%%%%%%%%%%%%% Applying Kalman Filter (for data2) %%%%%%%%%%%%%%%

% LPF gyroscope data
Ts = 0.01; % sec
f_cutoff = 5; % Hz
tau = 1 / 2/pi / f_cutoff; % time constant
alpha = Ts / (tau + Ts);

N = length(time_s_data2_shifted);

filtered_gyro_x_data2 = zeros(N, 1);
filtered_gyro_y_data2 = zeros(N, 1);
filtered_gyro_z_data2 = zeros(N, 1);
filtered_gyro_x_data2(1) = gyroX_rads_body_data2(1); % rad/s
filtered_gyro_y_data2(1) = gyroY_rads_body_data2(1); % rad/s
filtered_gyro_z_data2(1) = gyroZ_rads_body_data2(1); % rad/s

for i = 2:1:N
    filtered_gyro_x_data2(i) = alpha * gyroX_rads_body_data2(i) + (1 - alpha) * filtered_gyro_x_data2(i-1);
    filtered_gyro_y_data2(i) = alpha * gyroY_rads_body_data2(i) + (1 - alpha) * filtered_gyro_y_data2(i-1);
    filtered_gyro_z_data2(i) = alpha * gyroZ_rads_body_data2(i) + (1 - alpha) * filtered_gyro_z_data2(i-1);
end

dt_data2 = diff(time_s_data2_shifted);
angular_acc_x_data2 = diff(gyroX_rads_body_data2) ./ dt_data2;
angular_acc_y_data2 = diff(gyroY_rads_body_data2) ./ dt_data2;
angular_acc_z_data2 = diff(gyroZ_rads_body_data2) ./ dt_data2;
filtered_angular_acc_x_data2 = diff(filtered_gyro_x_data2) ./ dt_data2;
filtered_angular_acc_y_data2 = diff(filtered_gyro_y_data2) ./ dt_data2;
filtered_angular_acc_z_data2 = diff(filtered_gyro_z_data2) ./ dt_data2;

figure;
subplot(3, 1, 1);
hold on
plot(time_s_data2_shifted, gyroX_rads_body_data2, 'DisplayName', 'gyro_x')
plot(time_s_data2_shifted, filtered_gyro_x_data2, 'DisplayName', 'filtered-gyro_x')
legend show
title("Angular Rate: Original vs. Filtered (data2)")

subplot(3, 1, 2);
hold on
plot(time_s_data2_shifted, gyroY_rads_body_data2, 'DisplayName', 'gyro_y')
plot(time_s_data2_shifted, filtered_gyro_y_data2, 'DisplayName', 'filtered-gyro_y')
legend show

subplot(3, 1, 3);
hold on
plot(time_s_data2_shifted, gyroZ_rads_body_data2, 'DisplayName', 'gyro_z')
plot(time_s_data2_shifted, filtered_gyro_z_data2, 'DisplayName', 'filtered-gyro_z')
legend show

% Plot angular acceleration obtained by finite difference
figure; 
subplot(3, 1, 1);
hold on
plot(time_s_data2_shifted(1:end-1), angular_acc_x_data2, 'DisplayName', 'angular-acc_x')
plot(time_s_data2_shifted(1:end-1), filtered_angular_acc_x_data2, 'DisplayName', 'filtered-angular-acc_x')
legend show
title("Angular Acceleration: Original vs. Filtered (data2)")

subplot(3, 1, 2);
hold on
plot(time_s_data2_shifted(1:end-1), angular_acc_y_data2, 'DisplayName', 'angular-acc_y')
plot(time_s_data2_shifted(1:end-1), filtered_angular_acc_y_data2, 'DisplayName', 'filtered-angular-acc_y')
legend show

subplot(3, 1, 3);
hold on
plot(time_s_data2_shifted(1:end-1), angular_acc_z_data2, 'DisplayName', 'angular-acc_z')
plot(time_s_data2_shifted(1:end-1), filtered_angular_acc_z_data2, 'DisplayName', 'filtered-angular-acc_z')
legend show

r_IMU = [0.236287; 0; 0]; % m

angular_acc_x_data2 = [0; angular_acc_x_data2];
angular_acc_y_data2 = [0; angular_acc_y_data2];
angular_acc_z_data2 = [0; angular_acc_z_data2];

acc_ms2_imu_data2 = [accX_ms2_body_data2'; accY_ms2_body_data2'; accZ_ms2_body_data2'];
angular_rate_data2 = [gyroX_rads_body_data2'; gyroY_rads_body_data2'; gyroZ_rads_body_data2'];
angular_acc_data2 = [angular_acc_x_data2'; angular_acc_y_data2'; angular_acc_z_data2'];

acc_ms2_cg_data2 = zeros(3, N);

for i = 1:1:N
    acc_ms2_cg_data2(:, i) = acc_ms2_imu_data2(:, i) - cross(angular_acc_data2(:, i), r_IMU) - cross(angular_rate_data2(:, i), cross(angular_rate_data2(:, i), r_IMU));
end

figure;

subplot(3, 1, 1)
title("ACC felt by IMU and ACC of CG (data2)")
hold on
plot(time_s_data2_shifted, accX_ms2_body_data2, 'DisplayName', 'IMU ACC X')
plot(time_s_data2_shifted, acc_ms2_cg_data2(1, :), 'DisplayName', 'CG ACC X')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 2)
hold on
plot(time_s_data2_shifted, accY_ms2_body_data2, 'DisplayName', 'IMU ACC Y')
plot(time_s_data2_shifted, acc_ms2_cg_data2(2, :), 'DisplayName', 'CG ACC Y')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 3)
hold on
plot(time_s_data2_shifted, accZ_ms2_body_data2, 'DisplayName', 'IMU ACC Z')
plot(time_s_data2_shifted, acc_ms2_cg_data2(3, :), 'DisplayName', 'CG ACC Z')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

% Convert Acc from body to ENU frame
acc_cg_in_ENU_data2 = zeros(3, N);
for i=1:1:N
    q = quat_history_data2(:, i); % quaternion
    DCM_from_ENU_to_body = quat_2_dcm(q);
    DCM_from_body_to_ENU = transpose(DCM_from_ENU_to_body);
    acc_cg_in_ENU_data2(:, i) = DCM_from_body_to_ENU * acc_ms2_cg_data2(:, i);
end

figure;

subplot(3, 1, 1)
title("ACC in ENU frame (data2)")
hold on
plot(time_s_data2_shifted, acc_cg_in_ENU_data2(1, :), 'DisplayName', 'CG ACC E')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 2)
hold on
plot(time_s_data2_shifted, acc_cg_in_ENU_data2(2, :), 'DisplayName', 'CG ACC N')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

subplot(3, 1, 3)
hold on
plot(time_s_data2_shifted, acc_cg_in_ENU_data2(3, :), 'DisplayName', 'CG ACC U')
xlabel("time [sec]")
ylabel("acceleration [m/s2]")
grid on
legend show

% Apply Kalman Filter

u_data2 = acc_cg_in_ENU_data2;
z_data2 = h_baro_data2; 

N = length(z_data2);

x_data2 = zeros(6, N);

% initialization: x(1:3) -> pos U is replaced with the h_baro_data2(1), others are last estimates from data1
% x(4:6) -> velocities are obtained by assuming 11.038 sec free fall during
% FC shutdown, calculating from last estimated velocity from data 1

vE_ini = x(4, end);
vN_ini = x(5, end);
vU_ini = x(6, end) - g0 * 11.038;

x_data2(:, 1) = [164.418586986295; 64.7628863717792; 360.703905305781; vE_ini; vN_ini; vU_ini]; % initialize state vector 

P_data2 = zeros(6, 6, N);
P_data2(:, :, 1) = 10 * eye(6); % initialize state covariance

for i = 2:1:N
    inc_t = dt_data2(i-1);
    x_prev = x_data2(:, i-1);
    u_prev = u_data2(:, i-1);
    P_prev = P_data2(:, :, i-1);
    z_curr = z_data2(i-1);

    [x_curr, P_curr] = KalmanFilter(inc_t, x_prev, u_prev, P_prev, z_curr);
    x_data2(:, i) = x_curr;
    P_data2(:, :, i) = P_curr;
end

figure
hold on
plot(time_s_data2_shifted, x_data2(1, :), 'DisplayName', 'pos E (m)')
plot(time_s_data2_shifted, x_data2(2, :), 'DisplayName', 'pos N (m)')
plot(time_s_data2_shifted, x_data2(3, :), 'DisplayName', 'pos U (m)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
title("Estimated Position in ENU frame (data 2)")
legend show

figure
hold on
plot(time_s_data2_shifted, x_data2(4, :), 'DisplayName', 'vel E (m/s)')
plot(time_s_data2_shifted, x_data2(5, :), 'DisplayName', 'vel N (m/s)')
plot(time_s_data2_shifted, x_data2(6, :), 'DisplayName', 'vel U (m/s)')
xlabel("time [sec]")
ylabel("velocity [m/s]")
grid on
title("Estimated Velocity in ENU frame (data 2)")
legend show

figure
plot(time_s_data2_shifted, vecnorm(x_data2(4:6, :)), 'DisplayName', 'speed (m/s)')
xlabel("time [sec]")
ylabel("speed [m/s]")
grid on
title("Estimated Speed (data2)")
legend show

figure
plot3(x_data2(1, :), x_data2(2, :), x_data2(3, :), 'DisplayName', 'Flight Path', 'LineWidth', 1.5)
xlabel("pos E [m]")
ylabel("pos N [m]")
zlabel("pos U [m]")
grid on
title("Flight Path (data2)")
legend show
axis equal

%% Compare with GPS lat lon result
N = length(lat_deg_data2);
r_ECEF_data2 = zeros(3, N); % m, ECEF position of the vehicle
r_rel_ECEF_data2 = zeros(3, N); % m, relative position of the vehicle wrt origin, in ECEF
pos_ENU_from_gps_data2 = zeros(3, N); % m, ENU position directly obtained from GPS lla

for i=1:1:N
    r_ECEF_data2(:, i) = lla2ecef([lat_deg_data2(i), lon_deg_data2(i), alt_msl_m_data2(i) + geoid_m_data2(i)])'; % m
    r_rel_ECEF_data2(:, i) = r_ECEF_data2(:, i) - ls_pos_ECEF;
    pos_ENU_from_gps_data2(:, i) = DCM_ECEF_to_ENU * r_rel_ECEF_data2(:, i);
end

figure
hold on
plot(time_s_data2_shifted, pos_ENU_from_gps_data2(1, :), 'DisplayName', 'pos E (m)')
plot(time_s_data2_shifted, pos_ENU_from_gps_data2(2, :), 'DisplayName', 'pos N (m)')
plot(time_s_data2_shifted, pos_ENU_from_gps_data2(3, :), 'DisplayName', 'pos U (m)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
title("GPS measured position in ENU frame (data2)")
legend show

% Compare GPS measurement and Kalman Filtered result
figure;

subplot(3, 1, 1)
title("Comparison of GPS measurement and KF results (data2)")
hold on
plot(time_s_data2_shifted, x_data2(1, :), 'DisplayName', 'posE (KF)')
plot(time_s_data2_shifted, pos_ENU_from_gps_data2(1, :), 'DisplayName', 'posE (GPS)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
legend show
legend('Location','northeastoutside')

subplot(3, 1, 2)
hold on
plot(time_s_data2_shifted, x_data2(2, :), 'DisplayName', 'posN (KF)')
plot(time_s_data2_shifted, pos_ENU_from_gps_data2(2, :), 'DisplayName', 'posN (GPS)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
legend show
legend('Location','northeastoutside')

subplot(3, 1, 3)
hold on
plot(time_s_data2_shifted, x_data2(3, :), 'DisplayName', 'posU (KF)')
plot(time_s_data2_shifted, pos_ENU_from_gps_data2(3, :), 'DisplayName', 'posU (GPS)')
xlabel("time [sec]")
ylabel("position [m]")
grid on
legend show
legend('Location','northeastoutside')


%% Combine data 1 and data2

figure
hold on
plot(time_s(42272:end), x(1, :), 'DisplayName', 'pos E (m)', 'Color', 'R')
plot(time_s(42272:end), x(2, :), 'DisplayName', 'pos N (m)', 'Color', 'G')
plot(time_s(42272:end), x(3, :), 'DisplayName', 'pos U (m)', 'Color', 'B')

plot(time_s_data2_shifted, x_data2(1, :), 'DisplayName', 'pos E (m)', 'Color', 'R')
plot(time_s_data2_shifted, x_data2(2, :), 'DisplayName', 'pos N (m)', 'Color', 'G')
plot(time_s_data2_shifted, x_data2(3, :), 'DisplayName', 'pos U (m)', 'Color', 'B')
xlabel("time [sec]")
ylabel("position [m]")
grid on
title("Estimated Position in ENU frame (Combined)")
legend show

figure
hold on
plot(time_s(42272:end), x(4, :), 'DisplayName', 'vel E (m/s)', 'Color', 'R')
plot(time_s(42272:end), x(5, :), 'DisplayName', 'vel N (m/s)', 'Color', 'G')
plot(time_s(42272:end), x(6, :), 'DisplayName', 'vel U (m/s)', 'Color', 'B')

plot(time_s_data2_shifted, x_data2(4, :), 'DisplayName', 'vel E (m/s)', 'Color', 'R')
plot(time_s_data2_shifted, x_data2(5, :), 'DisplayName', 'vel N (m/s)', 'Color', 'G')
plot(time_s_data2_shifted, x_data2(6, :), 'DisplayName', 'vel U (m/s)', 'Color', 'B')
xlabel("time [sec]")
ylabel("velocity [m/s]")
grid on
title("Estimated Velocity in ENU frame (Combined)")
legend show

figure
hold on
plot(time_s(42272:end), vecnorm(x(4:6, :)), 'DisplayName', 'speed (m/s)')
plot(time_s_data2_shifted, vecnorm(x_data2(4:6, :)), 'DisplayName', 'speed (m/s)')
xlabel("time [sec]")
ylabel("speed [m/s]")
grid on
title("Estimated Speed (Combined)")
legend show

figure
plot3(x_data2(1, :), x_data2(2, :), x_data2(3, :), 'DisplayName', 'Flight Path', 'LineWidth', 1.5)
hold on
plot3(x(1, :), x(2, :), x(3, :), 'DisplayName', 'Flight Path', 'LineWidth', 1.5)

xlabel("pos E [m]")
ylabel("pos N [m]")
zlabel("pos U [m]")
grid on
title("Flight Path (Combined)")
legend show
axis equal

%% Calculate Euler Angle
RPY = quat2eul([q4, q1, q2, q3], 'ZYX'); % rad
roll = RPY(:, 3) * 180/pi;
pitch = RPY(:, 2) * 180/pi;
yaw = RPY(:, 1) * 180/pi;

figure;
subplot(3, 1, 1)
title("Euler Angle during ascent")
hold on
plot(time_s, roll, 'DisplayName', 'Roll')
xlabel("time [sec]")
ylabel("angle [deg]")
grid on
legend show
legend('Location','northeastoutside')
xlim([528, 540])

subplot(3, 1, 2)
hold on
plot(time_s, pitch, 'DisplayName', 'Pitch')
xlabel("time [sec]")
ylabel("angle [deg]")
grid on
legend show
legend('Location','northeastoutside')
xlim([528, 540])

subplot(3, 1, 3)
hold on
plot(time_s, yaw, 'DisplayName', 'Yaw')
xlabel("time [sec]")
ylabel("angle [deg]")
grid on
legend show
legend('Location','northeastoutside')
xlim([528, 540])


%% Flight Trajectory Animation
% animate_trajectory_ENU(time_s(42272:end), x(1:3, :), quat_history_data1(:, 42272:end))
% animate_combined(flight_merged_time_s(42272:end), x(1:3, :), flight_merged_quat_history(:, 42272:end), [1, 1, 1]) % similar to viewpoint of the recorded video)
animate_combined(flight_merged_time_s(42272:end), x(1:3, :), flight_merged_quat_history(:, 42272:end), [-0.1, -1, -0.3]) % similar to viewpoint of the recorded video)