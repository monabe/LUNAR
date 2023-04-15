%%  Script to upload/measure/calculate knee flexion from a simple/stationary knee flexion test.
% Algorithms by Song et al. 2022.
% REFER to the paper:  S. Y. Song, Y. Pei, and E. T. Hsiao-Wecksler, “Estimating Relative Angles Using Two Inertial Measurement Units Without Magnetometers,” IEEE Sens. J., 2021.
% Date: April 2023
% Group 13 Capstone

%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intializing 

% Clear all variables, close all figures
clear all; clc; close all;

% Define data type (filtered vs. raw) for processing and plotting
data_type = 'raw'; 

% Set sampling frequency of IMU in Hz
sampling_freq = 40;

%% Add function folders to path

% CHANGE FILE PATH FOR YOUR COMPUTER
addpath('C:\Users\monab\Documents\MATLAB\LUNAR_KneeFlexion\Function Files from Song et al. 2022')  

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.SelectedVariableNames = ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
shank = readtable("C:\Users\monab\Documents\MATLAB\LUNAR_KneeFlexion\shank.TXT", opts);

% Convert to output type
shank = table2array(shank);

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.SelectedVariableNames = ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
thigh = readtable("C:\Users\monab\Documents\MATLAB\LUNAR_KneeFlexion\thigh.TXT", opts);

% Convert to output type
thigh = table2array(thigh);

% Clear temporary variables
clear opts

%% Saving accel and gyro data
% LUNAR sensors are in this order: Gyro-Accel-Mag

%when initially time-syncing these, make these range from 0 to the full
%length of the array
x0 = 250; %start of data of interest
x1 = 450; %end of data of interest

accel_1_x_raw = thigh(x0:x1,4); % in g's
accel_1_y_raw = thigh(x0:x1,5);
accel_1_z_raw = thigh(x0:x1,6);

accel_2_x_raw = shank(x0:x1,4);
accel_2_y_raw = shank(x0:x1,5);
accel_2_z_raw = shank(x0:x1,6);

gyro_1_x_raw = thigh(x0:x1,1); % in (deg/s)
gyro_1_y_raw = thigh(x0:x1,2);
gyro_1_z_raw = thigh(x0:x1,3);

gyro_2_x_raw = shank(x0:x1,1);
gyro_2_y_raw = shank(x0:x1,2);
gyro_2_z_raw = shank(x0:x1,3);

%% Plotting RAW accel and gyro data for thigh and shank

figure(1)
subplot(2,2,1)
plot(accel_1_x_raw)
hold on;
plot(accel_1_y_raw)
plot(accel_1_z_raw)
hold off 
title("Raw Acceration (Thigh)")
subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
legend("AccelX", "AccelY", "AccelZ")
xlabel("Samples")
ylabel("m/s^2")

subplot(2,2,2)
plot(gyro_1_x_raw)
hold on;
plot(gyro_1_y_raw)
plot(gyro_1_z_raw)
hold off 
title('Raw Angular Rate (Thigh)')
subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
legend("GyroX", "GyroY", "GyroZ")
xlabel("Samples")
ylabel("Degrees/second")

subplot(2,2,3)
plot(accel_2_x_raw)
hold on;
plot(accel_2_y_raw)
plot(accel_2_z_raw)
hold off 
title("Raw Acceration (Shank)")
subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
legend("AccelX", "AccelY", "AccelZ")
xlabel("Samples")
ylabel("m/s^2")

subplot(2,2,4)
plot(gyro_2_x_raw)
hold on;
plot(gyro_2_y_raw)
plot(gyro_2_z_raw)
hold off 
title('Raw Angular Rate (Shank)')
subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
legend("GyroX", "GyroY", "GyroZ")
xlabel("Samples")
ylabel("Degrees/second")

%% Filter Raw Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filtering Accelerometer Data (low pass filter)
    % Apply low pass filter for accelerometer data
    lpf_filter_order = 4;
    lpf_cutoff_freq = 4; % in Hz                     

    [bLP, aLP] = butter(lpf_filter_order,lpf_cutoff_freq/(sampling_freq/2));
    
    accel_1_x_filt = filtfilt(bLP, aLP, accel_1_x_raw);
    accel_1_y_filt = filtfilt(bLP, aLP, accel_1_y_raw);
    accel_1_z_filt = filtfilt(bLP, aLP, accel_1_z_raw);

    accel_2_x_filt = filtfilt(bLP, aLP, accel_2_x_raw);
    accel_2_y_filt = filtfilt(bLP, aLP, accel_2_y_raw);
    accel_2_z_filt = filtfilt(bLP, aLP, accel_2_z_raw);
    
 %% Filtering Gyroscope Data (high pass filter)
    hpf_filter_order = 4;
    hpf_cutoff_freq = 0.07;

    type = 'high';

    [bHP, aHP] = butter(hpf_filter_order, hpf_cutoff_freq/(sampling_freq/2), type); 

    % Apply filter. Use filtfilt to minimize any lag introduced by the filter 
    gyro_1_x_filt = filtfilt(bLP, aLP, gyro_1_x_raw);
    gyro_1_y_filt = filtfilt(bLP, aLP, gyro_1_y_raw);
    gyro_1_z_filt = filtfilt(bLP, aLP, gyro_1_z_raw);

    gyro_2_x_filt = filtfilt(bLP, aLP, gyro_2_x_raw);
    gyro_2_y_filt = filtfilt(bLP, aLP, gyro_2_y_raw);
    gyro_2_z_filt = filtfilt(bLP, aLP, gyro_2_z_raw);

%% Plotting FILTERED accel and gyro data for thigh and shank

    figure(2)
    subplot(2,2,1)
    plot(accel_1_x_filt)
    hold on;
    plot(accel_1_y_filt)
    plot(accel_1_z_filt)
    hold off 
    title("Filtered Acceration (Thigh)")
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    legend("AccelX", "AccelY", "AccelZ")
    xlabel("Samples")
    ylabel("m/s^2")

    subplot(2,2,2)
    plot(gyro_1_x_filt)
    hold on;
    plot(gyro_1_y_filt)
    plot(gyro_1_z_filt)
    hold off 
    title('Filtered Angular Rate (Thigh)')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    legend("GyroX", "GyroY", "GyroZ")
    xlabel("Samples")
    ylabel("Degrees/second")

    subplot(2,2,3)
    plot(accel_2_x_filt)
    hold on;
    plot(accel_2_y_filt)
    plot(accel_2_z_filt)
    hold off 
    title("Filtered Acceration (Shank)")
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    legend("AccelX", "AccelY", "AccelZ")
    xlabel("Samples")
    ylabel("m/s^2")

    subplot(2,2,4)
    plot(gyro_2_x_filt)
    hold on;
    plot(gyro_2_y_filt)
    plot(gyro_2_z_filt)
    hold off 
    title('Filtered Angular Rate (Shank)')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    legend("GyroX", "GyroY", "GyroZ")
    xlabel("Samples")
    ylabel("Degrees/second")

%% Compute Angle from Raw Data using Various Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1 [Accelerometer Inclincation - Pitch] Knee Angle
    
    rotation_type = 'pitch';

    % offset needed since the encoder started 90deg from IMU 
    ac_offset = 0; 
    
    % AC angle of IMU 1 
    ac_angle_1 = compute_accel_inclination_angle(accel_1_x_filt, accel_1_y_filt, accel_1_z_filt, 'pitch', 'stationary', ac_offset);    

    % AC angle of IMU 2 
    ac_angle_2 = compute_accel_inclination_angle(accel_2_x_filt, accel_2_y_filt, accel_2_z_filt, 'pitch', 'moving', ac_offset);

    % AC angle of IMU 1 - IMU 2
    ac_angle = ac_angle_1 - ac_angle_2;
  
    figure(3)
    subplot(3,2,1)
    plot(ac_angle);
    title('[Accelerometer Inclincation - Pitch] Knee Angle')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min') 
    
%% 1.2 [Accelerometer Inclincation - Roll] Knee Angle

    rotation_type = 'roll';

    % offset needed since the encoder started 90deg from IMU 
    ac_offset = 0; 
    
    % AC angle of IMU 1 
    ac_angle_1 = compute_accel_inclination_angle(accel_1_x_filt, accel_1_y_filt, accel_1_z_filt, 'pitch', 'stationary', ac_offset);    

    % AC angle of IMU 2 
    ac_angle_2 = compute_accel_inclination_angle(accel_2_x_filt, accel_2_y_filt, accel_2_z_filt, 'pitch', 'moving', ac_offset);

    % AC angle of IMU 1 - IMU 2
    ac_angle = ac_angle_1 - ac_angle_2;
  
    subplot(3,2,2)
    plot(ac_angle);
    title('[Accelerometer Inclincation - Roll] Knee Angle')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min') 
 
    %% 2.1 [Complementary filter - Pitch] Knee Angle
    rotation_type = 'pitch';
    dt = 1/sampling_freq;
    ac_offset = 0;
    
     if strcmp(rotation_type, 'pitch') == 1 || strcmp(rotation_type, 'roll') == 1 
        cf_offset = 0;
        gamma = 0.11; % tuning parameter for CF.  
         
         % CF angle of IMU 1
        cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x_filt, accel_1_y_filt, accel_1_z_filt, gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type, 'stationary', ac_offset);

        % CF angle of IMU 2 
        cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x_filt, accel_2_y_filt, accel_2_z_filt, gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type, 'moving', ac_offset);

    elseif strcmp(rotation_type, 'yaw') == 1 
         gamma = 0.0000001; % set to small value to rely mostly on gyroscopic data for yaw rotations
         cf_offset = 90;
         
         % CF angle of IMU 1
         cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x_filt, accel_1_y_filt, accel_1_z_filt, gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type, 'stationary', ac_offset);

         % CF angle of IMU 2 
         cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x_filt, accel_2_y_filt, accel_2_z_filt, gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type, 'moving', ac_offset);
    end
    
    
    % CF angle of IMU 2 - IMU 1
    cf_angle_pitch = -(cf_angle_2 - cf_angle_1 + cf_offset);
    
    subplot(3,2,3)
    plot(cf_angle_pitch);
    title('[Complementary Filter - Pitch] Knee Angle')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min') 
    
    %% 2.2 [Complementary filter - Roll] Knee Angle
    rotation_type = 'roll';
    dt = 1/sampling_freq;
    ac_offset = 0;
    
     if strcmp(rotation_type, 'pitch') == 1 || strcmp(rotation_type, 'roll') == 1 
        cf_offset = 0;
        gamma = 0.11; % tuning parameter for CF.  
         
         % CF angle of IMU 1
        cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x_filt, accel_1_y_filt, accel_1_z_filt, gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type, 'stationary', ac_offset);

        % CF angle of IMU 2 
        cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x_filt, accel_2_y_filt, accel_2_z_filt, gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type, 'moving', ac_offset);

    elseif strcmp(rotation_type, 'yaw') == 1 
         gamma = 0.0000001; % set to small value to rely mostly on gyroscopic data for yaw rotations
         cf_offset = 90;
         
         % CF angle of IMU 1
         cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x_filt, accel_1_y_filt, accel_1_z_filt, gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type, 'stationary', ac_offset);

         % CF angle of IMU 2 
         cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x_filt, accel_2_y_filt, accel_2_z_filt, gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type, 'moving', ac_offset);
     end
    
    subplot(3,2,4)

    % CF angle of IMU 2 - IMU 1
    cf_angle_roll = -(cf_angle_2 - cf_angle_1 + cf_offset);
    
    plot(cf_angle_roll);
    title('[Complementary Filter - Roll] Knee Angle')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min') 
    
    %% 3.1 [Gyroscope Integration - Pitch] Knee Angle
      rotation_type = 'pitch';
      filter_type = 'filtered';
    dt = 1/sampling_freq;
     % offset needed since the encoder started 90deg from IMU 
    gi_offset = 1.2; 
    if strcmp(rotation_type, 'pitch') == 1 || strcmp(rotation_type, 'roll') == 1 
        % GI angle of IMU 1 
        gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type);
        % GI angle of IMU 2 
        gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type);

    elseif strcmp(rotation_type, 'yaw') == 1       
         if strcmp(filter_type, 'filtered') == 1 
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type);
         else
             % below is just for reference to demonstrate that if gyro data is not high-pass filtered, the GI angle drifts significantly
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type);
        end
    end
   
    
    % GI angle of IMU 2 - IMU 1
    gi_angle = gi_angle_2 - gi_angle_1 + gi_offset;
    
    subplot(3,2,5)
    plot(rad2deg(gi_angle))
    title('[Gyroscope Integration - Pitch] Knee Angle')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min') 
    
   %% 3.2 [Gyroscope Integration - Roll] Knee Angle
    rotation_type = 'roll';
    filter_type = 'filtered';
    
    dt = 1/sampling_freq;
     % offset needed since the encoder started 90deg from IMU 
    gi_offset = 1.2; 
    if strcmp(rotation_type, 'pitch') == 1 || strcmp(rotation_type, 'roll') == 1 
        % GI angle of IMU 1 
        gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type);
        % GI angle of IMU 2 
        gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type);

    elseif strcmp(rotation_type, 'yaw') == 1       
         if strcmp(filter_type, 'filtered') == 1 
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type);
         else
             % below is just for reference to demonstrate that if gyro data is not high-pass filtered, the GI angle drifts significantly
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_filt, gyro_1_y_filt, gyro_1_z_filt, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_filt, gyro_2_y_filt, gyro_2_z_filt, dt, rotation_type);
        end
    end
   
    % GI angle of IMU 2 - IMU 1
    gi_angle = gi_angle_2 - gi_angle_1 + gi_offset;
    
    subplot(3,2,6)
    plot(rad2deg(gi_angle))
    title('[Gyroscope Integration - Roll] Knee Angle')
    subtitle('{\it 1 Knee Flexion Cycle to 90^{o}, 13 sec, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min') 