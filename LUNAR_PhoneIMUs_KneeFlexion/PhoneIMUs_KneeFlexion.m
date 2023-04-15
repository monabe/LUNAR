%%  Script to upload/measure/calculate knee flexion from a simple/stationary knee flexion test.
% Algorithms by Song et al. 2022, who have shared open-source code in their
% paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9888780 
% Date: April 2023
% Group 13 Capstone


%% Add function folders to path

% CHANGE FILE PATH FOR YOUR COMPUTER
addpath('C:\Users\monab\Documents\MATLAB\LUNAR_PhoneIMUs_KneeFlexion\Function Files from Song et al. 2022') 

%% Import the accel_shank data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Xms2", "Yms2", "Zms2"];
opts.SelectedVariableNames = ["Xms2", "Yms2", "Zms2"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
accelshank = readtable("C:\Users\monab\Documents\MATLAB\LUNAR_PhoneIMUs_KneeFlexion\accel_shank.csv", opts);

accelshank = table2array(accelshank);

%% Import the accel_thigh data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "AccelerationXms2", "AccelerationYms2", "AccelerationZms2"];
opts.SelectedVariableNames = ["AccelerationXms2", "AccelerationYms2", "AccelerationZms2"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
accelthigh = readtable("C:\Users\monab\Documents\MATLAB\LUNAR_PhoneIMUs_KneeFlexion\accel_thigh.csv", opts);

accelthigh = table2array(accelthigh);

%% Import the gyro_shank data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Xrads", "Yrads", "Zrads"];
opts.SelectedVariableNames = ["Xrads", "Yrads", "Zrads"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
gyroshank = readtable("C:\Users\monab\Documents\MATLAB\LUNAR_PhoneIMUs_KneeFlexion\gyro_shank.csv", opts);

gyroshank = table2array(gyroshank);

%% Import the gyro_thigh data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "GyroscopeXrads", "GyroscopeYrads", "GyroscopeZrads"];
opts.SelectedVariableNames = ["GyroscopeXrads", "GyroscopeYrads", "GyroscopeZrads"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
gyrothigh = readtable("C:\Users\monab\Documents\MATLAB\LUNAR_PhoneIMUs_KneeFlexion\gyro_thigh.csv", opts);

gyrothigh = table2array(gyrothigh);

%% Saving accel and gyro data
accel_1_x = accelthigh(2:380,1); %2nd column due to heading on 1st row
accel_1_y = accelthigh(2:380,2);
accel_1_z = accelthigh(2:380,3);

accel_2_x = accelshank(2:380,1);
accel_2_y = accelshank(2:380,2);
accel_2_z = accelshank(2:380,3);

gyro_1_x_raw = gyrothigh(2:380,1);
gyro_1_y_raw = gyrothigh(2:380,2);
gyro_1_z_raw = gyrothigh(2:380,3);

gyro_2_x_raw = gyroshank(2:380,1);
gyro_2_y_raw = gyroshank(2:380,2);
gyro_2_z_raw = gyroshank(2:380,3);


%% Plotting raw accel and gyro data for thigh and shank

figure(1)
subplot(2,2,1)
plot(accel_1_x)
hold on;
plot(accel_1_y)
plot(accel_1_z)
hold off 
title("Raw Acceration (Shank)")
subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
legend("AccelX", "AccelY", "AccelZ")
xlabel("Samples")
ylabel("m/s^2")

subplot(2,2,2)
plot(gyro_1_x_raw)
hold on;
plot(gyro_1_y_raw)
plot(gyro_1_z_raw)
hold off 
title('Raw Angular Rate (Shank)')
subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
legend("GyroX", "GyroY", "GyroZ")
xlabel("Samples")
ylabel("Degrees/second")

subplot(2,2,3)
plot(accel_2_x)
hold on;
plot(accel_2_y)
plot(accel_2_z)
hold off 
title("Raw Acceration (Thigh)")
subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
legend("AccelX", "AccelY", "AccelZ")
xlabel("Samples")
ylabel("m/s^2")

subplot(2,2,4)
plot(gyro_2_x_raw)
hold on;
plot(gyro_2_y_raw)
plot(gyro_2_z_raw)
hold off 
title('Raw Angular Rate (Thigh)')
subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
legend("GyroX", "GyroY", "GyroZ")
xlabel("Samples")
ylabel("Degrees/second")

%% Sampling frequency

sampling_freq = 100; %for this test and these phones

%% 1.1 [Accelerometer Inclincation - Pitch] Knee Angle
    % Note that for yaw, the AC angles are zeros since it cannot be
    % computed for yaw. 
    
    % offset needed since the encoder started 90deg from IMU 
    ac_offset = 0; 
    
    % AC angle of IMU 1 
    ac_angle_1 = compute_accel_inclination_angle(accel_1_x, accel_1_y, accel_1_z, 'pitch', 'stationary', ac_offset);    

    % AC angle of IMU 2 
    ac_angle_2 = compute_accel_inclination_angle(accel_2_x, accel_2_y, accel_2_z, 'pitch', 'moving', ac_offset);

    % AC angle of IMU 2 - IMU 1
    ac_angle = ac_angle_2 - ac_angle_1;
  
    figure(2)
    subplot(3,2,1)
    plot(ac_angle);
    title('[Accelerometer Inclincation - Pitch] Knee Angle')
    subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min')    
        
%% 1.2 [Accelerometer Inclincation - Roll] Knee Angle
    % Note that for yaw, the AC angles are zeros since it cannot be
    % computed for yaw. 
    
    rotation_type = 'roll';

    % offset needed since the encoder started 90deg from IMU 
    ac_offset = 0; 
    
    % AC angle of IMU 1 
    ac_angle_1 = compute_accel_inclination_angle(accel_1_x, accel_1_y, accel_1_z, 'pitch', 'stationary', ac_offset);    

    % AC angle of IMU 2 
    ac_angle_2 = compute_accel_inclination_angle(accel_2_x, accel_2_y, accel_2_z, 'pitch', 'moving', ac_offset);

    % AC angle of IMU 2 - IMU 1
    ac_angle = ac_angle_2 - ac_angle_1;
  
    subplot(3,2,2)
    plot(ac_angle);
    title('[Accelerometer Inclincation - Roll] Knee Angle')
    subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
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
        cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x, accel_1_y, accel_1_z, gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type, 'stationary', ac_offset);

        % CF angle of IMU 2 
        cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x, accel_2_y, accel_2_z, gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type, 'moving', ac_offset);

    elseif strcmp(rotation_type, 'yaw') == 1 
         gamma = 0.0000001; % set to small value to rely mostly on gyroscopic data for yaw rotations
         cf_offset = 90;
         
         % CF angle of IMU 1
         cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x, accel_1_y, accel_1_z, gyro_1_x, gyro_1_y, gyro_1_z, dt, rotation_type, 'stationary', ac_offset);

         % CF angle of IMU 2 
         cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x, accel_2_y, accel_2_z, gyro_2_x, gyro_2_y, gyro_2_z, dt, rotation_type, 'moving', ac_offset);
    end
    
    
    % CF angle of IMU 2 - IMU 1
    cf_angle_pitch = cf_angle_2 - cf_angle_1 + cf_offset;
    
    subplot(3,2,3)
    plot(cf_angle_pitch);
    title('[Complementary Filter - Pitch] Knee Angle')
    subtitle('{\it Knee Flexion to 90 deg, Phone IMUs}')
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
        cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x, accel_1_y, accel_1_z, gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type, 'stationary', ac_offset);

        % CF angle of IMU 2 
        cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x, accel_2_y, accel_2_z, gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type, 'moving', ac_offset);

    elseif strcmp(rotation_type, 'yaw') == 1 
         gamma = 0.0000001; % set to small value to rely mostly on gyroscopic data for yaw rotations
         cf_offset = 90;
         
         % CF angle of IMU 1
         cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x, accel_1_y, accel_1_z, gyro_1_x, gyro_1_y, gyro_1_z, dt, rotation_type, 'stationary', ac_offset);

         % CF angle of IMU 2 
         cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x, accel_2_y, accel_2_z, gyro_2_x, gyro_2_y, gyro_2_z, dt, rotation_type, 'moving', ac_offset);
    end
    
    
    % CF angle of IMU 2 - IMU 1
    cf_angle_roll = cf_angle_2 - cf_angle_1 + cf_offset;
    
    subplot(3,2,4)
    plot(cf_angle_roll);
    title('[Complementary Filter - Roll] Knee Angle')
    subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
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
        gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type);
        % GI angle of IMU 2 
        gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type);

    elseif strcmp(rotation_type, 'yaw') == 1       
         if strcmp(filter_type, 'filtered') == 1 
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type);
         else
             % below is just for reference to demonstrate that if gyro data is not high-pass filtered, the GI angle drifts significantly
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type);
        end
    end
   
    
    % GI angle of IMU 2 - IMU 1
    gi_angle = gi_angle_2 - gi_angle_1 + gi_offset;
    
    subplot(3,2,5)
    plot(rad2deg(gi_angle))
    title('[Gyroscope Integration - Pitch] Knee Angle')
    subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")

    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min')
    
    %% 3.2 [Gyroscope Integration - Roll] Knee Angle
    rotation_type = 'roll';
    dt = 1/sampling_freq;
     % offset needed since the encoder started 90deg from IMU 
    gi_offset = 1.2; 
    if strcmp(rotation_type, 'pitch') == 1 || strcmp(rotation_type, 'roll') == 1 
        % GI angle of IMU 1 
        gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type);
        % GI angle of IMU 2 
        gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type);

    elseif strcmp(rotation_type, 'yaw') == 1       
         if strcmp(filter_type, 'filtered') == 1 
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type);
         else
             % below is just for reference to demonstrate that if gyro data is not high-pass filtered, the GI angle drifts significantly
            gi_angle_1 = compute_gyro_integration_angle(gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type);
            gi_angle_2 = compute_gyro_integration_angle(gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type);
        end
    end
   
    
    % GI angle of IMU 2 - IMU 1
    gi_angle = gi_angle_2 - gi_angle_1 + gi_offset;
    
    subplot(3,2,6)
    plot(rad2deg(gi_angle))
    title('[Gyroscope Integration - Roll] Knee Angle')
    subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
    xlabel("Samples")
    ylabel("Degrees")
    
    ymax = 90;
    ymin = 0;
    yline(ymax,'--','Max')
    yline(ymin,'--','Min')
