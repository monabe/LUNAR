%%  Script to upload/measure/calculate knee flexion from a simple/stationary knee flexion test.
% Algorithms by Song et al. 2022, who have shared open-source code in their
% paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9888780 

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
shank04101626 = readtable("C:\Users\monab\Documents\MATLAB\studioApril10_sensors_3kneeflexions_thighaligned\shank_04101626.TXT", opts);

%% Convert to output type
shank = table2array(shank04101626);

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
thigh04101629 = readtable("C:\Users\monab\Documents\MATLAB\studioApril10_sensors_3kneeflexions_thighaligned\thigh_04101629.TXT", opts);

%% Convert to output type
thigh = table2array(thigh04101629);

%% Clear temporary variables
clear opts

%% Saving accel and gyro data
% LUNAR sensors are in this order: Gyro-Accel-Mag

x0 = 250;

x1 = 450;

% x1 = length(shank); %change this to be thigh or shank, depending on which is shortest

accel_1_x = thigh(x0:x1,4).*9.81; %converting from g's to m/s^2
accel_1_y = thigh(x0:x1,5).*9.81;
accel_1_z = thigh(x0:x1,6).*9.81;

accel_2_x = shank(x0:x1,4).*9.81;
accel_2_y = shank(x0:x1,5).*9.81;
accel_2_z = shank(x0:x1,6).*9.81;

gyro_1_x_raw = thigh(x0:x1,1);
gyro_1_y_raw = thigh(x0:x1,2);
gyro_1_z_raw = thigh(x0:x1,3);

gyro_2_x_raw = shank(x0:x1,1);
gyro_2_y_raw = shank(x0:x1,2);
gyro_2_z_raw = shank(x0:x1,3);

%% Plotting raw accel and gyro data for thigh and shank

figure(1)
subplot(2,2,1)
plot(accel_1_x)
hold on;
plot(accel_1_y)
plot(accel_1_z)
hold off 
title("Raw Acceration (Thigh)")
subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
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
subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
legend("GyroX", "GyroY", "GyroZ")
xlabel("Samples")
ylabel("Degrees/second")

subplot(2,2,3)
plot(accel_2_x)
hold on;
plot(accel_2_y)
plot(accel_2_z)
hold off 
title("Raw Acceration (Shank)")
subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
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
subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
legend("GyroX", "GyroY", "GyroZ")
xlabel("Samples")
ylabel("Degrees/second")

%% Sampling frequency

sampling_freq = 40;

%% 1.1 [Accelerometer Inclincation - Pitch] Knee Angle
    % Note that for yaw, the AC angles are zeros since it cannot be
    % computed for yaw. 
    
    rotation_type = 'pitch';

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
    subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
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
    subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
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
    subtitle('{\it Knee Flexion to 90 deg, LUNAR Sensors}')
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
         cf_angle_1 = compute_complementary_filter_angle(gamma, accel_1_x, accel_1_y, accel_1_z, gyro_1_x_raw, gyro_1_y_raw, gyro_1_z_raw, dt, rotation_type, 'stationary', ac_offset);

         % CF angle of IMU 2 
         cf_angle_2 = compute_complementary_filter_angle(gamma, accel_2_x, accel_2_y, accel_2_z, gyro_2_x_raw, gyro_2_y_raw, gyro_2_z_raw, dt, rotation_type, 'moving', ac_offset);
     end
    
    subplot(3,2,4)

    % CF angle of IMU 2 - IMU 1
    cf_angle_roll = cf_angle_2 - cf_angle_1 + cf_offset;
    
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