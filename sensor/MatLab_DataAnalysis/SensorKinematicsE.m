

function SensorKinematicsE(TaskTitle, SubjectNumInt, TrialNumInt, Activity, plot_mode)

%Activity: 1 Block, 2 Drink (Soda Can), 3 Water (Water Pouring)

%% Clearing Workspace
%clc; clear;close all;
close all;

%% Import Libraries
% addpath('Quaternions');

%% Initialize Reference/Output Varibles

% Do Not Delete this section. This will make sure the csv output is as
% correct as possible should the script or test fail.
errorVelocity = 0;
errorScript = 0;
salimu = 0;
mtimu = 0;
ImaxVel = 0;
corrvaluex= 0; RMSE = 0; peakdetThres = 0;
jerkimu_dim =0;
jerkimu_dim_log =0;
ImaxVel = 0; ImeanVelocity =0;
ImaxAccel = 0;  ImeanAccel = 0;   stdAccel = 0;
TimeFromOnsetToLastPeak = 0;   TimeFromOnsetToFirstPeak =0;
TimeFromOnsetToHighestPeak = 0;   TimeFromOnsetToLowestPeak = 0;

%% File Name Parameters:

SubjectNum = int2str(SubjectNumInt);

TrialNum = int2str(TrialNumInt);

%% IMU Data Entry (Paste Data or Read From File)
% Load IMU File
data = csvread(strcat(TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv'), 1,0); % IMU file name will have two '_'
% data = csvread(strcat(TaskTitle,'_',SubjectNum,'_', TrialNum,'.txt'), 1,0); % IMU file name will have two '_'

%If necessary trim the data for better results:
%Example:
%data = data(11:(end-30),:);  For all Trials:  Eliminates the first 10 frames and the last 30 frames.

%Example: Trim Select Subject(s) and Trial(s)
% if(   ((SubjectNumInt==5) && (TrialNumInt ==13)) || ( ((SubjectNumInt==15) && (TrialNumInt ==12)))   )
%     data = data(11:(end-70),:);
% end


try %try/catch will run the catch block when an error occurs in the try block.
    
    
    
    %% IMU Data pre-processing on Import
    samplePeriod = 1/100;
    
    g = 9.81;
    accX = data(:,1)/g;
    accY = data(:,2)/g;
    accZ = data(:,3)/g+1;
    
    gyrX = data(:,4);
    gyrY = data(:,5);
    gyrZ = data(:,6);
    
    
    time = 0:samplePeriod:(length(data(:,1))-1)*samplePeriod;
    startTime = 0;
    stopTime = time(end);
    
    % -------------------------------------------------------------------------
    % Manually frame data
    
    indexSel = find(sign(time-startTime)+1, 1) : find(sign(time-stopTime)+1, 1);
    time = time(indexSel);
    gyrX = gyrX(indexSel, :);
    gyrY = gyrY(indexSel, :);
    gyrZ = gyrZ(indexSel, :);
    accX = accX(indexSel, :);
    accY = accY(indexSel, :);
    accZ = accZ(indexSel, :);
    
    % -------------------------------------------------------------------------
    % Detect stationary periods
    
    % Compute accelerometer magnitude
    acc_mag = sqrt(accX.*accX + accY.*accY + accZ.*accZ);
    
    % HP filter accelerometer data
    filtCutOff = 0.001;
    [b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'high');
    acc_magFilt = filtfilt(b, a, acc_mag);
    
    % Compute absolute value
    acc_magFilt = abs(acc_magFilt);
    
    % LP filter accelerometer data
    filtCutOff = 2;
    [b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'low');
    acc_magFilt = filtfilt(b, a, acc_magFilt);
    
    
    %    Threshold detection
    for h = 0.025:0.002:0.1
        stationary = acc_magFilt < h;
        if (   (mean(stationary((end-20):end))) < 1  )  || ( (mean(stationary(1:10)))<1  )
            
        else
            break;
        end
    end
    
    stationary(1:5) = 1;
    w = find(stationary<1,1);
    r = find(stationary<1,1,'last');
    mt_imu_with_stationary = r-w;
    
    if Activity <=3
        stationary(w:r) = acc_magFilt(w:r) < 0.018;
    else 
        stationary(w:r) = acc_magFilt(w:r) < 0.04;
    end
    
    if plot_mode == 1
        % -------------------------------------------------------------------------
        % Plot data raw sensor data and stationary periods
        f1 = figure(1);
        f1.Visible = 'off';
        f1.Name = 'Raw IMU Sensor Data';
        
        ax(1) = subplot(2,1,1);
        hold on;
        plot(time, gyrX, 'r');
        plot(time, gyrY, 'g');
        plot(time, gyrZ, 'b');
        title('Gyroscope');
        xlabel('Time (s)');
        ylabel('Angular velocity (^\circ/s)');
        legend('X', 'Y', 'Z');
        hold off;
        ax(2) = subplot(2,1,2);
        hold on;
        plot(time, accX, 'r');
        plot(time, accY, 'g');
        plot(time, accZ, 'b');
        plot(time, acc_magFilt, ':k');
        plot(time, stationary, 'k', 'LineWidth', 2);
        title('Accelerometer');
        xlabel('Time (s)');
        ylabel('Acceleration (g)');
        legend('X', 'Y', 'Z', 'Filtered', 'Stationary');
        hold off;
        linkaxes(ax,'x');
    end
    % -------------------------------------------------------------------------
    % Compute orientation
    
    quat = zeros(length(time), 4);
    %AHRSalgorithm = AHRS('SamplePeriod', 1/100, 'Kp', 1, 'KpInit', 1);
    AHRSalgorithm = AHRS('SamplePeriod', 1/100, 'Kp', 5, 'KpInit', 0.1);
    
    % Initial convergence
    initPeriod = 0.2;
    indexSel = 1 : find(sign(time-(time(1)+initPeriod))+1, 1);
    for i = 1:200
        AHRSalgorithm.UpdateIMU([0 0 0], [mean(accX(indexSel)) mean(accY(indexSel)) mean(accZ(indexSel))]);
    end
    
    % For all data
    for t = 1:length(time)
        if(stationary(t))
            AHRSalgorithm.Kp = 0.5;
        else
            AHRSalgorithm.Kp = 0;
        end
        AHRSalgorithm.UpdateIMU(deg2rad([gyrX(t) gyrY(t) gyrZ(t)]), [accX(t) accY(t) accZ(t)]);
        quat(t,:) = AHRSalgorithm.Quaternion;
    end
    
    % -------------------------------------------------------------------------
    % Compute translational accelerations
    
    % Rotate body accelerations to Earth frame
    acc = quaternRotate([accX accY accZ], quaternConj(quat));
    
    % % Remove gravity from measurements
    % acc = acc - [zeros(length(time), 2) ones(length(time), 1)];     % unnecessary due to velocity integral drift compensation
    
    % Convert acceleration measurements to m/s/s
    acc = acc * 9.81;
    
    % -------------------------------------------------------------------------
    % Compute translational velocities
    
    acc(:,3) = acc(:,3) - 9.81;
    
    
    % Integrate acceleration to yield velocity
    vel = zeros(size(acc));
    for t = 2:length(vel)
        vel(t,:) = vel(t-1,:) + acc(t,:) * samplePeriod;
        if(stationary(t) == 1)
            vel(t,:) = [0 0 0];     % force zero velocity when foot stationary
        end
    end
    
    % Compute integral drift during non-stationary periods
    velDrift = zeros(size(vel));
    stationaryStart = find([0; diff(stationary)] == -1);
    stationaryEnd = find([0; diff(stationary)] == 1);
    if length(stationaryStart)<length(stationaryEnd)
        stationaryStart = [1; stationaryStart(1:end)];
    end
    
    
    for i = 1:numel(stationaryEnd)
        driftRate = vel(stationaryEnd(i)-1, :) / (stationaryEnd(i) - stationaryStart(i));
        enum = 1:(stationaryEnd(i) - stationaryStart(i));
        drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
        velDrift(stationaryStart(i):stationaryEnd(i)-1, :) = drift;
    end
    
    % Remove integral drift
    vel = vel - velDrift;
    
    
    if plot_mode == 1
        % Plot translational velocity
        %figure('Position', [9 39 900 300], 'NumberTitle', 'off', 'Name', 'Velocity');
        f3 = figure(3);
        f3.Visible = 'off';
        f3.Name = 'IMU Velocity and Position';
        subplot (2,1,1); hold on;
        plot(time, acc(:,1), 'r');
        plot(time, acc(:,2), 'g');
        plot(time, acc(:,3), 'b');
        title('Acceleration');
        xlabel('Time (s)');
        ylabel('Acceleration (m/s/s)');
        legend('X', 'Y', 'Z');
        hold off;
        
        subplot (2,1,2)
        hold on;
        plot(time, vel(:,1), 'r');hold on;
        plot(time, vel(:,2), 'g');hold on;
        plot(time, vel(:,3), 'b');hold on;
        title('Velocity');
        xlabel('Time (s)');
        ylabel('Velocity (m/s)');
        legend('X', 'Y', 'Z');
        hold off;
        
    end
    
    % -------------------------------------------------------------------------
    % Compute translational position
    
    % Integrate velocity to yield position
    pos = zeros(size(vel));
    for t = 2:length(pos)
        pos(t,:) = pos(t-1,:) + vel(t,:) * samplePeriod;    % integrate velocity to yield position
    end
    
    % -------------------------------------------------------------------------
    % Plot 3D  trajectory
    
    % % Remove stationary periods from data to plot
    posPlot = pos;
    quatPlot = quat;
    
    % Extend final sample to delay end of animation
    extraTime = 20;
    onesVector = ones(extraTime*(1/samplePeriod), 1);
    posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
    quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];
    
    %% Assigning Filtered IMU Data to new variable.
    
    IMU_vx = vel(:,1);
    IMU_vy = vel(:,2);
    IMU_vz = vel(:,3);
    
    IMU_ax = acc(:,1);
    IMU_ay = acc(:,2);
    IMU_az = acc(:,3);
    
    
    %% Data import, Domenico's Filter, and length correction happen above.
    
    
    IMUA =  sqrt(IMU_ax.^2 + IMU_ay.^2 + IMU_az.^2);
    
    if plot_mode == 1
        f69 = figure(69);
        f69.Visible = 'off';
        f69.Name = 'Sensor Data: Accelerations';
        plot(IMUA)
        title('Acceleration');
        legend('IMU');
    end
    
    
    IMUV = (sqrt(IMU_vx.^2 + IMU_vy.^2 +IMU_vz.^2));
    
    
    if IMUV(10) > 0.1 || IMUV(end) > 0.1
        errorVelocity = 1;
        strcat('Errort: Possible Bad Velocity Signals Calculated for Subject: ', int2str(SubjectNumInt), ' Trial: ', int2str(TrialNumInt)) %Prints the error message to the command window.
        
    end
    
    
    peakdetThres = 0.1;
    
    peakdetThres4 = 0;
    peakdetThres3 = 0;
    
    for i = 1:70
        %Detect Peaks in Time Series (wmc)
        %Vicon
        
        %     [pks, locs] = findpeaks(ViconV,1,'MinPeakHeight',peakdetThres);
        %     maxtab_ViconV = [locs pks];
        
        %Detect Peaks in Time Series (wmc)
        %IMU
        [maxtab_IMUV, mintab_IMUV] = peakdet(IMUV, peakdetThres);
        %     [pks, locs] = findpeaks(IMUV,1,'MinPeakHeight',peakdetThres);
        %     maxtab_IMUV = [locs pks];
        
        if length(maxtab_IMUV)>=3 && i<70
            
            if (length(maxtab_IMUV)==4)
                peakdetThres4 = peakdetThres;
            end
            
            if (length(maxtab_IMUV)==3)
                peakdetThres3 = peakdetThres;
            end
            
            peakdetThres = peakdetThres*1.05;
        else
            if i > 1
                peakdetThres = peakdetThres/1.05;
            end
            break;
        end
    end
    
    
    if peakdetThres4 >0
        peakdetThres=peakdetThres4;
    elseif peakdetThres3 >0
        peakdetThres=peakdetThres3;
    end
    
    [maxtab_IMUV, mintab_IMUV] = peakdet(IMUV, peakdetThres);
    
    
    
    if plot_mode == 1
        
        f68 = figure(68);
        set(f68,'Position',[320 385 560 420]); % CHANGED
        f68.Visible = 'off';
        f68.Name = 'Sensor Data: Velocities';
        hold on;
        
    end
    
    %IMU
    temp = strcat(strcat('Subject: ', SubjectNum),strcat(' Trial: ',TrialNum)); % CHANGED
    if plot_mode == 1
        plot (IMUV, 'b'); hold on
        plot(maxtab_IMUV(:,1), maxtab_IMUV(:,2), 'k*');hold on; title(strcat(temp,' Velocity')); % CHANGED
    end
    
    
    MinDuration = 50;
    
    
    %Onset
    if any(IMUV(maxtab_IMUV(1,1):-1:1)<maxtab_IMUV(2,2)*0.01)
        for onset_IMUV = maxtab_IMUV(1,1):-1:1
            
            if onset_IMUV - MinDuration > 1
                minZeroVelocityFrom =  onset_IMUV - MinDuration;
            else
                minZeroVelocityFrom = 1;
            end
            
            
            if (IMUV(onset_IMUV)<maxtab_IMUV(2,2)*0.01 ) && (mean(IMUV(minZeroVelocityFrom:onset_IMUV))<0.1)
                break
            end
        end
    else
        for onset_IMUV = maxtab_IMUV(1,1):-1:1
            
            if onset_IMUV - MinDuration > 1
                minZeroVelocityFrom =  onset_IMUV - MinDuration;
            else
                minZeroVelocityFrom = 1;
            end
            
            
            if (IMUV(onset_IMUV)< maxtab_IMUV(2,2)*0.05) && (mean(IMUV(minZeroVelocityFrom:onset_IMUV))<0.1)
                break
            end
        end
    end
    
    
    MinDuration = 75;
    
    
    
   
    %Offset of IMU
    if any(IMUV(maxtab_IMUV(length(maxtab_IMUV),1):1:end)<maxtab_IMUV(length(maxtab_IMUV),2)*0.01)
        for offset_IMUV = (maxtab_IMUV(length(maxtab_IMUV),1):1:length(IMUV)-1)
            
            
            
            if offset_IMUV + MinDuration < length(IMUV)
                minZeroVelocityUntil =  offset_IMUV + MinDuration;
            else
                minZeroVelocityUntil = length(IMUV);
            end
            
            
            
            if (IMUV(offset_IMUV)< maxtab_IMUV(length(maxtab_IMUV),2)*0.01) || ((mean(IMUV(offset_IMUV:minZeroVelocityUntil))<0.02) && (IMUV(offset_IMUV)< 0.005))
                break
            end
        end
    elseif any(IMUV(maxtab_IMUV(length(maxtab_IMUV),1):1:size(IMUV)<maxtab_IMUV(length(maxtab_IMUV),2)*0.05))
        for offset_IMUV = (maxtab_IMUV(length(maxtab_IMUV),1):1:length(IMUV)-1)
            
            if offset_IMUV + MinDuration < length(IMUV)
                minZeroVelocityUntil =  offset_IMUV + MinDuration;
            else
                minZeroVelocityUntil = length(IMUV);
            end
            
            
            
            if (IMUV(offset_IMUV)< maxtab_IMUV(length(maxtab_IMUV),2)*0.05) || ((mean(IMUV(offset_IMUV:minZeroVelocityUntil))<0.02) && (IMUV(offset_IMUV)< 0.005))
                break
            end
        end
    else
        [s, loc_offset_IMUV] = min(IMUV(maxtab_IMUV(length(maxtab_IMUV),1):1:length(IMUV)-1));
        offset_IMUV = loc_offset_IMUV + maxtab_IMUV(length(maxtab_IMUV),1) -1;
    end
    
    if plot_mode == 1
        plot(onset_IMUV, IMUV(onset_IMUV), 'r*');
        plot(offset_IMUV, IMUV(offset_IMUV), 'r*');
        legend('IMU', 'IMU Peak' , 'IMU Onset', 'IMU Offset');
    end
    
    
    %% Trim Time Series
    %Plot Overall
    trimI = IMUV(onset_IMUV:offset_IMUV);
    
    if plot_mode == 1
        f70 = figure(70);
         set(f70,'Position',[881 385 560 420]);
        f70.Name = 'Trimmed Resultant Data';
        f70.Visible = 'off';
        hold on; plot(trimI,'b')
        title(strcat(temp,' Resultant Velocity')); % CHANGED
        legend('IMU');
    end
    
    
    %Correlation and RMSE
    %     corrx = corrcoef(trimV,trimI);
    %     corrvaluex = corrx (1,2);
    %
    %     RMSE = sqrt(mean((trimV-trimI).^2));
    
    
    
    %Trim all axes for IMU
    trimIx = IMU_vx(onset_IMUV:offset_IMUV);
    trimIy = IMU_vy(onset_IMUV:offset_IMUV);
    trimIz = IMU_vz(onset_IMUV:offset_IMUV);
    
    
    %----------------------------Kinematic Variables-------------------
    %Movement Time
    mtimu = (offset_IMUV-onset_IMUV); %value depends on sampling rate
    
    %Movement Smoothness
    salimu = func_SAL(trimIx, trimIy, trimIz);
    [jerkimu_dim,jerkimu_dim_log] = func_jerk_dim(trimIx, trimIy, trimIz);
    
    
    %Max Velocity
    ImaxVel = max(IMUV);
    
    %Max Acceleration
    ImaxAccel = max(IMUA);
    
    
    %Mean Acceleration
    ImeanAccel = mean(IMUA(onset_IMUV:offset_IMUV));
    
    
    %Mean Velocity
    ImeanVelocity = mean(IMUV(onset_IMUV:offset_IMUV));
    
    
    %Variability of Mean Accel
    
    stdAccel = std2(IMUA(onset_IMUV:offset_IMUV));
    
    
    
    
    if Activity == 1
        
        try
            
            peakLoc = maxtab_IMUV(:,1);
            peakValues = maxtab_IMUV(:,2);
            
            temp = maxk(peakValues,3);
            
            
            %Name of variables are named in order of highest peak.
            peak1 = peakLoc(find(peakValues == temp(1), 1 ,'first'));
            peak2 = peakLoc(find(peakValues == temp(2), 1 ,'first'));
            peak3 = peakLoc(find(peakValues == temp(3), 1 ,'first'));
            
            
            %Time From Onset to Peaks (for Velocity)
            
            TimeFromOnsetToLastPeak = max([peak1 peak2 peak3])-onset_IMUV;
            TimeFromOnsetToFirstPeak = min([peak1 peak2 peak3])-onset_IMUV;
            
            
            TimeFromOnsetToHighestPeak =  peakLoc(find(peakValues == max(temp), 1 ,'first'))-onset_IMUV;
            TimeFromOnsetToLowestPeak =  peakLoc(find(peakValues == min(temp), 1 ,'first'))-onset_IMUV;
            
            
        catch ME
            
            strcat('Error: 3 Peaks Not Detected.  Subject: ', int2str(SubjectNumInt), ' Trial: ', int2str(TrialNumInt)) %Prints the error message to the command window.
            
        end
        
    elseif (Activity == 2) || (Activity == 3)
        
        
        try
            
            peakLoc = maxtab_IMUV(:,1);
            peakValues = maxtab_IMUV(:,2);
            
            temp = maxk(peakValues,4);
            
            
            %Name of variables are named in order of highest peak.
            peak1 = peakLoc(find(peakValues == temp(1), 1 ,'first'));
            peak2 = peakLoc(find(peakValues == temp(2), 1 ,'first'));
            peak3 = peakLoc(find(peakValues == temp(3), 1 ,'first'));
            peak4 = peakLoc(find(peakValues == temp(4), 1 ,'first'));
            
            
            %Time From Onset to Peaks (for Velocity)
            
            TimeFromOnsetToLastPeak = max([peak1 peak2 peak3 peak4])-onset_IMUV;
            TimeFromOnsetToFirstPeak = min([peak1 peak2 peak3 peak4])-onset_IMUV;
            
            
            TimeFromOnsetToHighestPeak =  peakLoc(find(peakValues == max(temp), 1 ,'first'))-onset_IMUV;
            TimeFromOnsetToLowestPeak =  peakLoc(find(peakValues == min(temp), 1 ,'first'))-onset_IMUV;
            
            
        catch ME
            
            strcat('Error: 4 Peaks Not Detected. 3 Peaks will be attempted.  Subject: ', int2str(SubjectNumInt), ' Trial: ', int2str(TrialNumInt)) %Prints the error message to the command window.
            
            
            try
                
                peakLoc = maxtab_IMUV(:,1);
                peakValues = maxtab_IMUV(:,2);
                
                temp = maxk(peakValues,3);
                
                
                %Name of variables are named in order of highest peak.
                peak1 = peakLoc(find(peakValues == temp(1), 1 ,'first'));
                peak2 = peakLoc(find(peakValues == temp(2), 1 ,'first'));
                peak3 = peakLoc(find(peakValues == temp(3), 1 ,'first'));
                
                
                %Time From Onset to Peaks (for Velocity)
                
                TimeFromOnsetToLastPeak = max([peak1 peak2 peak3])-onset_IMUV;
                TimeFromOnsetToFirstPeak = min([peak1 peak2 peak3])-onset_IMUV;
                
                
                TimeFromOnsetToHighestPeak =  peakLoc(find(peakValues == max(temp), 1 ,'first'))-onset_IMUV;
                TimeFromOnsetToLowestPeak =  peakLoc(find(peakValues == min(temp), 1 ,'first'))-onset_IMUV;
                
                
            catch ME
                
                strcat('Error in script: 3 Peaks Not Detected.  Subject: ', int2str(SubjectNumInt), ' Trial: ', int2str(TrialNumInt)) %Prints the error message to the command window.
                
            end
        end
        
    end
    
    
    %Insert more variables of interest here.
    %Copy the outcome array (line 715) to line 733 of this file.
    %Add the variable to the declarations at the top of this file.
    %Also change the length (add 1 for each variable added) of the Master.m: line 47 (data3), line 70 (data3), line 89 (second data5)
    %Add the name of the variable
    
    outcome = [SubjectNumInt, TrialNumInt, salimu, mtimu,  jerkimu_dim,  jerkimu_dim_log, ImaxVel, ImeanVelocity, ImaxAccel,  ImeanAccel,   stdAccel, TimeFromOnsetToLastPeak,   TimeFromOnsetToFirstPeak,    TimeFromOnsetToHighestPeak,   TimeFromOnsetToLowestPeak, corrvaluex, RMSE, errorVelocity, errorScript];
    csvwrite(strcat('Output_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv'),outcome);
    
    
catch ME
    errorScript = 1;
    strcat('Error in script: Subject: ', int2str(SubjectNumInt), ' Trial: ', int2str(TrialNumInt)) %Prints the error message to the command window.
    outcome = [SubjectNumInt, TrialNumInt, salimu, mtimu,  jerkimu_dim,  jerkimu_dim_log, ImaxVel, ImeanVelocity, ImaxAccel,  ImeanAccel,   stdAccel, TimeFromOnsetToLastPeak,   TimeFromOnsetToFirstPeak,    TimeFromOnsetToHighestPeak,   TimeFromOnsetToLowestPeak, corrvaluex, RMSE, errorVelocity, errorScript];
    csvwrite(strcat('Output_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv'),outcome);
end


if plot_mode == 1
    f1.Visible='on';
    saveas(f1,strcat('f1_Raw_IMU_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.png'))
    f3.Visible='on';
    saveas(f3,strcat('f3_Velocity_and_Acceleration',TaskTitle,'_',SubjectNum,'_', TrialNum,'.png'))
    f68.Visible = 'on';
    saveas(f68,strcat('f68_Velocities_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.png'))
    f69.Visible='on';
    saveas(f69,strcat('f69_Acceleration_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.png'))
    f70.Visible = 'on';
    saveas(f70,strcat('f70_Resultant_Velocity_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.png'))
    
end

