

function SensorKinematicsAutoImport(TaskTitle, SubjectNumInt, TrialNumInt, Activity, trimStart, trimEnd, plot_mode, stationaryMidThreshold)
%% Clearing Workspace
%clc; clear;close all;

%% Import Libraries
addpath('Quaternions');

%% Initialize Reference/Output Varibles

errorStationary = 0;
errorVelocity = 0;
errorScript = 0;
salvicon = 0; salimu = 0;
mtvicon = 0; mtimu = 0;
VmaxVel = 0; ImaxVel = 0;
corrvaluex= 0; RMSE = 0; peakdetThres = 0;

%% File Name Parameters:

if SubjectNumInt <10
    SubjectNum = strcat('0',int2str(SubjectNumInt));
else
    SubjectNum = int2str(SubjectNumInt);
end

TrialNum = int2str(TrialNumInt);

%% IMU Data Entry (Paste Data or Read From File)
% Load IMU File
data = csvread(strcat('Data/',TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv'), 1,0); % IMU file name will have two '_'
data = data(max([trimStart 1]):(end-trimEnd),:);

if(   ((SubjectNumInt==21) && (TrialNumInt ==37))    || ((SubjectNumInt==31) && (TrialNumInt >= 61))   || ((SubjectNumInt==12) && (TrialNumInt == 35))   || ((SubjectNumInt==14) && (TrialNumInt == 14))  || ((SubjectNumInt==18) && (TrialNumInt == 43))  )
    data = a;%Intentionally Triggering an Error For Bad Trial. (Bad Sync or Bad Vicon)
end


if(  ((SubjectNumInt==1) && (TrialNumInt ==9))    )
    data = data(281:end,:);
end

if(   ((SubjectNumInt==1) && (TrialNumInt ==6))    )
    data = data(211:end,:);
end




%% Vicon Data Entry (Paste Data or Read From File)
% Load Vicon File
data1 = csvread(strcat('Data/',TaskTitle,SubjectNum,'_', TrialNum,'.csv'), 5,2); % Vicon file name will have one '_'
data1 = data1(max([trimStart 1]):(end-trimEnd),:);


if(    ((SubjectNumInt==1) && (TrialNumInt ==9))    )
    data1 = data1(281:end,:);
end

if(   ((SubjectNumInt==1) && (TrialNumInt ==6))    )
    data1 = data1(211:end,:);
end


% figure(401)
% clf;
% plot(100/(length(data)-1)*(0:(length(data)-1)),abs(fft(data(1:end,:))));
% 
% 
% figure(402)
% clf;
% plot(100/(916-275)*(0:(916-275)),abs(fft(data(275:916,:))));


try %try/catch will run the catch block when an error occurs in the try block.
    %% Pre-Processing Filter for Drink and Water
    
%     if (Activity == 2) || (Activity == 3)
%         % LP filter accelerometer data
%         filtCutOff = 10;
%         [b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
%         data = filtfilt(b, a, data);
%     end
    %% Automatic Trimming of IMU
    
    
    g = 9.81;
    
    axTemp = (data(:,1)/g);
    ayTemp = (data(:,2)/g);
    azTemp = (data(:,3)/g+1);
    
    gyrXTemp = data(:,4);
    gyrYTemp = data(:,5);
    gyrZTemp = data(:,6);
    
    ac1_mag = sqrt(  (data(:,1)/g).*(data(:,1)/g)   +   (data(:,2)/g).*(data(:,2)/g)  +   (data(:,3)/g+1).*(data(:,3)/g+1)  );
    trimTempL1 = max([1 min([find( azTemp<0.7, 1,'first') find(abs(gyrXTemp)>10,1,'first') find(abs(gyrYTemp)>10,1,'first') find(abs(gyrZTemp)>10,1,'first') find( abs(ac1_mag-1) > 0.1 , 1,'first')])  ]);
    trimLeft = max([1 find( abs(ac1_mag(1:trimTempL1)-1) < 0.025 , 1,'last')-100]);
    
    trimTemp1 = min([length(data)  max( [find( azTemp<0.7, 1,'last') find(abs(gyrXTemp)>5,1,'last') find(abs(gyrYTemp)>5,1,'last') find(abs(gyrZTemp)>5,1,'last')])]);
    trimTemp2 = min([find( abs(ac1_mag(trimTemp1:end)-1)>0.1,1,'last')+trimTemp1 length(data)] );
    trimTemp3 = min([find( abs(ac1_mag(trimTemp2:end)-1)<0.02,1,'last')+trimTemp2 length(data)] );
    trimRight = min([find( abs(ac1_mag((trimTemp2+100):trimTemp3)-1)<0.0025,1,'last')+trimTemp2 length(data)] );
    
    
    if plot_mode == 1
        f300 = figure(300);
        subplot(2,1,1)
        plot( ac1_mag , 'k')
        hold on;
        plot( data(:,1)/g, 'r');
        plot( data(:,2)/g, 'g');
        plot( data(:,3)/g+1, 'b');
        plot(trimLeft, ac1_mag(trimLeft), 'g*', 'MarkerSize',20);
        plot(trimRight, ac1_mag(trimRight), 'r*', 'MarkerSize',20);
        plot(trimTempL1, ac1_mag(trimTempL1), 'ko' , 'MarkerSize',20);
        plot(trimTemp1, ac1_mag(trimTemp1), 'bo' , 'MarkerSize',20);
        plot(trimTemp2, ac1_mag(trimTemp2), 'go' , 'MarkerSize',20);
        plot(trimTemp3, ac1_mag(trimTemp3), 'ro' , 'MarkerSize',20);
        legend('Raw A Resultant Magnitude', 'Raw Ax', 'Raw Ay', 'Raw Az',  'Left Trim' , 'Right Trim', 'TempL1','TempR1', 'TempR2', 'TempR3')
        
        subplot(2,1,2)
        hold on;
        plot(abs(ac1_mag-1))
        legend('abs(A Resultant Magnitude)-1')
        
        
    end
    
    data = data(trimLeft:trimRight,:);
    data1 = data1(trimLeft:min([trimRight length(data1)]),:);
    
    
    %% IMU Data pre-processing on Import
    samplePeriod = 1/100;
    
    g = 9.81;
    accX = data(:,1)/g;
    accY = data(:,2)/g;
    accZ = data(:,3)/g+1;
    
    gyrX = data(:,4);
    gyrY = data(:,5);
    gyrZ = data(:,6);
    
    % Magnetometer(:,1) = data(:,7)*10/1000;
    % Magnetometer(:,2) = data(:,8)*10/1000;
    % Magnetometer(:,3) = data(:,9)*10/1000;% Convert mGauss to Gauss
    
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
    %% Detect stationary periods
    
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
    for stationaryStart = 0.025:0.002:0.1
        stationary = acc_magFilt < stationaryStart;
        if ( (mean(stationary(1:10)))<1  )
            
        else
            break;
        end
    end
    
    if plot_mode ==1
        stationaryStart
    end
    
    w = min([find(stationary<1,1) find(abs(gyrX)>10,1,'first') find(abs(gyrY)>10,1,'first') find(abs(gyrZ)>10,1,'first')]);
    
    for stationaryStop = 0.025:0.002:0.1
        stationary = acc_magFilt < stationaryStop;
        if (   (mean(stationary((end-20):end))) < 1  )
            
        else
            break;
        end
    end
    
    if plot_mode ==1
        stationaryStop
    end
    
    stationary(1:w) = 1;
    r = max([find(stationary<1,1, 'last') find(abs(gyrX)>10,1,'last') find(abs(gyrY)>10,1,'last') find(abs(gyrZ)>10,1,'last')]);
    %r = find(stationary<1,1,'last');
    mt_imu_with_stationary = r-w;
    
    
    
    
    if Activity == 1
        stationary(w:r) = acc_magFilt(w:r) < stationaryMidThreshold;
    elseif Activity == 2
        stationary(w:r) = acc_magFilt(w:r) < stationaryMidThreshold;
    elseif Activity == 3
        stationary(w:r) = acc_magFilt(w:r) < stationaryMidThreshold;
    else
        strcat('Subject: ',SubjectNum,' Trial: ', TrialNum,'   Activity Error')
    end
    
    
    if plot_mode == 1
        % -------------------------------------------------------------------------
        % Plot data raw sensor data and stationary periods
        %figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Sensor Data');
        f1 = figure(1);
        clf;
        f1.Visible = 'off';
        f1.Name = 'Sensor Data';
        %f1.Position = [9 39 900 600];
        
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
        clf;
        f3.Visible = 'off';
        f3.Name = 'Velocity and Position';
        %f3.Position = [9 39 900 300];
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
    
    IMU_x = pos(1,:);
    IMU_y = pos(2,:);
    IMU_z = pos(3,:);
    
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
    
    %% Assigning Filtered IMU Data to new variable. Also Vicon Data is assigned to appropriate variables.
    
    IMU_vx1 = -vel(:,1);
    IMU_vy1 = -vel(:,2);
    IMU_vz1 = vel(:,3);
    
    IMU_ax1 = -acc(:,1);
    IMU_ay1 = -acc(:,2);
    IMU_az1 = acc(:,3);
    
    %Vicon Data
    Vicon_x = data1(:,1)/1000;
    Vicon_y = data1(:,2)/1000;
    Vicon_z = data1(:,3)/1000;
    
    %Derivative Operator on Vicon Data
    %velocity and acceleration
    Vicon_Vx1 = centdiff(Vicon_x,0.01);
    Vicon_Ax1 = centdiff(Vicon_Vx1,0.01);
    Vicon_Vy1 = centdiff(Vicon_y,0.01);
    Vicon_Ay1 = centdiff(Vicon_Vy1,0.01);
    Vicon_Vz1 = centdiff(Vicon_z,0.01);
    Vicon_Az1 = centdiff(Vicon_Vz1,0.01);
    Vicon_x = Vicon_x - Vicon_x(1);
    Vicon_y = Vicon_y - Vicon_y(1);
    Vicon_z = Vicon_z - Vicon_z(1);
    
    selena_x = figure(142857)
    subplot(2,1,1)
    plot(Vicon_x);
    title("Vicon X");
    subplot(2,1,2)
    plot(IMU_x);
    title("IMU X");
    selena_y = figure(428571)
    subplot(2,1,1)
    plot(Vicon_y);
    title("Vicon Y");
    subplot(2,1,2)
    plot(IMU_y);
    title("IMU Y");
    selena_z = figure(285714)
    subplot(2,1,1)
    plot(Vicon_y);
    title("Vicon Y");
    subplot(2,1,2)
    plot(IMU_y);
    title("IMU Y");
    
    %% Offset Adjustment for Vicon
    offset = 0; %Offset Observed Visually
    
    %% Implementing the Offset (if corrected visually) and Fixing the Length
    
    if offset >= 0
        
        if length(IMU_ax1)>length(Vicon_Ax1)
            IMU_ax = IMU_ax1(1:length(Vicon_Ax1)-offset);
            Vicon_Ax = Vicon_Ax1(1+offset:length(Vicon_Ax1));
            
            IMU_ay = IMU_ay1(1:length(Vicon_Ay1)-offset);
            Vicon_Ay = Vicon_Ay1(1+offset:length(Vicon_Ay1));
            
            IMU_az = IMU_az1(1:length(Vicon_Az1)-offset);
            Vicon_Az = Vicon_Az1(1+offset:length(Vicon_Az1));
            
            IMU_vx = IMU_vx1(1:length(Vicon_Vx1)-offset);
            Vicon_Vx = Vicon_Vx1(1+offset:length(Vicon_Vx1));
            
            IMU_vy = IMU_vy1(1:length(Vicon_Vy1)-offset);
            Vicon_Vy = Vicon_Vy1(1+offset:length(Vicon_Vy1));
            
            IMU_vz = IMU_vz1(1:length(Vicon_Vz1)-offset);
            Vicon_Vz = Vicon_Vz1(1+offset:length(Vicon_Vz1));
        else
            IMU_ax = IMU_ax1(1:length(IMU_ax1)-offset);
            Vicon_Ax = Vicon_Ax1(1+offset:length(IMU_ax1));
            
            IMU_ay = IMU_ay1(1:length(IMU_ay1)-offset);
            Vicon_Ay = Vicon_Ay1(1+offset:length(IMU_ay1));
            
            IMU_az = IMU_az1(1:length(IMU_az1)-offset);
            Vicon_Az = Vicon_Az1(1+offset:length(IMU_az1));
            
            
            
            IMU_vx = IMU_vx1(1:length(IMU_vx1)-offset);
            Vicon_Vx = Vicon_Vx1(1+offset:length(IMU_vx1));
            
            IMU_vy = IMU_vy1(1:length(IMU_vy1)-offset);
            Vicon_Vy = Vicon_Vy1(1+offset:length(IMU_vy1));
            
            IMU_vz = IMU_vz1(1:length(IMU_vz1)-offset);
            Vicon_Vz = Vicon_Vz1(1+offset:length(IMU_vz1));
        end
    else
        offset = -offset;
        if length(IMU_ax1)>length(Vicon_Ax1)
            IMU_ax = IMU_ax1(1+offset:length(Vicon_Ax1));
            Vicon_Ax = Vicon_Ax1(1:length(Vicon_Ax1)-offset);
            
            IMU_ay = IMU_ay1(1+offset:length(Vicon_Ay1));
            Vicon_Ay = Vicon_Ay1(1:length(Vicon_Ay1)-offset);
            
            IMU_az = IMU_az1(1+offset:length(Vicon_Az1));
            Vicon_Az = Vicon_Az1(1:length(Vicon_Az1)-offset);
            
            
            
            IMU_vx = IMU_vx1(1+offset:length(Vicon_Vx1));
            Vicon_Vx = Vicon_Vx1(1:length(Vicon_Vx1)-offset);
            
            IMU_vy = IMU_vy1(1+offset:length(Vicon_Vy1));
            Vicon_Vy = Vicon_Vy1(1:length(Vicon_Vy1)-offset);
            
            IMU_vz = IMU_vz1(1+offset:length(Vicon_Vz1));
            Vicon_Vz = Vicon_Vz1(1:length(Vicon_Vz1)-offset);
        else
            IMU_ax = IMU_ax1(1+offset:length(IMU_ax1));
            Vicon_Ax = Vicon_Ax1(1:length(IMU_ax1)-offset);
            
            IMU_ay = IMU_ay1(1+offset:length(IMU_ay1));
            Vicon_Ay = Vicon_Ay1(1:length(IMU_ay1)-offset);
            
            IMU_az = IMU_az1(1+offset:length(IMU_az1));
            Vicon_Az = Vicon_Az1(1:length(IMU_az1)-offset);
            
            
            IMU_vx = IMU_vx1(1+offset:length(IMU_vx1));
            Vicon_Vx = Vicon_vx1(1:length(IMU_vx1)-offset);
            
            IMU_vy = IMU_vy1(1+offset:length(IMU_vy1));
            Vicon_Vy = Vicon_Vy1(1:length(IMU_vy1)-offset);
            
            IMU_vz = IMU_vz1(1+offset:length(IMU_vz1));
            Vicon_Vz = Vicon_Vz1(1:length(IMU_vz1)-offset);
        end
        
    end
    
    %% Data import, Domenico's Filter, and length correction happen above.
    
    if plot_mode == 1
        f69 = figure(69);
        clf;
        f69.Name = 'Sensor Data: Accelerations';
        f69.Visible = 'off';
        plot(sqrt(Vicon_Ax.^2 + Vicon_Ay.^2 +Vicon_Az.^2))
        hold on;
        plot(sqrt(IMU_ax.^2 + IMU_ay.^2 + IMU_az.^2))
        title('Acceleration');
        legend('Vicon','IMU');
    end
    
    
    ViconV = (sqrt(Vicon_Vx.^2 + Vicon_Vy.^2 +Vicon_Vz.^2));
    IMUV = (sqrt(IMU_vx.^2 + IMU_vy.^2 +IMU_vz.^2));
    

    
    if IMUV(10) > 0.1 || IMUV(end) > 0.1
        errorVelocity = 1;
        %     error_trial = strcat(TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv')
        %     error_type = 'Velocity Bad'
    end
    
    
    % Change peakdet threshold until the same number of peaks are detected by
    % Vicon and IMU
    % if waterpouring == 1
    %     peakdetThres = 0.2;
    % else
    %     peakdetThres = 0.4;
    % end
    
    %
    %
    % for i = 1:10
    %     %Detect Peaks in Time Series (wmc)
    %     %Vicon
    %     [maxtab_ViconV, mintab_ViconV] = peakdet(ViconV, peakdetThres);
    %
    %     %Detect Peaks in Time Series (wmc)
    %     %IMU
    %     [maxtab_IMUV, mintab_IMUV] = peakdet(IMUV, peakdetThres);
    %
    %
    %     if length(maxtab_ViconV) == length(maxtab_IMUV)
    %         break;
    %     else
    %         peakdetThres = peakdetThres*0.9;
    %     end
    %
    %     if peakdetThres < 0.1
    %        break;
    %     end
    %
    % end
    
    
    % Change peakdet threshold until the same number of peaks are detected by
    % Vicon and IMU
    % if waterpouring == 1
    %     peakdetThres = 0.1;
    % else
    %     peakdetThres = 0.2;
    %
    % end
    
    peakdetThres = 0.1;
    
    peakdetThres4 = 0;
    peakdetThres3 = 0;
    
    for i = 1:70
        %Detect Peaks in Time Series (wmc)
        %Vicon
        [maxtab_ViconV, mintab_ViconV] = peakdet(ViconV, peakdetThres);
        
        %Detect Peaks in Time Series (wmc)
        %IMU
        [maxtab_IMUV, mintab_IMUV] = peakdet(IMUV, peakdetThres);

        
        if length(maxtab_ViconV)>=3 && length(maxtab_IMUV)>=3 && i<70
            
            if (length(maxtab_ViconV)==4 && length(maxtab_IMUV)==4)
                peakdetThres4 = peakdetThres;
            end
            
            if (length(maxtab_ViconV)==3 && length(maxtab_IMUV)==3)
                peakdetThres3 = peakdetThres;
            end
            
            peakdetThres = peakdetThres*1.05;
            
        elseif length(maxtab_ViconV)== length(maxtab_IMUV)
            break;
        else
            if i > 1
                peakdetThres = peakdetThres/1.05;
            end
            break;
        end
    end
    
    if Activity == 1
        if (peakdetThres3 >0 && peakdetThres4 >0) 
            peakdetThres = max([peakdetThres4 peakdetThres3]);
        elseif peakdetThres3 >0
            peakdetThres=peakdetThres3;
        elseif peakdetThres4 >0
            peakdetThres=peakdetThres4;
        end
    else
        if peakdetThres4 >0
            peakdetThres=peakdetThres4;
        elseif peakdetThres3 >0
            peakdetThres=peakdetThres3;
        end
    end
    
    [maxtab_ViconV, mintab_ViconV] = peakdet(ViconV, peakdetThres);
    [maxtab_IMUV, mintab_IMUV] = peakdet(IMUV, peakdetThres);
    
    
    
    if plot_mode == 1
        
        f68 = figure(68);
        clf;
        f68.Visible = 'off';
        f68.Name = 'Resultant Velocities';
        hold on;
        plot (ViconV, 'k'); hold on
        plot(maxtab_ViconV(:,1), maxtab_ViconV(:,2), 'k*');hold on; title('Resultant Velocity')
        
    end
    
    MinDuration = 50;
    
    %Onset
    if any(ViconV(maxtab_ViconV(1,1):-1:1)<maxtab_ViconV(2,2)*0.02)
        for onset_ViconV = maxtab_ViconV(1,1):-1:1
            
            if onset_ViconV - MinDuration > 1
                minZeroVelocityFrom =  onset_ViconV - MinDuration;
            else
                minZeroVelocityFrom = 1;
            end
            
            if (ViconV(onset_ViconV)<maxtab_ViconV(2,2)*0.02)  && mean(ViconV(minZeroVelocityFrom:onset_ViconV))<0.05
                break
            end
        end
    else
        for onset_ViconV = maxtab_ViconV(1,1):-1:1
            if onset_ViconV - MinDuration > 1
                minZeroVelocityFrom =  onset_ViconV - MinDuration;
            else
                minZeroVelocityFrom = 1;
            end
            if (ViconV(onset_ViconV)< maxtab_ViconV(2,2)*0.05)  && mean(ViconV(minZeroVelocityFrom:onset_ViconV))<0.05
                break
            end
        end
    end
    
    
    MinDuration = 75;
    %Offset of Vicon
    if any(ViconV(maxtab_ViconV(length(maxtab_ViconV),1):1:end)<maxtab_ViconV(length(maxtab_ViconV),2)*0.02)
        for offset_ViconV = (maxtab_ViconV(length(maxtab_ViconV),1):1:length(ViconV)-1)
            if offset_ViconV + MinDuration < length(ViconV)
                minZeroVelocityUntil =  offset_ViconV + MinDuration;
            else
                minZeroVelocityUntil = length(ViconV);
            end
            ViconV_diff = centdiff(ViconV,1/100);
            if (   ((ViconV(offset_ViconV)< maxtab_ViconV(length(maxtab_ViconV),2)*0.02) || (ViconV(offset_ViconV)< 0.05))  &&   (mean(ViconV(offset_ViconV:minZeroVelocityUntil))<0.1)  )
                break
            end
        end
    elseif any(ViconV(maxtab_ViconV(length(maxtab_ViconV),1):1:size(ViconV)<maxtab_ViconV(length(maxtab_ViconV),2)*0.05))
        for offset_ViconV = (maxtab_ViconV(length(maxtab_ViconV),1):1:length(ViconV)-1)
            if offset_ViconV + MinDuration < length(ViconV)
                minZeroVelocityUntil =  offset_ViconV + MinDuration;
            else
                minZeroVelocityUntil = length(ViconV);
            end
            
            if (   ((ViconV(offset_ViconV)< maxtab_ViconV(length(maxtab_ViconV),2)*0.05) || (ViconV(offset_ViconV)< 0.05))  &&   (mean(ViconV(offset_ViconV:minZeroVelocityUntil))<0.1)  )
                break
            end
        end
    else
        [s, loc_offset_ViconV] = min(ViconV(maxtab_ViconV(length(maxtab_ViconV),1):1:length(ViconV)-1));
        offset_ViconV = loc_offset_ViconV + maxtab_ViconV(length(maxtab_ViconV),1) -1;
    end
    
    
    if plot_mode == 1
        
        plot(onset_ViconV, ViconV(onset_ViconV), 'g*');
        plot(offset_ViconV, ViconV(offset_ViconV), 'g*');
    end
    
    %IMU
    
    if plot_mode == 1
        plot (IMUV, 'b'); hold on
        plot(maxtab_IMUV(:,1), maxtab_IMUV(:,2), 'k*');hold on; title(' Velocity')
        
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
            
            if (IMUV(onset_IMUV)<maxtab_IMUV(2,2)*0.01 ) && (mean(IMUV(minZeroVelocityFrom:onset_IMUV))<0.05)
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
            
            
            if (IMUV(onset_IMUV)< maxtab_IMUV(2,2)*0.05) && (mean(IMUV(minZeroVelocityFrom:onset_IMUV))<0.05)
                break
            end
        end
    end

    
    MinDuration = 75;
    %Offset of IMU
    if any(IMUV(maxtab_IMUV(length(maxtab_IMUV),1):1:end)<maxtab_IMUV(length(maxtab_IMUV),2)*0.02)
        for offset_IMUV = (maxtab_IMUV(length(maxtab_IMUV),1):1:length(IMUV)-1)
            if offset_IMUV + MinDuration < length(IMUV)
                minZeroVelocityUntil =  offset_IMUV + MinDuration;
            else
                minZeroVelocityUntil = length(IMUV);
            end
            IMUV_diff = centdiff(IMUV,1/100);
            if (   ((IMUV(offset_IMUV)< maxtab_IMUV(length(maxtab_IMUV),2)*0.02) || (IMUV(offset_IMUV)< 0.05))  &&   (mean(IMUV(offset_IMUV:minZeroVelocityUntil))<0.05)  )
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
            
            if (   ((IMUV(offset_IMUV)< maxtab_IMUV(length(maxtab_IMUV),2)*0.05) || (IMUV(offset_IMUV)< 0.05))  &&   (mean(IMUV(offset_IMUV:minZeroVelocityUntil))<0.05)   )
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
        legend('Vicon', 'Vicon Peak' , 'Vicon Onset', 'Vicon Offset', 'IMU', 'IMU Peak' , 'IMU Onset', 'IMU Offset');
    end
    
    
    %% Trim Time Series
    %Plot Overall
    trimV = ViconV(onset_ViconV:offset_ViconV);
    trimI = IMUV(onset_ViconV:offset_ViconV);
    
    
    if plot_mode == 1
        f70 = figure(70);
        clf;
        f70.Visible = 'off';
        f70.Name = 'Trimmed Resultant Data';
        plot(trimV,'k'), hold on; plot(trimI,'b')
        title('Resultant Velocity');
        legend('Vicon','IMU');
    end
    
    
    %Correlation and RMSE
    corrx = corrcoef(trimV,trimI);
    corrvaluex = corrx (1,2);
    
    RMSE = sqrt(mean((trimV-trimI).^2));
    
    %Trim all axes for both IMU and Vicon
    trimVx = Vicon_Vx(onset_ViconV:offset_ViconV);
    trimVy = Vicon_Vy(onset_ViconV:offset_ViconV);
    trimVz = Vicon_Vz(onset_ViconV:offset_ViconV);
    
    trimVx_accel = Vicon_Ax(onset_ViconV:offset_ViconV);
    trimVy_accel = Vicon_Ay(onset_ViconV:offset_ViconV);
    trimVz_accel = Vicon_Az(onset_ViconV:offset_ViconV);
    
    trimIx = IMU_vx(onset_ViconV:offset_ViconV);
    trimIy = IMU_vy(onset_ViconV:offset_ViconV);
    trimIz = IMU_vz(onset_ViconV:offset_ViconV);
    
    trimIx_accel = IMU_ax(onset_ViconV:offset_ViconV);
    trimIy_accel = IMU_ay(onset_ViconV:offset_ViconV);
    trimIz_accel = IMU_az(onset_ViconV:offset_ViconV);
    
    
    if plot_mode == 1
        f101=figure (101);
        clf;
        f101.Name = 'Trimmed Velocity Data';
        f101.Visible = 'off';
        
        subplot (3,1,1)
        plot (trimVx,'k'); hold on; plot (trimIx,'b');
        title('Velocities')
        legend ('Vicon X','IMU X');
        
        
        subplot (3,1,2)
        plot (trimVy,'k'); hold on; plot (trimIy,'b');
        legend ('Vicon Y','IMU Y');
        
        subplot (3,1,3)
        plot (trimVz,'k'); hold on; plot (trimIz,'b');
        legend ('Vicon Z','IMU Z');
        
        
        f102=figure (102);
        clf;
        f102.Name = 'Trimmed Acceleration Data';
        f102.Visible = 'off';
        
        subplot (3,1,1)
        plot (trimVx_accel,'k'); hold on; plot (trimIx_accel,'b');
        title('Accelerations')
        legend ('Vicon X','IMU X');
              
        subplot (3,1,2)
        plot (trimVy_accel,'k'); hold on; plot (trimIy_accel,'b');
        legend ('Vicon Y','IMU Y');
        
        subplot (3,1,3)
        plot (trimVz_accel,'k'); hold on; plot (trimIz_accel,'b');
        legend ('Vicon Z','IMU Z');
        
        
    end
    
    %----------------------------Kinematic Variables-------------------
    %Movement Time
    mtvicon = (offset_ViconV-onset_ViconV); %value depends on sampling rate
    mtimu = (offset_IMUV-onset_IMUV); %value depends on sampling rate
    mtdelta = abs(mtvicon -mtimu);
    if(plot_mode ==1)
        mtdelta
    end
  
    
    
    %Movement Smoothness
    salvicon = func_SAL(trimVx, trimVy, trimVz);
    [jerkvicon_dim,jerkvicon_dim_log] = func_jerk_dim(trimVx, trimVy, trimVz);
    
    salimu = func_SAL(trimIx, trimIy, trimIz);
    [jerkimu_dim,jerkimu_dim_log] = func_jerk_dim(trimIx, trimIy, trimIz);
    
    %Max Velocity
    VmaxVel = max(ViconV);
    ImaxVel = max(IMUV);
    
    %Insert more variables of interest here.
    %Copy the outcome array (line 715) to line 733 of this file.
    %Add the variable to the declarations at the top of this file.
    %Also change the length (add 1 for each variable added) of the Master.m: line 47 (data3), line 70 (data3), line 89 (second data5)
    %Add the name of the variable
    outcome = [SubjectNumInt, TrialNumInt, salvicon, salimu, mtvicon, mtimu, VmaxVel, ImaxVel, corrvaluex, RMSE, peakdetThres, errorStationary, errorVelocity, errorScript];
    
    csvwrite(strcat('Output/','Output_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv'),outcome);
    
    
    if plot_mode == 1
        f200 = figure(200);
        f200.Visible = 'off';
        f200.Name = 'Trimmed Velocity Data';
        plot(onset_ViconV:1:offset_ViconV,ViconV(onset_ViconV:offset_ViconV) ,'k')
        hold on;
        plot(onset_IMUV:1:offset_IMUV,IMUV(onset_IMUV:offset_IMUV),'b')
        
        legend ('Vicon Resultant Velocity', 'IMU Resultant Velocity');
    end
    
    
    
catch ME
    errorScript = 1;
    outcome = [SubjectNumInt, TrialNumInt, salvicon, salimu, mtvicon, mtimu, VmaxVel, ImaxVel, corrvaluex, RMSE, peakdetThres, errorStationary, errorVelocity, errorScript];
    csvwrite(strcat('Output/','Output_',TaskTitle,'_',SubjectNum,'_', TrialNum,'.csv'),outcome);
end


if plot_mode == 1
    
    %     f1 = figure(1);
    %     f1.Position = [1374 -137 560 420];
    %     f300 = figure(300);
    %     f300.Position = [2488 -163 700 942];
    %     f200 = figure(200);
    %     f200.Position = [1930 -137 560 420];
    %     f69 = figure(69);
    %     f69.Position = [1366 -141 560 420];
    %     f68 = figure(68);
    %     f68.Position = [1927 357 560 420];
    
    
    f1 = figure(1);
%     f1.Position = [1281 285 560 420];
    f300 = figure(300);
%     f300.Position =  [1850 285 700 942];
    f200 = figure(200);
%     f200.Position = [2500 -279 700 942];
    f68 = figure(68);
%     f68.Position = [1281 -258 560 420];
    
    
    
    f1.Visible='on';
    f3.Visible='on';
    f68.Visible = 'on';
    f69.Visible='on';
    f70.Visible = 'on';
    f101.Visible='on';
    f102.Visible='on';
    f200.Visible = 'on';
    f300.Visible = 'on';
    
    
%     f101.Visible = 'off';
%     f102.Visible = 'off';
%     f70.Visible = 'off';
%     f3.Visible = 'off';
%     f69.Visible='off';
    
    
end

