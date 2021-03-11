
%% Importing Raw Data
clear; close all; clc;

dataIn = readtable('GondarTrialsDay1.csv');
TaskTitle = 'Gondar';
SubNum = table2array(dataIn(:,1));
TrialNum = table2array(dataIn(:,2));
Hand = table2array(dataIn(:,3));
Impaired = table2array(dataIn(:,4));
Activity = table2array(dataIn(:,5));

MasterStatus = 'Processing Trials'


for i = 1:length(SubNum)
    
    try
        ActivityInt = 0;
        
        if Activity(i) =="Block"
            ActivityInt = 1;
        elseif Activity(i) =="Drink"
            ActivityInt = 2;
        elseif Activity(i) =="Water"
             ActivityInt = 3;
        elseif Activity(i) =="BicepCurl"
            ActivityInt = 4;
        end
        
        SensorKinematicsE(TaskTitle, SubNum(i), TrialNum(i),ActivityInt , 0);
        
        %         disp('Press Any Key to Continue');
        %         pause;
        %         disp('Now Continuing');
        
    catch ME
        strcat(TaskTitle,' Subject: ', int2str(SubNum(i)), ' Trial: ', int2str(TrialNum(i)), '   IMU File Problem (Not Found, Blank, etc)' )
        
        output = [SubNum(i) TrialNum(i) zeros(1,17)];
        csvwrite(strcat('Output_',TaskTitle,'_',int2str(SubNum(i)),'_', int2str(TrialNum(i)),'.csv'),output);
        
        %         disp('Press Any Key to Continue');
        %         pause;
        %         disp('Now Continuing');
        
    end
    clearvars -except TaskTitle SubNum TrialNum Hand Activity Impaired;
    
    close all;
end

%% Combining CSVs

close all;

data3 = zeros(1,19);

for i = 1:length(SubNum)
    data4 = csvread(strcat('Output_',TaskTitle,'_',int2str(SubNum(i)),'_',int2str(TrialNum(i)),'.csv'));
    data3 = [data3; data4];
end

data5 = array2table(data3(2:size(data3,1), :) ,'VariableNames', {'SubjectNumber' 'TrialNumber' 'SalImu' 'MtImu' 'jerkimu_dim' 'jerkimu_dim_log' 'ImuMaxVel' 'ImeanVelocity' 'ImaxAccel' 'ImeanAccel' 'IstdAccel'  'TimeFromOnsetToLastPeak'   'TimeFromOnsetToFirstPeak'    'TimeFromOnsetToHighestPeak'   'TimeFromOnsetToLowestPeak'  'corrvaluex' 'RMSE' 'ErrorVelocity' 'ErrorScript'});
outputCombined = [data5(:, 1:2) array2table(Activity,'VariableNames', {'Activity'}) array2table(Hand, 'VariableNames',{'Hand'}) array2table(Impaired,'VariableNames', {'Impaired'}) data5(:,3:19)];
writetable(outputCombined , strcat('CombinedOutput_',TaskTitle,'_',datestr(now, 'mm-dd-YYYY HH-MM PM'),'.csv' ));

MasterStatus = 'Processing Complete'
