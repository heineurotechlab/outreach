clear all;
close all;
filename = 'Fixed_units_data_set/Fixed_units_data_set/serial_20210117_001209.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('TranFixed0117_1_1.csv',data);
% close all; clear; SensorKinematicsE('Sen', 12,3, 2, 1);
% Replace 'Sen' with the Task Title you used.
% Replace 12 with the Subject Number.
% Replace 3 with the Trial Number.
% Replace 2 with the Activity Code (1 for Block, 2 for Drink(soda can), 3 for Water pouring).
% Leave the 1. This will keep the figures and plots enabled.
SensorKinematicsE('TranFixed0117',1,1, 1, 1);

close all;
filename = 'Fixed_units_data_set/Fixed_units_data_set/serial_20210117_001237.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('TranFixed0117_1_2.csv',data);
SensorKinematicsE('TranFixed0117',1,2, 1, 1);

close all;
filename = 'Fixed_units_data_set/Fixed_units_data_set/serial_20210117_001256.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('TranFixed0117_1_3.csv',data);
SensorKinematicsE('TranFixed0117',1,3, 1, 1);

close all;
filename = 'Fixed_units_data_set/Fixed_units_data_set/serial_20210117_001313.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('TranFixed0117_1_4.csv',data);
SensorKinematicsE('TranFixed0117',1,4, 1, 1);

close all;
filename = 'Fixed_units_data_set/Fixed_units_data_set/serial_20210117_001345.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('TranFixed0117_1_5.csv',data);
SensorKinematicsE('TranFixed0117',1,5, 1, 1);
