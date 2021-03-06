clear all;
close all;
filename = 'Raw_data_sets/Raw_data_sets/raw_16b_0.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran1224raw_1_0.csv',data);
% close all; clear; SensorKinematicsE('Sen', 12,3, 2, 1);
% Replace 'Sen' with the Task Title you used.
% Replace 12 with the Subject Number.
% Replace 3 with the Trial Number.
% Replace 2 with the Activity Code (1 for Block, 2 for Drink(soda can), 3 for Water pouring).
% Leave the 1. This will keep the figures and plots enabled.
SensorKinematicsE('Tran1224raw',1,0, 1, 1);

clear all;
close all;
filename = 'Raw_data_sets/Raw_data_sets/raw_16b_1.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran1224raw_1_1.csv',data);
SensorKinematicsE('Tran1224raw',1,1, 1, 1);

clear all;
close all;
filename = 'Raw_data_sets/Raw_data_sets/raw_16b_2.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran1224raw_1_2.csv',data);
SensorKinematicsE('Tran1224raw',1,2, 1, 1);

clear all;
close all;
filename = 'Raw_data_sets/Raw_data_sets/raw_16b_3.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran1224raw_1_3.csv',data);
SensorKinematicsE('Tran1224raw',1,3, 1, 1);

clear all;
close all;
filename = 'Raw_data_sets/Raw_data_sets/raw_16b_4.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran1224raw_1_4.csv',data);
SensorKinematicsE('Tran1224raw',1,4, 1, 1);
