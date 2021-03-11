clear all;
close all;
filename = 'Data_sets_12.2.20/serial_20201222_130309.csv';
data = load(filename);  
% data = data_raw(:,2:7);
csvwrite('Tran1222_1_1.csv',data);

% close all; clear; SensorKinematicsE('Sen', 12,3, 2, 1);
% Replace 'Sen' with the Task Title you used.
% Replace 12 with the Subject Number.
% Replace 3 with the Trial Number.
% Replace 2 with the Activity Code (1 for Block, 2 for Drink(soda can), 3 for Water pouring).
% Leave the 1. This will keep the figures and plots enabled.
SensorKinematicsE('Tran1222',1,1, 1, 1);

clear all;
close all;
filename = 'Data_sets_12.2.20/serial_20201222_130247.csv';
data = load(filename);  
csvwrite('Tran1222_1_2.csv',data);
SensorKinematicsE('Tran1222',1,2, 1, 1);

clear all;
close all;
filename = 'Data_sets_12.2.20/serial_20201222_130218.csv';
data = load(filename);  
csvwrite('Tran1222_1_3.csv',data);
SensorKinematicsE('Tran1222',1,3, 1, 1);

clear all;
close all;
filename = 'Data_sets_12.2.20/serial_20201222_130142.csv';
data = load(filename);  
csvwrite('Tran1222_1_4.csv',data);
SensorKinematicsE('Tran1222',1,4, 1, 1);

clear all;
close all;
filename = 'Data_sets_12.2.20/serial_20201222_130203.csv';
data = load(filename);  
csvwrite('Tran1222_1_5.csv',data);
SensorKinematicsE('Tran1222',1,5, 1, 1);

close all;