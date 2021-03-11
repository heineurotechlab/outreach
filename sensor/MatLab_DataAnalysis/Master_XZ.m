clear all;
close all;
filename = 'Tran_data/serial_20201209_095935.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_1.csv',data);

% close all; clear; SensorKinematicsE('Sen', 12,3, 2, 1);
% Replace 'Sen' with the Task Title you used.
% Replace 12 with the Subject Number.
% Replace 3 with the Trial Number.
% Replace 2 with the Activity Code (1 for Block, 2 for Drink(soda can), 3 for Water pouring).
% Leave the 1. This will keep the figures and plots enabled.
SensorKinematicsE('Tran',1,1, 1, 1);

clear all;
close all;
filename = 'Tran_data/serial_20201209_095918.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_2.csv',data);
SensorKinematicsE('Tran',1,2, 1, 1);

clear all;
close all;
filename = 'Tran_data/serial_20201209_095902.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_3.csv',data);
SensorKinematicsE('Tran',1,3, 1, 1);

clear all;
close all;
filename = 'Tran_data/serial_20201209_095840.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_4.csv',data);
SensorKinematicsE('Tran',1,4, 1, 1);

clear all;
close all;
filename = 'Tran_data/serial_20201209_095821.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_5.csv',data);
SensorKinematicsE('Tran',1,5, 1, 1);

clear all;
close all;
filename = 'Tran_data/serial_20201209_095803.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_6.csv',data);
SensorKinematicsE('Tran',1,6, 1, 1);

clear all;
close all;
filename = 'Tran_data/serial_20201209_095742.txt';
data_raw = load(filename);  
data = data_raw(:,2:7);
csvwrite('Tran_1_7.csv',data);
SensorKinematicsE('Tran',1,7, 1, 1);