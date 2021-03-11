clear all;
close all;
%% Import Libraries
addpath('../');

addpath('../Quaternions');
% close all; clear; SensorKinematicsE('Sen', 12,3, 2, 1);
% Replace 'Sen' with the Task Title you used.
% Replace 12 with the Subject Number.
% Replace 3 with the Trial Number.
% Replace 2 with the Activity Code (1 for Block, 2 for Drink(soda can), 3 for Water pouring).
% Leave the 1. This will keep the figures and plots enabled.
SensorKinematicsE('Gondar',1,14, 1, 1);
SensorKinematicsE('Gondar',1,15, 1, 1);
SensorKinematicsE('Gondar',1,16, 1, 1);