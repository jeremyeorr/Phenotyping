clear all
close all

% Add paths
addpath('C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\New analysis_Campana')
choices = {'Veup', 'LG', 'Var', 'Vpassive', 'Vactive', 'GG analysis', 'GG cal', 'Pcrit', 'Alt Vact', 'Cancel'};

% Load organized data file
SaveDir = ['C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\New analysis_Campana\Patient_data']; 
cd(SaveDir)
[fname, pathname, Findex] = uigetfile('*.mat', 'Select File'); % load file
cd(pathname)
name = strcat(pathname,fname);
load(fname); % Loads file
x = x';

GG = Data.GG; 
Fs = Fs.GG;
clear Data Time

%% Calibrate GG
indxGGcal = find(strcmp(pick, choices{7}));
DCremoveWindow = 0.1; %seconds
SmoothWindow = 0.2; %seconds
if ~isempty(GG)
    tempGG = GG(ceil(x(indxGGcal,1)/Fs): floor(x(indxGGcal,2)/Fs));
    timeGG = [ceil(x(indxGGcal,1)/Fs)*Fs: Fs: floor(x(indxGGcal,2)/Fs)*Fs]';
    [scale, offset] = CalibrateGG(Fs, timeGG, tempGG, DCremoveWindow, SmoothWindow);
end
title('Press any key to continue')
pause 
close all
save GGcals scale offset DCremoveWindow SmoothWindow
%% Use scale and offset and user input of gains to scale GG accordingly

% note single(GG) means that we have rounded to back to 7 decimal places
[t, Gain, GGscaled] = ScaleGG(Fs, single(GG), scale, offset, DCremoveWindow,SmoothWindow);

save GGnew GGscaled t Gain
