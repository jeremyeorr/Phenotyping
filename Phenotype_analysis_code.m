%% Phenotype traits analysis
% created Brad Edwards, PhD (August 2014)
% based on original code by Lisa M. Campana, PhD (December 2012)

% Loads in data file with labeled times you want to analyze  
% Calls functions Veup.m, Varousal.m, Vzero.m, Uag.m, LoopGain.m, PcritCalc.m, GGanalysis.m
%% Add paths, load spike data and save as Matlab file for analysis
clear all
close all
clc
% Add paths
%addpath('C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode')
addpath('C:\Users\j1orr\Dropbox\Code\MATLab\Phenotyping')
%Traits events codes and labels - if defaults setting are used in Spike
% 1 - Veup
% 2 - Vpass
% 3 - Var
% 4 - Vact
% 5 - LG
% 6 - GGcal
% 7 - Pcrit
% 8 - GGr

%Change RouteDir to the location of IdentifyTimePeriodsofInterest.m 
%RouteDir = 'C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode';
RouteDir = 'C:\Users\j1orr\Dropbox\Code\MATLab\Phenotyping';
addpath(RouteDir)

%Specify Directories where you want to load and then save data

SaveDir = ['S:\RESEARCH\STUDIES\3. DATA\Endophenotyping studies\Phenotyping night'];

LoadDir = ['S:\RESEARCH\STUDIES\3. DATA\Endophenotyping studies\Phenotyping night'];


%User inputs the number of files 
[n] = inputdlg('Number of Files'); 

%Load Data file(s)
cd(LoadDir)
[Data, Fs, Time] = LoadData(str2num(n{1}));

%change directory to where you would like to save data
CurrentDir = cd;
temp=findstr(CurrentDir, '\');
folder = CurrentDir(temp(end)+1:end);
if ~exist([SaveDir, '\', folder],'dir')
    mkdir(SaveDir,folder);
end
cd([SaveDir, '\', folder])
clear temp folder 

%Saves the imported data with Traits scored (saves importing again)
saveFile = 'Imported_data.mat';
save(saveFile, 'Data','Fs','Time','LoadDir','RouteDir','CurrentDir','n');

%% Load already saved data
%Specify Directories where you want to load already saved data

clear all 
close all
clc
addpath('C:\Users\j1orr\Dropbox\Code\MATLab\Phenotyping')

Loadsaveddata = ['S:\RESEARCH\STUDIES\3. DATA\Endophenotyping studies\Phenotyping night'];
[fname, pathname, Findex] = uigetfile('*.mat', 'Select File'); % load file
cd(pathname) 
name = strcat(pathname,fname);
load(fname); % Loads file!
if exist('GGnew.mat') %Loads GGnew file is already exists
   load('GGnew.mat')
else
    warning('No scaled GG was imported')
end
clear Loadsaveddata pathname fname name Findex

%% Calibrating the GGrecordings
GG = Data.GG; 
FsGG = Fs.GG;
TimeGG = Time.Vdot;

% Find location of calibration
indxGGcal = find(Data.Traits1Evts==6);
if ~isempty(GG)
    tempGG = GG(ceil(Data.Traits1EvtTimeStart(indxGGcal,1)/FsGG): floor(Data.Traits1EvtTimeEnd(indxGGcal,1)/FsGG));
    timeGG = [ceil(Data.Traits1EvtTimeStart(indxGGcal,1)/FsGG)*FsGG: FsGG: floor(Data.Traits1EvtTimeEnd(indxGGcal,1)/FsGG)*FsGG]';
    [scale, offset] = CalibrateGG(FsGG, timeGG, tempGG);
end
title('Press any key to continue')
pause 
close all
save GGcals scale offset 

% Use scale and offset and user input of gains to scale GG accordingly
% saves new GGsignal expressed as %max throughout night
[t, Gain, GGscaled] = ScaleGG(FsGG, GG, scale, offset);
save GGnew GGscaled t Gain
clear GG FsGG TimeGG indxGGcal scale offset timeGG tempGG

%% Veupnea analysis
% Find all locations were you want to analyze Veupnea
indxVeup=find(Data.Traits1Evts==1); 

%Loop over all windows of interest
for i = 1:length(indxVeup);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxVeup(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxVeup(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxVeup(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxVeup(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxVeup(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxVeup(i),1));% finding start and end of any arousals that occured
    
    %create temp GG signal over time period of interset
    if exist('GGscaled')
         tempGG = GGscaled(ceil(Data.Traits1EvtTimeStart(indxVeup(i),1)/Fs.GG): floor(Data.Traits1EvtTimeEnd(indxVeup(i),1)/Fs.GG));
         timeGG = [ceil(Data.Traits1EvtTimeStart(indxVeup(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits1EvtTimeEnd(indxVeup(i),1)/Fs.GG)*Fs.GG]';
    else
        %warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    %create temp Pepi, PO2, PCO2 over time period of interest
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    if isempty(Data.PO2)
        pO2temp = [];
    else
        pO2temp= Data.PO2(a);
    end
    if isempty(Data.PCO2)
        pCO2temp = [];
    else
        pCO2temp= Data.PCO2(a);
    end
    
    %Call Veup.m and calculate variables over the time
    %period of interest
    [Veupnea.VdotE(i), Veupnea.Cpap(i), Veupnea.Pepi(i), Veupnea.PO2(i), Veupnea.PCO2(i),Veupnea.GGphasic(i), Veupnea.GGtonic(i), Veupnea.GGpeak(i),Veupnea.Ti(i),Veupnea.Te(i),Veupnea.Ttot(i),Veupnea.VT(i)]...
        = Veup(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, pO2temp, pCO2temp, timeGG, tempGG);
    Veupnea.Time(i,:) = Data.Traits1EvtTimeStart(indxVeup(i),:);

    subplot(8,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVeup))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp pCO2temp pO2temp
 
end

%Automatically determining if its oxygen or placebo night 
if nanmean(Veupnea.PO2)>0 & nanmean(Veupnea.PO2)<180
    night = 'Baseline';
elseif nanmean(Veupnea.PO2)>180 
% (if PO2 is over %180 assume it is the treatment night)
    night = 'Treatment';
else
    night = 'unknown'; %unknown
end

% Calculating average values for summary data
Veupnea_averageddata.VdotE=nanmean(Veupnea.VdotE);
Veupnea_averageddata.Cpap=nanmean(Veupnea.Cpap);
Veupnea_averageddata.Pepi=-nanmean(Veupnea.Pepi);%added -ve sign to make Pepi a neagative value
Veupnea_averageddata.PO2=nanmean(Veupnea.PO2);
Veupnea_averageddata.PCO2=nanmean(Veupnea.PCO2);
Veupnea_averageddata.GGphasic=nanmean(Veupnea.GGphasic);
Veupnea_averageddata.GGtonic=nanmean(Veupnea.GGtonic);
Veupnea_averageddata.GGpeak=nanmean(Veupnea.GGpeak);
Veupnea_averageddata.Ti=nanmean(Veupnea.Ti);
Veupnea_averageddata.Te=nanmean(Veupnea.Te);
Veupnea_averageddata.Ttot=nanmean(Veupnea.Ttot);
Veupnea_averageddata.VT=nanmean(Veupnea.VT);
Veupnea_meandata=[nanmean(Veupnea.VdotE) nanmean(Veupnea.Cpap) -nanmean(Veupnea.Pepi) nanmean(Veupnea.PO2) nanmean(Veupnea.PCO2) nanmean(Veupnea.GGphasic) nanmean(Veupnea.GGtonic) nanmean(Veupnea.GGpeak) nanmean(Veupnea.Ti) nanmean(Veupnea.Te) nanmean(Veupnea.Ttot) nanmean(Veupnea.VT) numel(Veupnea.VdotE)]

%Save variables to patient folder
save Veupnea Veupnea night Veupnea_meandata

%% Varousal analysis
% Find all locations were you want to analyze Varousal
indxVar=find(Data.Traits1Evts==3);

%Loop over all windows of interest
for i = 1:length(indxVar);
    try
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxVar(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxVar(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxVar(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxVar(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxVar(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxVar(i),1));
    
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    if exist('GGscaled')
         tempGG = GGscaled(ceil(Data.Traits1EvtTimeStart(indxVar(i),1)/Fs.GG): floor(Data.Traits1EvtTimeEnd(indxVar(i),1)/Fs.GG));
         timeGG = [ceil(Data.Traits1EvtTimeStart(indxVar(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits1EvtTimeEnd(indxVar(i),1)/Fs.GG)*Fs.GG]';
    else
        %warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    [Var.VdotE(i), Var.Pepi(i), Var.GGphasic(i), Var.GGtonic(i), Var.GGpeak(i), Var.Ti(i), Var.Te(i), Var.Ttot(i), Var.VT(i), Var.AvgCPAP(i), Var.tarousal(i)]...
        = Varousal(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
    Var.Time(i,:) = Data.Traits1EvtTimeStart(indxVar(i),:);
    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVar))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp
    catch me
        disp(me.message)
    end
end

%Deletes all Pepi values that are not at least 2cmH2O higher than its
%resting value
Pepithreshold=nanmean(Veupnea.Pepi)+2

% Calculating average values for summary data
Var_averageddata.VdotE=nanmean(Var.VdotE)
Var_averageddata.Pepi=-nanmean(Var.Pepi(Var.Pepi>Pepithreshold)) % throws out Pepi's less than threshold
Var_averageddata.GGphasic=nanmean(Var.GGphasic)
Var_averageddata.GGtonic=nanmean(Var.GGtonic)
Var_averageddata.GGpeak=nanmean(Var.GGpeak)
Var_averageddata.Ti=nanmean(Var.Ti)
Var_averageddata.Te=nanmean(Var.Te)
Var_averageddata.Ttot=nanmean(Var.Ttot)
Var_averageddata.VT=nanmean(Var.VT)
Var_averageddata.AvgCPAP=nanmean(Var.AvgCPAP)
Var_meandata=[Var_averageddata.VdotE Var_averageddata.Pepi Var_averageddata.GGphasic Var_averageddata.GGtonic Var_averageddata.GGpeak Var_averageddata.Ti Var_averageddata.Te Var_averageddata.Ttot Var_averageddata.VT Var_averageddata.AvgCPAP numel(Var.VdotE)]
Var_mediandata=[nanmedian(Var.VdotE) -nanmedian(Var.Pepi(Var.Pepi>Pepithreshold)) nanmedian(Var.GGphasic) nanmedian(Var.GGtonic) nanmedian(Var.GGpeak) nanmedian(Var.Ti) nanmedian(Var.Te) nanmedian(Var.Ttot) nanmedian(Var.VT) nanmedian(Var.AvgCPAP)] 
save Var Var Var_meandata Var_mediandata Pepithreshold

%% Vpassive analysis
% Find all locations were you want to analyze Vpassive
indxV0=find(Data.Traits1Evts==2);

%Loop over all windows of interest
for i = 1:length(indxV0);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxV0(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxV0(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxV0(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxV0(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxV0(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxV0(i),1));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    if exist('GGscaled')
         tempGG = GGscaled(ceil(Data.Traits1EvtTimeStart(indxV0(i),1)/Fs.GG): floor(Data.Traits1EvtTimeEnd(indxV0(i),1)/Fs.GG));
         timeGG = [ceil(Data.Traits1EvtTimeStart(indxV0(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits1EvtTimeEnd(indxV0(i),1)/Fs.GG)*Fs.GG]';
    else
        %warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    [V0.VdotE(i), V0.Pepi(i), V0.GGphasic(i), V0.GGtonic(i), V0.GGpeak(i), V0.Ti(i), V0.Te(i), V0.Ttot(i), V0.VT(i), V0.CPAPdrop(i)]...
        = Vzero(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
    V0.Time(i,:) = Data.Traits1EvtTimeStart(indxV0(i),:);
    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxV0))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp
end

V0_meandata=[mean(V0.VdotE) -nanmean(V0.Pepi) nanmean(V0.GGphasic) nanmean(V0.GGtonic) nanmean(V0.GGpeak) numel(V0.VdotE)]
save V0 V0 V0_meandata
%% Vpassive analysis MARK 2 - only need to use when no Pepi AND zero flow
% Find all locations were you want to analyze Vpassive
indxV0=find(Data.Traits1Evts==2);

%Loop over all windows of interest
for i = 1:length(indxV0);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxV0(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxV0(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxV0(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxV0(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxV0(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxV0(i),1));

    [V0.VdotE(i), V0.Ti(i), V0.Te(i), V0.Ttot(i), V0.VT(i), V0.CPAPdrop(i)]...
        = VzeroALT(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)]);
    V0.Time(i,:) = Data.Traits1EvtTimeStart(indxV0(i),:);
    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxV0))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* 
end

V0_meandata=[mean(V0.VdotE) numel(V0.VdotE)]
save V0 V0 V0_meandata
%% LG analysis
%This code produces 4 LG measurements - the first 2 averages breath 1 and 2
%and the last 2 uses only breath 1 to determine the ventilatory response. 
%It also produces LG calculations using Lisa's method for
%leak removal and another using Andrews method (see Line 94-95 of LoopGain
%m-file) 

indxLG = find(Data.Traits1Evts==5);

for i = 1:length(indxLG);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxLG(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxLG(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxLG(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxLG(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxLG(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxLG(i),1));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    [LG_2br.LGcampana(i), LG_2br.LGwellman(i), LG_1br.LGcampana(i), LG_1br.LGwellman(i), LG_2br.Vrespcampana(i), LG_1br.Vrespcampana(i), LG_1br.Vdistcampana(i), LG_2br.Vrespwellman(i), LG_1br.Vrespwellman(i), LG_1br.Vdistwellman(i), LG_1br.CPAPmin(i)]...
        = LoopGain(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, mean(Veupnea.VdotE));
    LG_1br.Time(i,:) = Data.Traits1EvtTimeStart(indxLG(i),:);

    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxLG))]; 'Press Enter to continue'})
    pause;
    close all
    clear a*
end
% Prints following information:
% 1. LGs determined using 2 breaths (Lisa & Andrew)
% 2. LGs determined using 1 breath (Lisa & Andrew)
% 3. Total number of measurements and number of non-Nan'd measurements
LGsummary=[nanmean(LG_2br.LGcampana) nanmean(LG_2br.LGwellman) nanmean(LG_1br.LGcampana) nanmean(LG_1br.LGwellman) nanmean(LG_2br.Vrespcampana) nanmean(LG_1br.Vrespcampana) nanmean(LG_1br.Vdistcampana) nanmean(LG_2br.Vrespwellman) nanmean(LG_1br.Vrespwellman) nanmean(LG_1br.Vdistwellman) nanmean(LG_1br.CPAPmin) numel(LG_1br.LGwellman) sum(~isnan(LG_1br.LGwellman))] %shows LG results in workspace
LGsummarymedians=[nanmedian(LG_2br.LGcampana) nanmedian(LG_2br.LGwellman) nanmedian(LG_1br.LGcampana) nanmedian(LG_1br.LGwellman) nanmedian(LG_2br.Vrespcampana) nanmedian(LG_1br.Vrespcampana) nanmedian(LG_1br.Vdistcampana) nanmedian(LG_2br.Vrespwellman) nanmedian(LG_1br.Vrespwellman) nanmedian(LG_1br.Vdistwellman) nanmean(LG_1br.CPAPmin)] 
save LG LG_2br LG_1br LGsummary LGsummarymedians 

%% Vactive analysis
indxVact = find(Data.Traits1Evts==4);

for i = 1:length(indxVact);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxVact(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxVact(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxVact(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    if exist('GGscaled')
         tempGG = GGscaled(ceil(Data.Traits1EvtTimeStart(indxVact(i),1)/Fs.GG): floor(Data.Traits1EvtTimeEnd(indxVact(i),1)/Fs.GG));
         timeGG = [ceil(Data.Traits1EvtTimeStart(indxVact(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits1EvtTimeEnd(indxVact(i),1)/Fs.GG)*Fs.GG]';
    else
        warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    if isempty(Data.PO2)
        pO2temp = [];
    else
        pO2temp= Data.PO2(a);
    end
    [Vactive.VdotE(i), Vactive.Vdrive(i), Vactive.Ug(i), Vactive.Cpap(i), Vactive.Pepi(i), Vactive.GGphasic(i),Vactive.GGtonic(i), Vactive.GGpeak(i), Vactive.Ti(i), Vactive.Te(i), Vactive.Ttot(i), Vactive.VT(i), Vactive.Vdotpredrop(i), Vactive.Pepipredrop(i), Vactive.CPAPdrop(i)]=...
        Uag(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG, nanmean(Veupnea.VdotE), nanmean(LG_1br.LGwellman), nanmean(V0.VdotE), pO2temp);
    
    Vactive.Time(i,:) = Data.Traits1EvtTimeStart(indxVact(i),:);

    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVact))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp

end

save Vactive Vactive 
%{
%Comparing Vactive data with pepi and peakGG at arousal
temp1=Var_meandata(2); %mean Pepi at arousal
temp2=Var_meandata(5); %mean GGpeak at arousal
temp3=nanmean(Var.AvgCPAP)
figure(1)
ax1(1)=subplot(1,3,1); plot(Vactive.Cpap, Vactive.VdotE, 'b.', 'MarkerSize',20)
ylabel('VdotE (L/min)')
xlabel('active CPAP level(cmH2O)')
title({['CPAP at Ar: ', num2str(sprintf('%.1f',temp3))]})
ax1(2)=subplot(1,3,2); plot(Vactive.Cpap, Vactive.GGpeak, 'g.', 'MarkerSize',20)
ylabel('GGpeak (%max)')
xlabel('active CPAP level(cmH2O)')
title({['GGpeak at Ar: ', num2str(sprintf('%.1f',temp2))]})
ax1(3)=subplot(1,3,3); plot(Vactive.Cpap, Vactive.Pepi, 'r.', 'MarkerSize',20)
ylabel('Pepi(%cmH2O)')
xlabel('active CPAP level(cmH2O)')
title({['Pepi at Ar: ', num2str(sprintf('%.1f',temp1))]})
saveas(figure(1), 'Comparing_Vactive_with_Varousal');
clear temp1 temp2 temp3
%}

% Only using data that is from the lowest CPAP level
if length(Vactive.Cpap)>3
    tempCpap = round(Vactive.Cpap); 
    [minCpap] = min(tempCpap); 
    minInd = find(tempCpap==minCpap); 
    if length(minInd)>3
        CPAPmin = Vactive.Cpap(minInd);
        %Vactlow = mean(Vactive(minInd));
    else
        temp = find(tempCpap>minCpap);
        minInd = find(tempCpap==minCpap | tempCpap==min(tempCpap(temp)));
        CPAPmin = Vactive.Cpap(minInd);
    end
else
    minInd = [1:length(Vactive.Cpap)];
    CPAPmin = Vactive.Cpap;
end
VactiveVE = nanmean(Vactive.VdotE(minInd)); %previously VactiveLow
VactiveVdrive = nanmean(Vactive.Vdrive(minInd));
VactiveUAG = nanmean(Vactive.Ug(minInd));
VactiveGGpeak = nanmean(Vactive.GGpeak(minInd));
VactiveGGtonic = nanmean(Vactive.GGtonic(minInd));
VactivePepi = nanmean(Vactive.Pepi(minInd));
VactiveVdotpredrop = nanmean(Vactive.Vdotpredrop(minInd));
VactivePepipredrop = nanmean(Vactive.Pepipredrop(minInd));
CPAPminavg=nanmean(Vactive.Cpap(minInd));
Nodropsused=numel(minInd);
Totaldropsperformed=numel(Vactive.Cpap);
Vactive_summary=[VactiveVE CPAPminavg VactiveVdrive VactiveUAG VactiveGGpeak VactiveGGtonic -VactivePepi VactiveVdotpredrop -VactivePepipredrop Nodropsused Totaldropsperformed]
save Vactivesummary  Vactive_summary minInd CPAPmin
clear tempCpap minCpap minInd Nodropsused Totaldropsperformed VactiveVE VactiveVdrive VactiveUAG VactiveGGpeak VactiveGGtonic VactivePepi VactiveVdotpredrop VactivePepipredrop Nodropsused Totaldropsperformed

%% Vactive analysis no pepi/GG
indxVact = find(Data.Traits1Evts==4);
for i = 1:length(indxVact);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxVact(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxVact(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxVact(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    if isempty(Data.PO2)
        pO2temp = [];
    else
        pO2temp= Data.PO2(a);
    end
    
  % Allows Vactive to run whether or not LG has been measured (won't give
  % UAG however)- JEO
  
 
    if exist('LG', 'var');
    % note this uses wellman 1 breath LG to determine UAG
    [Vactive.VdotE(i), Vactive.Vdrive(i), Vactive.Ug(i), Vactive.Cpap(i), Vactive.Ti(i), Vactive.Te(i), Vactive.Ttot(i), Vactive.VT(i), Vactive.Vdotpredrop(i), Vactive.CPAPdrop(i)]=...
        Uag_Alt(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], nanmean(Veupnea.VdotE), nanmean(LG_1br.LGwellman), nanmean(V0.VdotE), pO2temp);
  
    else
    [Vactive.VdotE(i), Vactive.Vdrive(i), Vactive.Ug(i), Vactive.Cpap(i), Vactive.Ti(i), Vactive.Te(i), Vactive.Ttot(i), Vactive.VT(i), Vactive.Vdotpredrop(i), Vactive.CPAPdrop(i)]=...
      Uag_Alt(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], nanmean(Veupnea.VdotE), NaN, nanmean(V0.VdotE), pO2temp);
    end  
  
  Vactive.Time(i,:) = Data.Traits1EvtTimeStart(indxVact(i),:);

    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVact))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp

end

save Vactive Vactive 

% Only using data that is from the lowest CPAP level
if length(Vactive.Cpap)>3
    tempCpap = round(Vactive.Cpap); 
    [minCpap] = min(tempCpap); 
    minInd = find(tempCpap==minCpap); 
    if length(minInd)>3
        CPAPmin = Vactive.Cpap(minInd);
        %Vactlow = mean(Vactive(minInd));
    else
        temp = find(tempCpap>minCpap);
        minInd = find(tempCpap==minCpap | tempCpap==min(tempCpap(temp)));
        CPAPmin = Vactive.Cpap(minInd);
    end
else
    minInd = [1:length(Vactive.Cpap)];
    CPAPmin = Vactive.Cpap;
end
VactiveVE = nanmean(Vactive.VdotE(minInd)); %previously VactiveLow
VactiveVdrive = nanmean(Vactive.Vdrive(minInd));
VactiveUAG = nanmean(Vactive.Ug(minInd));
VactiveVdotpredrop = nanmean(Vactive.Vdotpredrop(minInd));
CPAPminavg=nanmean(Vactive.Cpap(minInd));
Nodropsused=numel(minInd);
Totaldropsperformed=numel(Vactive.Cpap);
Vactive_summary=[VactiveVE CPAPminavg VactiveVdrive VactiveUAG VactiveVdotpredrop  Nodropsused Totaldropsperformed]
save Vactivesummary  Vactive_summary minInd CPAPmin
clear tempCpap minCpap minInd Nodropsused Totaldropsperformed VactiveVE VactiveVdrive VactiveUAG VactiveVdotpredrop Nodropsused Totaldropsperformed
%% Vactive baseline analysis
indxVact = find(Data.Traits1Evts==4);
    
for i = 1:length(indxVact);
    a = find(Time.Vdot>=Data.Traits1EvtTimeStart(indxVact(i),1) & Time.Vdot<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    a2 = find(Time.Epochs>=Data.Traits1EvtTimeStart(indxVact(i),1)-30 & Time.Epochs<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    a3 = find(Data.EvtTimeEnd>=Data.Traits1EvtTimeStart(indxVact(i),1) & Data.EvtTimeStart<=Data.Traits1EvtTimeEnd(indxVact(i),1));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    if exist('GGscaled')
         tempGG = GGscaled(ceil(Data.Traits1EvtTimeStart(indxVact(i),1)/Fs.GG): floor(Data.Traits1EvtTimeEnd(indxVact(i),1)/Fs.GG));
         timeGG = [ceil(Data.Traits1EvtTimeStart(indxVact(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits1EvtTimeEnd(indxVact(i),1)/Fs.GG)*Fs.GG]';
    else
        %warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end

[Vactbaseline.VdotE(i), Vactbaseline.Pepi(i), Vactbaseline.GGphasic(i), Vactbaseline.GGtonic(i), Vactbaseline.GGpeak(i), Vactbaseline.Ti(i), Vactbaseline.Te(i), Vactbaseline.Ttot(i), Vactbaseline.VT(i), Vactbaseline.CPAPmin(i)]=...
    Vactivebaseline(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
  
    Vactbaseline.Time(i,:) = Data.Traits1EvtTimeStart(indxVact(i),:);
    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVact))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp

end
save Vactivebaseline Vactbaseline 

%% Pcrit analysis
indxPcrit=find(Data.Traits2Evts==7);
DetectThreshold=0.5; %Lisa had set to 0.75
Pcrit = struct('Cpap', [], 'VdotMax', [], 'VdotE', [], 'Time', [], 'pc', []);
for i = 1:length(indxPcrit)
        a = find(Time.Vdot>=Data.Traits2EvtTimeStart(indxPcrit(i),1) & Time.Vdot<=Data.Traits2EvtTimeEnd(indxPcrit(i),1));
        a2 = find(Time.Epochs>=Data.Traits2EvtTimeStart(indxPcrit(i),1)-30 & Time.Epochs<=Data.Traits2EvtTimeEnd(indxPcrit(i),1));
        a3 = find(Data.EvtTimeEnd>=Data.Traits2EvtTimeStart(indxPcrit(i),1) & Data.EvtTimeStart<=Data.Traits2EvtTimeEnd(indxPcrit(i),1));
        [Pcrit(i).Cpap, Pcrit(i).VdotMax, Pcrit(i).VdotE, Pcrit(i).Time, Pcrit(i).pc Vcrit(i).vc]...
            = PcritCalc(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], DetectThreshold);
        pause
        close all
        clear a*
end

%Plot Pmask-VE graph
if ~isempty(indxPcrit)
    fh= figure;
    hold on
    Pall = [];
    VdotEall = [];
    for i = 1:size(Pcrit,2)
        ind = find(Pcrit(i).VdotE>0);
        VdotEtemp = Pcrit(i).VdotE(ind);
        PmaskTemp = Pcrit(i).Cpap(ind);
        Pall = [Pall; PmaskTemp];
        VdotEall = [VdotEall; VdotEtemp];
    end
    plot(Pall, VdotEall,'ko')
    ylabel('VdotE (L/min)')
    xlabel('Pmask (cmH2O)')

    p = polyfit(Pall, VdotEall, 1);
    xtemp = [0:max(Pall)+1];
    ytemp = xtemp.*p(1)+p(2);
    plot(xtemp, ytemp, 'k--');
    Vcritall=p(2);

    title({'VdotE average' ; ['VdotE at 0 pressure: ', num2str(p(2))]})
    saveas(fh, 'VcritAll');
    save Pcrit Pcrit Vcrit VdotEall Pall Vcritall
end

%Replot Pmask-Vpeak graph with all data
if ~isempty(indxPcrit)
    fh= figure;
    hold on
    Pall2 = [];
    VdotMaxall = [];
    for i = 1:size(Pcrit,2)
        ind = find(Pcrit(i).VdotMax>0);
        VdotMaxtemp = Pcrit(i).VdotMax(ind);
        PmaskTemp2 = Pcrit(i).Cpap(ind);
        Pall2 = [Pall2; PmaskTemp2];
        VdotMaxall = [VdotMaxall; VdotMaxtemp];
    end
    plot(Pall2, VdotMaxall,'k.')
    ylabel('VdotMax (L/min)')
    xlabel('Pmask (cmH2O)')

    p = polyfit(Pall2, VdotMaxall, 1);
    xtemp = [0:max(Pall2)+1];
    ytemp = xtemp.*p(1)+p(2);
    plot(xtemp, ytemp, 'k--')
    PcritAll=-p(2)/p(1);
    
    title({'VdotMax average' ; ['Pcrit: ', num2str(PcritAll)]})
    saveas(fh, 'PcritAll');
    save Pcrit Pcrit Vcrit Pall2 VdotMaxall PcritAll -append
end

% Easily export variabiles:
for i=1:length(Vcrit)
Vcrit_export(i) = Vcrit(i).vc;
end
for i=1:length(Pcrit)
Pcrit_export(i) = Pcrit(i).pc;
end
Vcrit_export
Pcrit_export

%% Modified Pcrit/Vcrit analysis (breakpoint analysis to include zero mask pressure)
global Pmaskdata VdotMaxdata VdotEdata fixedslope upperlimit
close all
Pmaskdata=[];
VdotMaxdata=[];
VdotEdata=[];
fixedslope=NaN
upperlimit.on=0
exclude_zeroflowdata=0 % zero indicates include zero data points
includevpassivedata=0 %1 indicates you want to add in the vpassive data
plotpoints=1

include_runs=[1 2] %manually select which runs to include. Eg. [1]-will graph run 1 [1 2]-will graph runs 1 and 2
for i=include_runs
tempsize=size(Pcrit(i).VdotMax);
N=tempsize(1)*tempsize(2);

Pmaskdata=[Pmaskdata reshape(Pcrit(i).Cpap,1,N)]; 
VdotMaxdata=[VdotMaxdata reshape(Pcrit(i).VdotMax,1,N)];
VdotEdata=[VdotEdata reshape(Pcrit(i).VdotE,1,N)];
end

VdotMaxdata(isnan(Pmaskdata))=[];
VdotEdata(isnan(Pmaskdata))=[];
Pmaskdata(isnan(Pmaskdata))=[];

VdotMaxdata(Pmaskdata==0)=[];
VdotEdata(Pmaskdata==0)=[];
Pmaskdata(Pmaskdata==0)=[];

fh=figure(1)
[~,~,PcritSEM,Pcritvalue,PVSlope1]=pcrit_model_runss(Pmaskdata,VdotMaxdata,exclude_zeroflowdata,plotpoints)
saveas(fh, 'PcritAllwithzeros');

if includevpassivedata==1
    VdotEdata=[VdotEdata V0.VdotE];
    Pmaskdata=[Pmaskdata V0.CPAPdrop];
else
end

fh=figure(2)
[VcritSEM,Vcritvalue,~,~,PVSlope2]=pcrit_model_runss(Pmaskdata,VdotEdata,exclude_zeroflowdata,plotpoints)
saveas(fh, 'VcritAllwithzeros');

ModPcrit=[Pcritvalue PcritSEM PVSlope1 Vcritvalue VcritSEM PVSlope2]
save ModifiedPcrit ModPcrit    

%% Measuring Vactive on Vcrit graph
%clear all, 
close all;
%load('Vactive.mat'), load('V0.mat'), load('Veupnea.mat'), load('Vactivebaseline.mat');
global Pmaskdata2 VdotEdata2 fixedslope upperlimit
Pmaskdata2=[];
VdotEdata2=[];
fixedslope=1; %use ModPcrit(6) if you want to fix slope=passive slope - NaN otherwise
upperlimit.on=0;
exclude_zeroflowdata=0; % zero indicates include zero data points
plotpoints=1;
includeVactdrops=1; % Use this setting to include Vactives from zero CPAP
redonebaselineVE=1; % Make=1 if I've remeasured Ve pre-Vactive drop
includeVarousaldata=1; %Make=1 if you want to include Varousal data
% currently does nothing - need to add the Varousal data.......

% Lisa code put zero for Ve before drop (if drop was zero) - this code
% excludes those zeros
if  redonebaselineVE==1
    nonzerodata=find(Vactbaseline.VdotE);
    tempVdotpredrop=Vactbaseline.VdotE(nonzerodata);
    tempCPAP=Vactbaseline.CPAPmin(nonzerodata);
else
    nonzerodata=find(Vactive.Vdotpredrop);
    tempVdotpredrop=Vactive.Vdotpredrop(nonzerodata);
    tempCPAP=Vactive.Cpap(nonzerodata);
end

if includeVactdrops==1
    VdotEdata2=[tempVdotpredrop'; Vactive.VdotE']';
    %runs=length(Vactive.VdotE);
    ptemp=Vactive.CPAPdrop;
    Pmaskdata2=[tempCPAP'; ptemp']'; 
%     ptemp=zeros([1,runs]);
%     Pmaskdata2=[tempCPAP'; ptemp']'; 
else
    VdotEdata2=tempVdotpredrop;
    Pmaskdata2=tempCPAP;
end

hgload('VcritAllwithzeros.fig');
[VcritSEM,Vcritvalue,~,~,PVSlope2]=pcrit_model_runss(Pmaskdata2,VdotEdata2,exclude_zeroflowdata,plotpoints)
plot(Pmaskdata2, VdotEdata2, 'r.', 'MarkerSize',18)
%hold on
%plot(Var.AvgCPAP, Var.VdotE, 'k.', 'MarkerSize',18)
saveas(figure(1),'Alternative Vactive.fig');

% % Revised upper airway gain 
% if exist('ArTHandGUA.mat')
%         load('ArTHandGUA.mat')
%         ta = ArousalThreshold;
%     else
%         prompt = 'What is the arousal threshold (L/min)? ';
%         ta = input(prompt);
% end
AllVactive=nanmean(Vactive.VdotE);
%uag_revised1 = (AllVactive-(V0_meandata(1)))/(ta-(Veupnea_meandata(1))); %[using average of all data(active) value]
%uag_revised2 = (Vcritvalue-(V0_meandata(1)))/(ta-(Veupnea_meandata(1))); %[using Vcrit(active) value]
%ModVactivecrit=[AllVactive Vcritvalue VcritSEM PVSlope2 uag_revised1 uag_revised2]
ModVactivecrit=[AllVactive Vcritvalue VcritSEM PVSlope2]
save ModVactivecrit ModVactivecrit VdotEdata2 Pmaskdata2

%% GG analysis
indxGG=find(Data.Traits2Evts==8);
clear Breaths
if exist('GGscaled')
    for i = 1:length(indxGG)
        a = find(Time.Vdot>=Data.Traits2EvtTimeStart(indxGG(i),1) & Time.Vdot<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        a2 = find(Time.Epochs>=Data.Traits2EvtTimeStart(indxGG(i),1)-30 & Time.Epochs<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        a3 = find(Data.EvtTimeEnd>=Data.Traits2EvtTimeStart(indxGG(i),1) & Data.EvtTimeStart<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        if isempty(Data.Pepi)
            peptemp = [];
        else
            peptemp= Data.Pepi(a);
        end
        tempGG = GGscaled(ceil(Data.Traits2EvtTimeStart(indxGG(i),1)/Fs.GG): floor(Data.Traits2EvtTimeEnd(indxGG(i),1)/Fs.GG));
        timeGG = [ceil(Data.Traits2EvtTimeStart(indxGG(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits2EvtTimeEnd(indxGG(i),1)/Fs.GG)*Fs.GG]';
        [Breaths(i)] = GGanalysis(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
        subplot(7,1,1)
        title({['Number ', num2str(i), ' of ', num2str(length(indxGG))]; 'Press Enter to continue'})
        pause;
        close all
         
        %Calculating GGphasic stuff
        I=~isnan(Breaths(i).VdotE)&~isnan(Breaths(i).GGpeak);
        [fit.VdotEandGGpeak(i,:)]=linreg(Breaths(i).VdotE(I),Breaths(i).GGpeak(I),1,'VE(L/min)','GGpeak(%max)');
        I=~isnan(Breaths(i).VdotE)&~isnan(Breaths(i).GGtonic);
        [fit.VdotEandGGtonic(i,:)]=linreg(Breaths(i).VdotE(I),Breaths(i).GGtonic(I),2,'VE(L/min)','GGtonic(%max)');
        I=~isnan(Breaths(i).Pepi)&~isnan(Breaths(i).GGtonic);
        [fit.PepiandGGtonic(i,:)]=linreg(Breaths(i).Pepi(I),Breaths(i).GGtonic(I),4,'Pepi(cmH2O)','GGtonic(%max)');
        I=~isnan(Breaths(i).Pepi)&~isnan(Breaths(i).GGpeak);
        [fit.PepiandGGpeak(i,:)]=linreg(Breaths(i).Pepi(I),Breaths(i).GGpeak(I),3,'Pepi(cmH2O)','GGpeak(%max)');
        I=~isnan(Breaths(i).Pepi)&~isnan(Breaths(i).VdotE);
        [fit.PepiandVdotE(i,:)]=linreg(Breaths(i).Pepi(I),Breaths(i).VdotE(I),5,'Pepi(cmH2O)','VE(L/min)');

        %Combining all data in one large matrix
    N=length(Breaths);
    Pepi_allbreaths=[];
    VdotE_allbreaths=[];
    GGpeak_allbreaths=[];
    GGtonic_allbreaths=[];
    for i=1:N
        Pepi_allbreaths=[Pepi_allbreaths; Breaths(1,i).Pepi];
        VdotE_allbreaths=[VdotE_allbreaths; Breaths(1,i).VdotE];
        GGpeak_allbreaths=[GGpeak_allbreaths; Breaths(1,i).GGpeak];
        GGtonic_allbreaths=[GGtonic_allbreaths; Breaths(1,i).GGtonic];
    end
        GGcombineddata=[Pepi_allbreaths VdotE_allbreaths GGpeak_allbreaths GGtonic_allbreaths];
        % Graphing and recording all Data combined
        I=~isnan(VdotE_allbreaths)&~isnan(GGpeak_allbreaths);
        [fit.VdotEandGGpeak_combined]=linreg(VdotE_allbreaths(I),GGpeak_allbreaths(I),1,'VE(L/min)','GGpeak(%max)');
        I=~isnan(VdotE_allbreaths)&~isnan(GGtonic_allbreaths);
        [fit.VdotEandGGtonic_combined]=linreg(VdotE_allbreaths(I),GGtonic_allbreaths (I),2,'VE(L/min)','GGtonic(%max)');
        I=~isnan(Pepi_allbreaths)&~isnan(GGtonic_allbreaths);
        [fit.PepiandGGtonic_combined]=linreg(Pepi_allbreaths(I),GGtonic_allbreaths(I),4,'Pepi(cmH2O)','GGtonic(%max)');
        I=~isnan(Pepi_allbreaths)&~isnan(GGpeak_allbreaths);
        [fit.PepiandGGpeak_combined]=linreg(Pepi_allbreaths(I),GGpeak_allbreaths(I),3,'Pepi(cmH2O)','GGpeak(%max)');
        I=~isnan(Pepi_allbreaths)&~isnan(VdotE_allbreaths);
        [fit.PepiandVdotE_combined]=linreg(Pepi_allbreaths(I),VdotE_allbreaths(I),5,'Pepi(cmH2O)','VE(L/min)');
        fh = gcf;
        %scrsz = get(0,'ScreenSize');
        %set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
        saveas(fh, ['GGanalysis_combined']);
        pause;
        close all;
        clear a*
    end

fit.PepiandGGpeak_combined
fit.PepiandGGtonic_combined

clear N Pepi_allbreaths VdotE_allbreaths GGpeak_allbreaths GGtonic_allbreaths
save GGresults Breaths fit GGcombineddata
else
    warning('Scaled GG signal is empty') 
end

%% Summarize and put data together
% Settings
usevcrit=1 % make =1 if you want to use Vcrit instead of Vpassive
prompt = 'What was the study intervention (use apostrophe)? ';
night = input(prompt);

if usevcrit==1
    Vpassivefinal=ModPcrit(4);
    Vactivefinal=ModVactivecrit(2);
else
    Vpassivefinal=nanmean(V0.VdotE);
    Vactivefinal=Vactive_summary(1);    
end
%loads MR data in summary
if exist('GGresults.mat')
         MR = fit.PepiandGGpeak_combined(1);
    else
        warning('No Muscle responsiveness data')
        MR = [];
end
subid = inputdlg('Enter subject Id');
[ArousalThreshold, UAG]=SummarizeData(subid{1}, night, MR, mean(Veupnea.VdotE), nanmean(Var.VdotE), Vactivefinal, Vpassivefinal, Vactive.Ug, -nanmean(LG_1br.LGwellman), Vactive.Cpap, -nanmean(Var.Pepi(Var.Pepi>Pepithreshold)))
save ArTHandGUA ArousalThreshold UAG night subid
%% Write Summary data to excel

WriteResultsToExcel(Veupnea, Var, V0, LG, Vactive, Pcrit, Breaths, fit)

%% GG analysis (graphing data individually)
indxGG=find(Data.Traits2Evts==8);
clear Breaths
if exist('GGscaled')
    for i = 1:length(indxGG)
        a = find(Time.Vdot>=Data.Traits2EvtTimeStart(indxGG(i),1) & Time.Vdot<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        a2 = find(Time.Epochs>=Data.Traits2EvtTimeStart(indxGG(i),1)-30 & Time.Epochs<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        a3 = find(Data.EvtTimeEnd>=Data.Traits2EvtTimeStart(indxGG(i),1) & Data.EvtTimeStart<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        if isempty(Data.Pepi)
            peptemp = [];
        else
            peptemp= Data.Pepi(a);
        end
        tempGG = GGscaled(ceil(Data.Traits2EvtTimeStart(indxGG(i),1)/Fs.GG): floor(Data.Traits2EvtTimeEnd(indxGG(i),1)/Fs.GG));
        timeGG = [ceil(Data.Traits2EvtTimeStart(indxGG(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits2EvtTimeEnd(indxGG(i),1)/Fs.GG)*Fs.GG]';
        [Breaths(i)] = GGanalysis(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
        subplot(7,1,1)
        title({['Number ', num2str(i), ' of ', num2str(length(indxGG))]; 'Press Enter to continue'})
        pause;
        close all
         
        %Calculating GGphasic stuff
        I=~isnan(Breaths(i).VdotE)&~isnan(Breaths(i).GGpeak);
        [fit.VdotEandGGpeak(i,:)]=linreg(Breaths(i).VdotE(I),Breaths(i).GGpeak(I),1,'VE(L/min)','GGpeak(%max)');
        I=~isnan(Breaths(i).VdotE)&~isnan(Breaths(i).GGtonic);
        [fit.VdotEandGGtonic(i,:)]=linreg(Breaths(i).VdotE(I),Breaths(i).GGtonic(I),2,'VE(L/min)','GGtonic(%max)');
        I=~isnan(Breaths(i).Pepi)&~isnan(Breaths(i).GGtonic);
        [fit.PepiandGGtonic(i,:)]=linreg(Breaths(i).Pepi(I),Breaths(i).GGtonic(I),4,'Pepi(cmH2O)','GGtonic(%max)');
        I=~isnan(Breaths(i).Pepi)&~isnan(Breaths(i).GGpeak);
        [fit.PepiandGGpeak(i,:)]=linreg(Breaths(i).Pepi(I),Breaths(i).GGpeak(I),3,'Pepi(cmH2O)','GGpeak(%max)');
        I=~isnan(Breaths(i).Pepi)&~isnan(Breaths(i).VdotE);
        [fit.PepiandVdotE(i,:)]=linreg(Breaths(i).Pepi(I),Breaths(i).VdotE(I),5,'Pepi(cmH2O)','VE(L/min)');
        %fh = gcf;
        %scrsz = get(0,'ScreenSize');
        %set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
        %saveas(fh, ['GGanalysis', num2str(round(Time.Vdot(a(1))))]);
        pause;
        %close all;
        clear a*
    end

    %Combining all data in one large matrix
    N=length(Breaths);
    Pepi_allbreaths=[];
    VdotE_allbreaths=[];
    GGpeak_allbreaths=[];
    GGtonic_allbreaths=[];
    for i=1:N
        Pepi_allbreaths=[Pepi_allbreaths; Breaths(1,i).Pepi];
        VdotE_allbreaths=[VdotE_allbreaths; Breaths(1,i).VdotE];
        GGpeak_allbreaths=[GGpeak_allbreaths; Breaths(1,i).GGpeak];
        GGtonic_allbreaths=[GGtonic_allbreaths; Breaths(1,i).GGtonic];
    end
        GGcombineddata=[Pepi_allbreaths VdotE_allbreaths GGpeak_allbreaths GGtonic_allbreaths];
        % Graphing and recording all Data combined
        I=~isnan(VdotE_allbreaths)&~isnan(GGpeak_allbreaths);
        [fit.VdotEandGGpeak_combined]=linreg(VdotE_allbreaths(I),GGpeak_allbreaths(I),1,'VE(L/min)','GGpeak(%max)');
        I=~isnan(VdotE_allbreaths)&~isnan(GGtonic_allbreaths);
        [fit.VdotEandGGtonic_combined]=linreg(VdotE_allbreaths(I),GGtonic_allbreaths (I),2,'VE(L/min)','GGtonic(%max)');
        I=~isnan(Pepi_allbreaths)&~isnan(GGtonic_allbreaths);
        [fit.PepiandGGtonic_combined]=linreg(Pepi_allbreaths(I),GGtonic_allbreaths(I),4,'Pepi(cmH2O)','GGtonic(%max)');
        I=~isnan(Pepi_allbreaths)&~isnan(GGpeak_allbreaths);
        [fit.PepiandGGpeak_combined]=linreg(Pepi_allbreaths(I),GGpeak_allbreaths(I),3,'Pepi(cmH2O)','GGpeak(%max)');
        I=~isnan(Pepi_allbreaths)&~isnan(VdotE_allbreaths);
        [fit.PepiandVdotE_combined]=linreg(Pepi_allbreaths(I),VdotE_allbreaths(I),5,'Pepi(cmH2O)','VE(L/min)');
        fh = gcf;
        scrsz = get(0,'ScreenSize');
        set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
        saveas(fh, ['GGanalysis_combined']);
        pause;
        close all;
        clear a*
    clear N Pepi_allbreaths VdotE_allbreaths GGpeak_allbreaths GGtonic_allbreaths
    save GGresults Breaths fit GGcombineddata
else
    warning('Scaled GG signal is empty')
end

%% Measuring Vactive with all CPAP data
%loads files needed
clear all 
close all
load('Vactive.mat')
load('Var.mat')
load('Veupnea.mat')
temptitle=nanmean(Var.AvgCPAP);
% Discards Ventilations higher than Veupnea
tempVar=nanmean(Var.VdotE(Var.VdotE<(Veupnea_meandata(1))));
tempVactive=nanmean(Vactive.VdotE);
tempVdrive=Vactive.Vdrive+(Veupnea_meandata(1));
mean(tempVdrive);

hFig=figure(3);
set(hFig, 'Position', [200 400 1250 350])
ax1(1)=subplot(1,4,2); plot(Vactive.Cpap, Vactive.VdotE, 'b.', 'MarkerSize',20)
ylabel('VdotE (L/min)')
xlabel('CPAP pre-active drop(cmH2O)')

ax1(2)=subplot(1,4,1); plot(Var.AvgCPAP, Var.VdotE, 'g.', 'MarkerSize',20)
ylabel('VdotE (L/min)')
xlabel('Varousal CPAP level(cmH2O)')
title({['Varousal CPAP level: ', num2str(sprintf('%.1f',temptitle))]})

ax1(3)=subplot(1,4,3); plot(tempVdrive, Vactive.VdotE, 'r.', 'MarkerSize',20)
ylabel('VdotE (L/min)')
xlabel('Vdrive (L/min)')

ax1(4)=subplot(1,4,4); plot( Vactive.Cpap, tempVdrive, 'k.', 'MarkerSize',20)
ylabel('Vdrive (L/min)')
xlabel('CPAP (cmH2O)')
saveas(figure(3), 'Comparing_Vactive_with_Varousal');

%% Calculation SD and COV
%Load data (only if doing retrospectively)
clear all
files = dir('*.mat'); for i = 1:numel(files) load(files(i).name); end

%Calculate SDs and COVs
VeupneaSD=nanstd(Veupnea.VdotE);
VeupneaCOV=VeupneaSD/(nanmean(Veupnea.VdotE));

VarousalSD=nanstd(Var.VdotE);
VarousalCOV=VarousalSD/(nanmean(Var.VdotE));

VpassiveSD=nanstd(V0.VdotE);
VpassiveCOV=VpassiveSD/(nanmean(V0.VdotE));

VactiveSD=nanstd(Vactive.VdotE(minInd));
VactiveCOV=VactiveSD/(nanmean(Vactive.VdotE(minInd)));

LGSD=nanstd(LG_1br.LGwellman);
LGCOV=LGSD/(nanmean(LG_1br.LGwellman));

%Report data for export
night
SDandCOVdata=[VeupneaSD VeupneaCOV VarousalSD VarousalCOV VpassiveSD VpassiveCOV VactiveSD VactiveCOV LGSD LGCOV]

%Save variables to patient folder
save SDandCOVdata SDandCOVdata
