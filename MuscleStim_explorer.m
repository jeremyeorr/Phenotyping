%% Muscle stimulation study analysis
% created Brad Edwards, PhD (August 2014)
% Loads in data file with labeled times you want to analyze  

%% Add paths, load spike data and save as Matlab file for analysis
clear all
close all
clc
% Add paths
addpath('C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode')

%Change RouteDir to the location of IdentifyTimePeriodsofInterest.m 
RouteDir = 'C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode';
addpath(RouteDir)

%Specify Directories where you want to load and then save data
%SaveDir = ['C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode\Patient_data']; 
SaveDir = ['C:\Users\bae6\Documents\MATLAB\Muscle Stimulation']; 
%LoadDir = ['J:\PEOPLE\FACULTY\EDWARDS\Brad_work\Data\Combination therapy\Trait nights\Scored files']; 
LoadDir = ['J:\PEOPLE\FACULTY\EDWARDS\Brad_work\Data\Muscle stimulation study\Scored data'];

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
addpath('C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode')
%Loadsaveddata = ['C:\Users\bae6\Documents\MATLAB\Combination therapy\Matlab files\Edwards_newcode\Patient_data']; 
Loadsaveddata = ['C:\Users\bae6\Documents\MATLAB\Muscle Stimulation']; 
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

%% Awake analysis
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
    [Veupnea.VdotE(i), Veupnea.Cpap(i), Veupnea.Pepi(i), Veupnea.PO2(i), Veupnea.PCO2(i),Veupnea.GGphasic(i), Veupnea.GGtonic(i), Veupnea.GGpeak(i),Veupnea.Ti(i),Veupnea.Te(i),Veupnea.Ttot(i),Veupnea.VT(i)]= Veup(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, pO2temp, pCO2temp, timeGG, tempGG);
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


%% Pcrit analysis
indxPcrit=find(Data.Traits2Evts==7);
DetectThreshold=0.5; %Lisa had set to 0.75
Pcrit = struct('Cpap', [], 'VdotMax', [], 'VdotE', [], 'Time', [], 'pc', []);
for i = 1:length(indxPcrit)
        a = find(Time.Vdot>=Data.Traits2EvtTimeStart(indxPcrit(i),1) & Time.Vdot<=Data.Traits2EvtTimeEnd(indxPcrit(i),1));
        a2 = find(Time.Epochs>=Data.Traits2EvtTimeStart(indxPcrit(i),1)-30 & Time.Epochs<=Data.Traits2EvtTimeEnd(indxPcrit(i),1));
        a3 = find(Data.EvtTimeEnd>=Data.Traits2EvtTimeStart(indxPcrit(i),1) & Data.EvtTimeStart<=Data.Traits2EvtTimeEnd(indxPcrit(i),1));
        [Pcrit(i).Cpap, Pcrit(i).VdotMax, Pcrit(i).VdotE, Pcrit(i).Time, Pcrit(i).pc Vcrit(i).vc] = PcritCalc(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], DetectThreshold);
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

%% Modified Pcrit (breakpoint analysis to include zero mask pressure)
global Pmaskdata VdotMaxdata VdotEdata fixedslope upperlimit
close all
Pmaskdata=[];
VdotMaxdata=[];
VdotEdata=[];
fixedslope=NaN
upperlimit.on=0
exclude_zeroflowdata=0 % zero idicates include zero data points
includevpassivedata=1 %1 indicates you want to add in the vpassive data
plotpoints=1

include_runs=[1 2 3] %manually select which runs to include. Eg. [1]-will graph run 1 [1 2]-will graph runs 1 and 2
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


