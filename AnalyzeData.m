%Lisa M. Campana, PhD
%December 2012

%Loads in data file with labeled times you want to analyze.  
%Calls functions Veup.m, Varousal.m, Vzero.m, Uag.m, LoopGain.m,
%PcritCalc.m, GGanalysis.m
%% Add paths, and load data file
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

%% Must have run CalandScaleGG.m before loading here
%load scaled and calilbrated GG signal
load GGnew  
%% Veupnea analyze

% Find all locations were you want to analyze Veupnea
%indxVeup = find(strcmp(pick, choices{1})); %
indxVeup=find(Traits1.codes==1); %

%Loop over all windows of interest
for i = 1:length(indxVeup);
    a = find(Time.Vdot>=x(indxVeup(i),1) & Time.Vdot<=x(indxVeup(i),2));
    a2 = find(Time.Epochs>=x(indxVeup(i),1)-30 & Time.Epochs<=x(indxVeup(i),2));
    a3 = find(Data.EvtTimeEnd>=x(indxVeup(i),1) & Data.EvtTimeStart<=x(indxVeup(i),2));
    
    %create temp GG signal over time period of interset
    if exist('GGscaled')
         tempGG = GGscaled(ceil(x(indxVeup(i),1)/Fs.GG): floor(x(indxVeup(i),2)/Fs.GG));
         timeGG = [ceil(x(indxVeup(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(x(indxVeup(i),2)/Fs.GG)*Fs.GG]';
    else
        warning('Scaled GG signal is empty')
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
    Veupnea.Time(i,:) = x(indxVeup(i),:);

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

%save variable
save Veupnea Veupnea night

%% Varousal analyze
% Find all locations were you want to analyze Varousal
indxVar = find(strcmp(pick, choices{3}));

%Loop over all windows of interest
for i = 1:length(indxVar);
    a = find(Time.Vdot>=x(indxVar(i),1) & Time.Vdot<=x(indxVar(i),2));
    a2 = find(Time.Epochs>=x(indxVar(i),1)-30 & Time.Epochs<=x(indxVar(i),2));
    a3 = find(Data.EvtTimeEnd>=x(indxVar(i),1) & Data.EvtTimeStart<=x(indxVar(i),2));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
     if exist('GGscaled')
         tempGG = GGscaled(ceil(x(indxVar(i),1)/Fs.GG): floor(x(indxVar(i),2)/Fs.GG));
         timeGG = [ceil(x(indxVar(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(x(indxVar(i),2)/Fs.GG)*Fs.GG]';
    else
        warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    [Var.VdotE(i), Var.Pepi(i), Var.GGphasic(i), Var.GGtonic(i), Var.GGpeak(i), Var.Ti(i), Var.Te(i), Var.Ttot(i), Var.VT(i), Var.AvgCPAP(i), Var.tarousal(i)]= Varousal(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
    Var.Time(i,:) = x(indxVar(i),:);
    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVar))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp
end

save Var Var 


%% Vpassive analyze

% Find all locations were you want to analyze Vpassive
indxV0 = find(strcmp(pick, choices{4}));

%Loop over all windows of interest
for i = 1:length(indxV0);
    a = find(Time.Vdot>=x(indxV0(i),1) & Time.Vdot<=x(indxV0(i),2));
    a2 = find(Time.Epochs>=x(indxV0(i),1)-30 & Time.Epochs<=x(indxV0(i),2));
    a3 = find(Data.EvtTimeEnd>=x(indxV0(i),1) & Data.EvtTimeStart<=x(indxV0(i),2));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    if exist('GGscaled')
         tempGG = GGscaled(ceil(x(indxV0(i),1)/Fs.GG): floor(x(indxV0(i),2)/Fs.GG));
         timeGG = [ceil(x(indxV0(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(x(indxV0(i),2)/Fs.GG)*Fs.GG]';
    else
        warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    [V0.VdotE(i), V0.Pepi(i), V0.GGphasic(i), V0.GGtonic(i), V0.GGpeak(i), V0.Ti(i), V0.Te(i), V0.Ttot(i), V0.VT(i)]= Vzero(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
    V0.Time(i,:) = x(indxV0(i),:);
    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxV0))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp
end

save V0 V0 

%% LG analyze

indxLG = find(strcmp(pick, choices{2}));

for i = 1:length(indxLG);
    a = find(Time.Vdot>=x(indxLG(i),1) & Time.Vdot<=x(indxLG(i),2));
    a2 = find(Time.Epochs>=x(indxLG(i),1)-30 & Time.Epochs<=x(indxLG(i),2));
    a3 = find(Data.EvtTimeEnd>=x(indxLG(i),1) & Data.EvtTimeStart<=x(indxLG(i),2));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
    [LG.Lisa(i), LG.Andrew(i)] = LoopGain(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, mean(Veupnea.VdotE));
    LG.Time(i,:) = x(indxLG(i),:);

    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxLG))]; 'Press Enter to continue'})
    pause;
    close all
    clear a*
end

save LG LG

%% Vactive analyze
indxVact = find(strcmp(pick, choices{5}));
% Vactive = [];
for i = 1:length(indxVact);
    a = find(Time.Vdot>=x(indxVact(i),1) & Time.Vdot<=x(indxVact(i),2));
    a2 = find(Time.Epochs>=x(indxVact(i),1)-30 & Time.Epochs<=x(indxVact(i),2));
    a3 = find(Data.EvtTimeEnd>=x(indxVact(i),1) & Data.EvtTimeStart<=x(indxVact(i),2));
    if isempty(Data.Pepi)
        peptemp = [];
    else
        peptemp= Data.Pepi(a);
    end
     if exist('GGscaled')
         tempGG = GGscaled(ceil(x(indxVact(i),1)/Fs.GG): floor(x(indxVact(i),2)/Fs.GG));
         timeGG = [ceil(x(indxVact(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(x(indxVact(i),2)/Fs.GG)*Fs.GG]';
    else
        warning('Scaled GG signal is empty')
        tempGG = [];
        timeGG = [];
    end
    [Vactive.VdotE(i), Vactive.Vdrive(i), Vactive.Ug(i), Vactive.Cpap(i), Vactive.Pepi(i), Vactive.GGphasic(i),Vactive.GGtonic(i), Vactive.GGpeak(i), Vactive.Ti(i), Vactive.Te(i), Vactive.Ttot(i), Vactive.VT(i)]=...
    Uag(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG, nanmean(Veupnea.VdotE), nanmean(LG.Andrew), nanmean(V0.VdotE));
    
    Vactive.Time(i,:) = x(indxVact(i),:);

    subplot(6,1,1)
    title({['Number ', num2str(i), ' of ', num2str(length(indxVact))]; 'Press Enter to continue'})
    pause;
    close all
    clear a* peptemp

end

save Vactive Vactive 

%% Summarize and put data together
subid = inputdlg('Enter subject Id');
SummarizeData(subid{1}, night, mean(Veupnea.VdotE), nanmean(Var.VdotE), Vactive.VdotE, nanmean(V0.VdotE), Vactive.Ug, -nanmean(LG.Andrew), Vactive.Cpap, nanmean(Var.Pepi))

%% Pcrit analysis

indxPcrit = find(strcmp(pick, choices{8}));
Pcrit = struct('Cpap', [], 'VdotMax', [], 'VdotE', [], 'Time', [], 'pc', []);
for i = 1:length(indxPcrit)
        a = find(Time.Vdot>=x(indxPcrit(i),1) & Time.Vdot<=x(indxPcrit(i),2));
        a2 = find(Time.Epochs>=x(indxPcrit(i),1)-30 & Time.Epochs<=x(indxPcrit(i),2));
        a3 = find(Data.EvtTimeEnd>=x(indxPcrit(i),1) & Data.EvtTimeStart<=x(indxPcrit(i),2));
        [Pcrit(i).Cpap, Pcrit(i).VdotMax, Pcrit(i).VdotE, Pcrit(i).Time, Pcrit(i).pc] = PcritCalc(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)]);
        pause
        close all
        clear a
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
    Vcrit=p(2);

    title({'VdotE average' ; ['VdotE at 0 pressure: ', num2str(p(2))]})
    saveas(fh, 'VdotEcrit');
    save Pcrit Pcrit VdotEall Pall Vcrit
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
    save Pcrit Pcrit Pall2 VdotMaxall PcritAll -append
end

%% Modified Pcrit (breakpoint analysis to include zero mask pressure)
global Pmaskdata VdotMaxdata VdotEdata

Pmaskdata=[];
VdotMaxdata=[];
VdotEdata=[];

include_runs=[1 2 3 4] %manually select which runs to include. Eg. [1]-will graph run 1 [1 2]-will graph runs 1 and 2

for i=include_runs
tempsize=size(Pcrit(i).VdotMax);
N=tempsize(1)*tempsize(2);

Pmaskdata=[Pmaskdata reshape(Pcrit(i).Cpap,1,N)]; 
VdotMaxdata=[VdotMaxdata reshape(Pcrit(i).VdotMax,1,N)];
VdotEdata=[VdotEdata reshape(Pcrit(i).VdotE,1,N)];
end 
    
VdotMaxdata(Pmaskdata==0)=[];
VdotEdata(Pmaskdata==0)=[];
Pmaskdata(Pmaskdata==0)=[];

pcrit_model_run()

%% GG analysis
indxGG = find(strcmp(pick, choices{6}));
if exist('GGscaled')
    for i = 1:length(indxGG)
        a = find(Time.Vdot>=x(indxGG(i),1) & Time.Vdot<=x(indxGG(i),2));
        a2 = find(Time.Epochs>=x(indxGG(i),1)-30 & Time.Epochs<=x(indxGG(i),2));
        a3 = find(Data.EvtTimeEnd>=x(indxGG(i),1) & Data.EvtTimeStart<=x(indxGG(i),2));
        if isempty(Data.Pepi)
            peptemp = [];
        else
            peptemp= Data.Pepi(a);
        end
        tempGG = GGscaled(ceil(x(indxGG(i),1)/Fs.GG): floor(x(indxGG(i),2)/Fs.GG));
        timeGG = [ceil(x(indxGG(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(x(indxGG(i),2)/Fs.GG)*Fs.GG]';
        [Breaths(i)] = GGanalysis(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
         subplot(7,1,1)
         title({['Number ', num2str(i), ' of ', num2str(length(indxGG))]; 'Press Enter to continue'})
          pause;
         close all
         
         %Calculating GGphasic stuff
        [fit.VdotEandGGphasic(i,:)]=linreg(Breaths(i).VdotE,Breaths(i).GGphasic,1,'VE(L/min)','GGphasic(%max)');
        [fit.VdotEandGGpeak(i,:)]=linreg(Breaths(i).VdotE,Breaths(i).GGpeak,2,'VE(L/min)','GGpeak(%max)');    
        [fit.VdotEandGGtonic(i,:)]=linreg(Breaths(i).VdotE,Breaths(i).GGtonic,3,'VE(L/min)','GGtonic(%max)'); 
        [fit.PepiandGGtonic(i,:)]=linreg(Breaths(i).Pepi,Breaths(i).GGtonic,6,'Pepi(cmH2O)','GGtonic(%max)'); 
        [fit.PepiandGGphasic(i,:)]=linreg(Breaths(i).Pepi,Breaths(i).GGphasic,4,'Pepi(cmH2O)','GGphasic(%max)'); 
        [fit.PepiandGGpeak(i,:)]=linreg(Breaths(i).Pepi,Breaths(i).GGpeak,5,'Pepi(cmH2O)','GGpeak(%max)'); 
        [fit.PepiandVdotE(i,:)]=linreg(Breaths(i).Pepi,Breaths(i).VdotE,7,'Pepi(cmH2O)','VE(L/min)'); 
        fh = gcf;
        scrsz = get(0,'ScreenSize');
        set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
        saveas(fh, ['GGanalysis', num2str(round(Time.Vdot(a(1))))]);      
         pause;
         close all;
            clear a*
    end
    
    save GGresults Breaths fit
else
    warning('Scaled GG signal is empty')

end
%% Write to excel
WriteResultsToExcel(Veupnea, Var, V0, LG, Vactive, Pcrit, Breaths, fit)