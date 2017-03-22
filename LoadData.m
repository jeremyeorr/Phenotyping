function [Data, Fs, Time] = LoadData(n)
%Written by Lisa Campana, PhD
%December 2012

%Loads MATLAB Data file exported by Spike
%Inputs: n - number of files to concatenate together (if data is split into more
%than 1 spike file)

%Outputs: 
%Data - Structure containing raw recorded data from Spike, SaO2,
%flow, tidal volume, O2, CO2, Pepi, Pmask, GG, Epochs, and Evt
%Fs - Structure containing the sampling rate of each channel
%Time - Structure containing time vector corresponding to flow channel.  If
%n>1 contains a variable Time.Break which lists the new start time of each
%successive file

if nargin<1
    n = 1; %set default to 1 file if not specified
end
LagTime = 30; %number of seconds to add between different files

%Create variables you want to extract
Data.SaO2 = []; %oxygen saturation
Data.Vdot = []; %flow
Data.Vt = []; %Tidal Volume
Data.PO2 = []; %O2
Data.PCO2 = []; %CO2
Data.Pepi = [];  %epiglottic pressure
Data.Pmask =[]; %mask pressure
Data.GG = []; % EMG GG
Data.Epochs = []; %sleep stages
Data.Evts=[]; %events
Data.Traits1Evts=[]; % 4 ventilations and LG events
Data.Traits2Evts=[]; % Pcrit and GG responsiveness events

% Load Files
for i = 1:n
    [fname, pathname, Findex] = uigetfile('*.mat', 'Select File'); % load file
    cd(pathname) 
    name = strcat(pathname,fname);
    load(fname); % Loads file!
    s = whos('-regexp', '_Ch'); % determines variable names for channels
    a = strfind(s(1).name,'_Ch');
    fname2 = s(1).name(1:a);

    clear a
    if i==1
        for k = 1:length(s)
            eval(['nametemp = ', s(k).name, '.title;']) 
            switch nametemp        %compares nametemp to different variable names
                case {'SaO2', 'Sao2','sao2'}
                    eval(['Data.SaO2 = ', s(k).name, '.values;']);
                    eval(['Fs.SaO2 = ', s(k).name, '.interval;']);
                    Time.SaO2 = single([Fs.SaO2:Fs.SaO2:length(Data.SaO2)*Fs.SaO2]');
                case {'Flow', 'flow'}
                    eval(['Data.Vdot = ', s(k).name, '.values;']);
                    eval(['Fs.Vdot = ', s(k).name, '.interval;']);
                    Time.Vdot = single([Fs.Vdot:Fs.Vdot:length(Data.Vdot)*Fs.Vdot]');
                case {'Vt','VT'}
                    eval(['Data.Vt = ', s(k).name, '.values;']);
                    eval(['Fs.Vt = ', s(k).name, '.interval;']);
                case {'Pmask','PMask'}
                    eval(['Data.Pmask = ', s(k).name, '.values;']);
                    eval(['Fs.Pmask = ', s(k).name, '.interval;']);
                case {'PCO2','CO2 anal', 'pCO2'}
                    eval(['Data.PCO2 = ', s(k).name, '.values;']);
                    eval(['Fs.PCO2 = ', s(k).name, '.interval;']);
                case {'pO2','PO2', 'O2 anal'}
                    eval(['Data.PO2 = ', s(k).name, '.values;']);
                    eval(['Fs.PO2 = ', s(k).name, '.interval;']);
                case {'Pepi','pepi', 'Epi', 'epi'}
                    eval(['Data.Pepi = ', s(k).name, '.values;']);
                    eval(['Fs.Pepi = ', s(k).name, '.interval;']);
                case {'EMGgg ra', 'GG', 'gg'}
                    eval(['Data.GG = ', s(k).name, '.values;']);
                    eval(['Fs.GG = ', s(k).name, '.interval;']);
                case 'Epochs'
                    eval(['Epochs = ', s(k).name, '.codes;']);
                    Data.Epochs = single(Epochs(:,1));
                    Time.Epochs = single([0:30:(length(Data.Epochs)-1)*30]');
                case {'New Evts','Evts'}
                    eval(['Evts = ', s(k).name, '.codes;']);
                    eval(['Evttimes = ', s(k).name, '.times;']);
                    eval(['Evtname = ', s(k).name, '.text;']);
                    Data.Evts = single(Evts(1:2:end,1));
                    Data.EvtTimeStart = Evttimes(1:2:end);
                    Data.EvtTimeEnd= Evttimes(2:2:end);
                    Data.Evtname = Evtname(1:2:end,:);
                case 'Traits1'
                    eval(['Traits1Evts = ', s(k).name, '.codes;']);
                    eval(['Traits1times = ', s(k).name, '.times;']);
                    eval(['Traits1name = ', s(k).name, '.text;']);
                    Data.Traits1Evts = single(Traits1Evts(1:2:end,1));
                    Data.Traits1EvtTimeStart = Traits1times(1:2:end);
                    Data.Traits1EvtTimeEnd= Traits1times(2:2:end);
                    Data.Traits1Evtname = Traits1name(1:2:end,:);
                case 'Traits2'
                    eval(['Traits2Evts = ', s(k).name, '.codes;']);
                    eval(['Traits2times = ', s(k).name, '.times;']);
                    eval(['Traits2name = ', s(k).name, '.text;']);
                    Data.Traits2Evts = single(Traits2Evts(1:2:end,1));
                    Data.Traits2EvtTimeStart = Traits2times(1:2:end);
                    Data.Traits2EvtTimeEnd= Traits2times(2:2:end);
                    Data.Traits2Evtname = Traits2name(1:2:end,:);        
            end
        end
        tend = Time.Vdot(end); %last time of file 1
        Time.Break = [];
        clear -regexp '_Ch'
    else %concatenates next file on end of first one
         for k = 1:length(s)
            eval(['nametemp = ', s(k).name, '.title;'])
            switch nametemp
                case {'SaO2', 'Sao2','sao2'}
                    eval(['tempSaO2 = ', s(k).name, '.values;']);
                    tempTime = [Fs.SaO2:Fs.SaO2:length(tempSaO2)*Fs.SaO2]'+tend+LagTime; %add last time of file1 + lag time to new time vector
                    Data.SaO2 = [Data.SaO2; tempSaO2];
                    Time.SaO2 = [Time.SaO2; tempTime];
                    
                case {'Flow', 'flow'}
                    eval(['tempVdot = ', s(k).name, '.values;']);
                    tempTime = [Fs.Vdot:Fs.Vdot:length(tempVdot)*Fs.Vdot]'+tend+LagTime;
                    Data.Vdot = [Data.Vdot; tempVdot];
                    Time.Vdot = [Time.Vdot; tempTime];
                    
                case {'Vt','VT'}
                    eval(['tempVt = ', s(k).name, '.values;']);
                    Data.Vt = [Data.Vt; tempVt];%                     
                    
                case {'Pmask','PMask'}
                    eval(['tempPmask = ', s(k).name, '.values;']);
                    Data.Pmask = [Data.Pmask; tempPmask];
                  
                case {'PCO2','CO2 anal'}
                    eval(['tempPCO2 = ', s(k).name, '.values;']);
                    Data.PCO2 = [Data.PCO2; tempPCO2];
                  
                case {'pO2','PO2', 'O2 anal'}
                    eval(['tempPO2 = ', s(k).name, '.values;']);
                    Data.PO2 = [Data.PO2; tempPO2];
                 
                case {'Pepi','pepi', 'Epi', 'epi'}
                    eval(['tempPepi = ', s(k).name, '.values;']);
                    Data.Pepi = [Data.Pepi; tempPepi];
                  
                case {'EMGgg ra', 'GG', 'gg'}
                    eval(['tempGG = ', s(k).name, '.values;']);
                    Data.GG = [Data.GG; tempGG];
                  
                case 'Epochs'
                    eval(['Epochs = ', s(k).name, '.codes;']);
                    tempEpochs = single(Epochs(:,1));
                    Data.Epochs = [Data.Epochs; tempEpochs];
                    tempTime = [0:30:(length(tempEpochs)-1)*30]'+tend+LagTime;
                    Time.Epochs = [Time.Epochs; tempTime];
                    
                case {'New Evts','Evts'}
                    eval(['Evts = ', s(k).name, '.codes;']);
                    eval(['tempTimes = ', s(k).name, '.times;']);
                    eval(['tempname = ', s(k).name, '.text;']);
                    
                    tempEvts = single(Evts(1:2:end,1));
                    tempTimes = tempTimes+tend+LagTime;
                    Data.Evts = [Data.Evts;tempEvts];
                    Data.EvtTimeStart = [Data.EvtTimeStart; tempTimes(1:2:end)];
                    Data.EvtTimeEnd= [Data.EvtTimeEnd; tempTimes(2:2:end)];
                    Data.Evtname = [Data.Evtname; tempname(1:2:end,:)];
               
                case 'Traits1'
                    eval(['Traits1Evts = ', s(k).name, '.codes;']);
                    eval(['tempTraits1times = ', s(k).name, '.times;']);
                    eval(['tempTraits1name = ', s(k).name, '.text;']);
                    
                    tempTraits1 = single(Traits1Evts(1:2:end,1));
                    tempTimes1 = tempTraits1times+tend+LagTime;
                    Data.Traits1Evts = [Data.Traits1Evts;tempTraits1];
                    Data.Traits1EvtTimeStart = [Data.Traits1EvtTimeStart; tempTimes1(1:2:end)];
                    Data.Traits1EvtTimeEnd= [Data.Traits1EvtTimeEnd; tempTimes1(2:2:end)];
                    Data.Traits1Evtname = [Data.Traits1Evtname; tempTraits1name(1:2:end,:)];
                    
                case 'Traits2'
                    eval(['Traits2Evts = ', s(k).name, '.codes;']);
                    eval(['tempTraits2times = ', s(k).name, '.times;']);
                    eval(['tempTraits2name = ', s(k).name, '.text;']);
                    
                    tempTraits2 = single(Traits2Evts(1:2:end,1));
                    tempTimes2 = tempTraits2times+tend+LagTime;
                    Data.Traits2Evts = [Data.Traits2Evts;tempTraits1];
                    Data.Traits2EvtTimeStart = [Data.Traits2EvtTimeStart; tempTimes2(1:2:end)];
                    Data.Traits2EvtTimeEnd= [Data.Traits2EvtTimeEnd; tempTimes2(2:2:end)];
                    Data.Traits2Evtname = [Data.Traits2Evtname; tempTraits2name(1:2:end,:)];  
            end
         end
         Time.Break(i-1) = tend+LagTime;
        tend = Time.Vdot(end);
    end
    clear temp* 
    clear -regexp '_Ch'
end
%% Warn user if any variables are empty
if isempty(Data.SaO2)
    warning('Empty SaO2')
end
if isempty(Data.Vdot)
    warning('Empty Flow')
end
if isempty(Data.Pmask)
    warning('Empty Pmask')
end
if isempty(Data.Pepi)
    warning('Empty Pepi')
end
if isempty(Data.PCO2)
    warning('Empty PCO2')
end
if isempty(Data.PO2)
    warning('Empty PO2')
end
if isempty(Data.GG)
    warning('Empty GG')
end
if isempty(Data.Epochs)
    warning('Empty Epoch')
end
if isempty(Data.Evts)
    warning('Empty Evt')
end

clear (strcat(fname2,'*'))

