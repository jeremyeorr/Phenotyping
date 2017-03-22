function [AvgVdotE, AvgCpap, DelPepi, AvgPO2, AvgPCO2, GGphasic, GGtonic, GGpeak, AvgTi, AvgTe, AvgTtot, AvgVT]= Veup(Fs, Time, Vdot, Pmask, stage, evts, Pepi, PO2, PCO2, TimeGG, GG)
% Script for getting Eupneic Ventilation on optimum CPAP
%Lisa M Campana, PhD
%December 2012

% Integrate the raw flow signal to get volume
vol =cumtrapz(Vdot*Fs);
%detrend volume signal to remove leak
vol = detrend(vol, 'linear');
%Calculate end inspiration and end expiration points, tidal volumes
[TidalVol,indVI,indVE] = peakVol(Time,vol,0.1);

% find the CPAP by looking at the Pmask at the start of inspiration.
cpap = Pmask(indVE);
cpap(end) = [];

if isempty(PCO2) & isempty(GG)
    numplots = 6;
elseif isempty(GG)
    numplots = 8;
elseif isempty(PCO2);
    numplots = 7;
else
    numplots=9;
end
% Plot figure
fh = figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(numplots,1,1); %Sleep stage
hold on;
w = find(stage(:,2) == 0);
s1 = find(stage(:,2) == 1);
s2 = find(stage(:,2) == 2);
s3 = find(stage(:,2) == 3);
rem = find(stage(:,2) ==5);
for i = 1:length(w)
    plot([stage(w(i),1), stage(w(i),1)+30], [stage(w(i),2), stage(w(i),2)],'k-', 'LineWidth', 5)
end
for i = 1:length(s1)
    plot([stage(s1(i),1), stage(s1(i),1)+30], [stage(s1(i),2), stage(s1(i),2)],'b-', 'LineWidth', 5)
end
for i = 1:length(s2)
    plot([stage(s2(i),1), stage(s2(i),1)+30], [stage(s2(i),2), stage(s2(i),2)],'g-', 'LineWidth', 5)
end
for i = 1:length(s3)
    plot([stage(s3(i),1), stage(s3(i),1)+30], [stage(s3(i),2), stage(s3(i),2)],'y-', 'LineWidth', 5)
end
for i = 1:length(rem)
    plot([stage(rem(i)), stage(rem(i),1)+30], [stage(rem(i),2), stage(rem(i),2)],'r-', 'LineWidth', 5)
end
for i = 1:size(evts,1)
    plot([evts(i,1), evts(i,2)],[6,6], 'm', 'LineWidth', 5);
end
ylabel('Sleep Stage')
xlim([Time(1) Time(end)])
ylim([0 6])

ax(2) = subplot(numplots,1,2); % Pmask
plot(Time, Pmask,'k')
hold on;
plot(Time(indVE(1:end-1)), cpap, 'r.-')
ylabel('Pmask (cmH2O')

ax(3) = subplot(numplots,1,3);%Pepi (if it exists)
if ~isempty(Pepi)
    plot(Time, Pepi,'k')
end
hold on
ylabel('Pepi (cmH2O)')

ax(4) = subplot(numplots,1,4); %Volume
plot(Time, vol, 'k')
hold on;
hI = plot(Time(indVI), vol(indVI),'r+');
hE = plot(Time(indVE), vol(indVE),'ro');
ylabel('Volume (L)')

ax(5) = subplot(numplots,1,5); %Flow
plot(Time, Vdot,'k')
hold on
plot([Time(1), Time(end)], [0,0],'r')
ylabel('Flow (L/s)')

if numplots>7
   ax(6) = subplot(numplots,1,6);
   plot(Time, PO2,'k')
   ylabel('PO2')
   axis([Time(1) Time(end) min(PO2)-1 max(PO2)]+1)
   hold on;
   ax(7) = subplot(numplots,1,7); 
   plot(Time, PCO2, 'k')
   hold on
      ylabel('PCO2')
end
if ~isempty(TimeGG)
    if numplots==7
        ggplot = 6;
    elseif numplots==9
        ggplot = 8;
    end

    ax(ggplot) = subplot(numplots,1,ggplot);
    plot(TimeGG, GG,'k'); hold on;
    ylabel('GG');
    hold on;
end

linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
% set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])

% Fix any volume points?
RedoVol = listdlg('PromptString','Fix any volumes?',...
    'SelectionMode','single', 'ListSize', [150 50],...
    'ListString',{'Yes','No'});

if RedoVol==1
    [indVI, indVE] = FixVolumes(Time, vol, indVI, indVE);
    TidalVol = vol(indVI)-vol(indVE);
    set(hI, 'Visible', 'off')
    set(hE, 'Visible', 'off')
    subplot(numplots,1,4)
    plot(Time(indVI), vol(indVI),'r+');
    plot(Time(indVE), vol(indVE),'ro');
end

% Determine Minute ventilation
VdotE = [];
for n = 1:length(indVI)-1;
    Ttot(n) = Time(indVE(n+1))-Time(indVE(n));
    Ti(n) = Time(indVI(n))-Time(indVE(n));
    Te(n) = Ttot(n)-Ti(n);
    VdotE(n) = (TidalVol(n))/(Ttot(n));
    
end
VdotE = VdotE.*60;
a = find(VdotE<0);
VdotE(a) = 0;
VdotE = VdotE';
clear n;

AvgCpap = mean(cpap);
AvgVdotE = mean(VdotE);
AvgVT = mean(TidalVol);
AvgTi = mean(Ti);
AvgTe = mean(Te);
AvgTtot = mean(Ttot);

%Plot Minute ventilation
figure(fh)
ax(numplots) = subplot(numplots,1,numplots);
plot(Time(indVI(1:end-1)), VdotE,'k.-')
hold on
plot([Time(indVI(1)), Time(indVI(end))], [AvgVdotE,AvgVdotE],'r')
text(double(Time(indVI(round(length(indVI)/2))))-1,round(AvgVdotE)+1, strcat('VdotE mean =', num2str(AvgVdotE,3)))
ylabel('VdotE (L/min)')

linkaxes(ax,'x')

% Find peaks/valleys of Pepi signal
if ~isempty(Pepi)
    %low pass filter pepi signal, cutoff 2hz
    cutoff=2;
    wn = cutoff/(1/(2*Fs)); 
    [den,num] = butter(2, wn);
    pfilt = filtfilt(den,num,Pepi); 
    
    cutoff=5;
    wn = cutoff/(1/(2*Fs)); 
    [den,num] = butter(2, wn);
    pfilt2 = filtfilt(den,num,Pepi);
   
    %compute 2nd derivative of smoothed pepi signal
    yp = diff(pfilt)/Fs;
    ypp = diff(yp)/Fs;
    ypp(end+1:end+2) = 0; 
    yp(end+1) = 0; 
    
    %find max and mins of smoothed pepi
    [mx,mn] = peakdet(pfilt, 1);
    if length(mx(:,1))<0.75*length(indVI)
        [mx,mn] = peakdet(pfilt, .5);
    end
    if mx(1,1)>mn(1,1); 
        mn(1,:) = [];
    end
    if length(mx(:,1))>length(mn(:,1))
        mx(end,:) = [];
    end
    
    %repick the peaks of pepi based on where the minimum of the 2nd derivative
    for i = 1:length(mx(:,1))
        [peakpepi(i), peakind(i)] = min(ypp(mx(i,1):mn(i,1)));
        peakind(i) = peakind(i)+mx(i,1)-1;
    end
    subplot(numplots,1,3)
    plot(Time(mn(:,1)), pfilt2(mn(:,1)),'r+')
    plot(Time(peakind), pfilt2(peakind),'ro' )
    DelPepi = mean(pfilt2(peakind)-pfilt2(mn(:,1)));
    
    KeepPepi = listdlg('PromptString','Pepi quality',...
                       'SelectionMode','single', 'ListSize', [150 50],...
                       'ListString',{'Keep','Do not use','Repick'});
    if KeepPepi==2
        DelPepi = nan;
    end
    if KeepPepi==3
        [mn, peakind] = FixPepi(Time, pfilt2, Vdot, mn, peakind);
        peakind'
        DelPepi = mean(pfilt2(peakind)-pfilt2(mn(:,1)));
        figure(fh)
        subplot(numplots,1,3)
        plot(Time(mn(:,1)), pfilt2(mn(:,1)),'g+');
        plot(Time(peakind), pfilt2(peakind),'go' );
    end
else
    DelPepi = nan;
end

figure(fh)
%Find miminums of P02 signal
if ~isempty(PO2)
    %filter PO2 to remove any ECG artifcat
    [b,a] = butter(2, 2/(125/2),'low'); %low pass butterworth filter, order 2, cutoff frequency = 2hz;
    pof = filtfilt(b,a,PO2); 
    subplot(numplots,1,6)
    plot(Time, pof,'r--');
    for i = 1:length(indVI)-1
        [mnPO2(i), mnIndPO2(i)] = min(pof(indVI(i):indVE(i+1))); %find minimum PO2 between end inspiration and start of next breath 
        mnIndPO2(i) = mnIndPO2(i)+indVI(i)-1;
    end
    a = find(PO2(mnIndPO2)<50); 
    PO2(mnIndPO2(a)) = nan; 
    clear a
    plot(Time(mnIndPO2), pof(mnIndPO2),'g+');
    AvgPO2 = mean(PO2(mnIndPO2));
else
    AvgPO2 = nan;
end

%Find maximums of PCO2
if ~isempty(PCO2)
    for i = 1:length(indVI)-1
        [mxPCO2(i), mxIndPCO2(i)] = max(PCO2(indVI(i):indVE(i+1))); %find minimum PO2 between end inspiration and start of next breath 
        mxIndPCO2(i) = mxIndPCO2(i)+indVI(i)-1;
    end
    subplot(numplots,1,7)
    plot(Time(mxIndPCO2), PCO2(mxIndPCO2),'g+');
    a = find(PCO2(mxIndPCO2)<25); 
    PCO2(mxIndPCO2(a)) = nan; 
    clear a
    AvgPCO2 = nanmean(PCO2(mxIndPCO2));    
else
    AvgPCO2 = nan;
end

%GG analysis
if ~isempty(GG)
    for i = 2:length(indVI)
        indMn = find(TimeGG>=Time(indVI(i-1)) & TimeGG<=Time(indVE(i)));
        indMx = find(TimeGG>=Time(indVE(i)) & TimeGG<=Time(indVI(i)));
      [mnGG(i-1), mnIndGG(i-1)] = min(GG(indMn));  
      [mxGG(i-1), mxIndGG(i-1)] = max(GG(indMx));
       mxIndGG(i-1) = mxIndGG(i-1) + indMx(1)-1;
       mnIndGG(i-1) = mnIndGG(i-1) + indMn(1)-1;
    end
      
    GGtonic = mean(mnGG); 
    GGpeak = mean(mxGG);
    GGphasic = mean(mxGG-mnGG); 
    
    subplot(numplots,1,ggplot);
    plot(TimeGG(mxIndGG), mxGG,'r+'); 
    plot(TimeGG(mnIndGG), mnGG, 'ro');
    plot([TimeGG(mnIndGG(1)), TimeGG(mnIndGG(end))], [GGtonic, GGtonic],'g')
    plot([TimeGG(mxIndGG(1)), TimeGG(mxIndGG(end))], [GGphasic, GGphasic],'b')
    plot([TimeGG(mxIndGG(1)), TimeGG(mxIndGG(end))], [GGpeak, GGpeak],'r')
else
    GGphasic = nan;
    GGtonic = nan;
    GGpeak = nan;
end