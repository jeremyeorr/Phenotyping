function [var, DelPepi, GGphasic, GGtonic, GGpeak, Ti_ar, Te_ar, Ttot_ar, VT_ar, AvgCPAP, tarousal] = Varousal(Fs, Time, Vdot, Pmask, stage, evts, Pepi, TimeGG, GG)
% Script for calculating ventilation before arousal from sleep
%Lisa M Campana, PhD
%December 2012

close all
% Integrate the raw flow signal to get volume
vol =cumtrapz(Vdot*Fs);
%detrend volume signal to remove leak
vol = detrend(vol, 'linear');
%Calculate end inspiration and end expiration points, tidal volumes
[TidalVol,indVI,indVE] = peakVol(Time,vol,0.1);

if isempty(GG)
    numplots = 6;
else
    numplots = 7;
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
ylabel('Pmask (cmH2O')

ax(3) = subplot(numplots,1,3);%Pepi (if it exists)
if ~isempty(Pepi)
    plot(Time, Pepi,'k')
    axis([Time(1) Time(end) min(Pepi) max(Pepi)])
end
hold on
ylabel('Pepi (cmH2O)')

ax(4) = subplot(numplots,1,4); %Volume
v1temp = plot(Time, vol, 'k');
hold on;
hI = plot(Time(indVI), vol(indVI),'r+');
hE = plot(Time(indVE), vol(indVE),'ro');
ylabel('Volume (L)')

ax(5) = subplot(numplots,1,5); %Flow
plot(Time, Vdot,'k')
hold on
plot([Time(1), Time(end)], [0,0],'r')
ylabel('Flow (L/s)')

ax(6) = subplot(numplots,1,6);
ylabel('VdotE (L/min')

if numplots==7
    ax(7) = subplot(numplots,1,7);
    plot(TimeGG, GG,'k');
    ylabel('GG')
    hold on;
end
linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
% set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])


%Fix any volume points??
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
if RedoVol == 0
    return
end
% Determine vdote
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

%Pick where arousal occurs
figure(fh);
subplot(numplots,1,1)
title('LEFT CLICK ON the arousal breath');
[tarousal,y] = ginput(1);
aind = find(Time(indVI)<tarousal,1,'last'); % last inspiratory index before arousal breath

%plot arousal location
subplot(numplots,1,1)
plot([tarousal,tarousal], [0 6],'r')

subplot(numplots,1,2)
plot([tarousal,tarousal], [min(Pmask) max(Pmask)],'r')

if ~isempty(Pepi)
subplot(numplots,1,3)
plot([tarousal,tarousal], [min(Pepi) max(Pepi)],'r')
end

subplot(numplots,1,4)
plot(Time(indVI(aind-4:aind)), vol(indVI(aind-4:aind)), 'g+');
plot(Time(indVE(aind-4:aind)), vol(indVE(aind-4:aind)), 'go');

plot([tarousal,tarousal], [min(vol) max(vol)],'r')

subplot(numplots,1,5)
plot([tarousal,tarousal], [min(Vdot) max(Vdot)],'r')

if numplots==7
    subplot(numplots,1,7);
    plot([tarousal,tarousal], [min(GG) max(GG)],'r')
end

%calculate ventilation 5 breaths before arousal
var = mean(VdotE(aind-4:aind));
Ti_ar = mean(Ti(aind-4:aind));
Te_ar = mean(Te(aind-4:aind));
Ttot_ar = mean(Ttot(aind-4:aind));
VT_ar = mean(TidalVol(aind-4:aind));
cpap = Pmask(indVE); 
AvgCPAP = mean(cpap(aind-4:aind));

ARbreaths = VdotE(aind-4:aind);
if max(ARbreaths)-min(ARbreaths) > 0.5*max(ARbreaths)
    var = NaN;
end

%Plot arousal breaths
subplot(numplots,1,6);
plot(Time(indVI(1:end-1)), VdotE,'k.-')
hold on
if ~isnan(var)
    plot(Time(indVI(aind-4:aind)), VdotE(aind-4:aind),'r*');
    plot([Time(indVI(aind-4)),Time(indVI(aind))], [var var],'r');
    text(double(Time(indVI(aind-4)))-1,mean(VdotE)+1, strcat('Var mean =', num2str(var,3)))
else
    text(double(Time(indVI(aind-4)))-1,mean(VdotE)+1, 'Var mean = nan')
end
ytemp = ylim(ax(6));
plot([tarousal,tarousal], ytemp,'r')
ylabel('VdotE (L/min)')
if ~isnan(var)
    %calculate pepi delta right before arousal
    % Find peaks/valleys of Pepi signal
    if ~isempty(Pepi)
        %low pass filter pepi signal, cutoff 2hz
        wn = 2/(1/(2*Fs));
        [den,num] = butter(2, wn);
        pfilt = filtfilt(den,num,Pepi);

        yp = diff(pfilt)/Fs;
        ypp = diff(yp)/Fs;
        ypp(end+1:end+2) = 0;
        yp(end+1) = 0;

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

        for i = 1:length(mx(:,1))
            [peakpepi(i), peakind(i)] = min(ypp(mx(i,1):mn(i,1)));
            peakind(i) = peakind(i)+mx(i,1)-1;
        end

        [b,bI] = min(abs(Time(mn(:,1))-Time(indVI(aind))));

        subplot(numplots,1,3)
        plot(Time(mn(:,1)), Pepi(mn(:,1)),'ro');
        plot(Time(peakind), Pepi(peakind),'r+');
        h2a = plot(Time(peakind(bI)), Pepi(peakind(bI)),'g+');
        h2b = plot(Time(mn(bI,1)), Pepi(mn(bI,1)),'go');

        KeepPepi = listdlg('PromptString','Repick Pepi max and min?',...
            'SelectionMode','single', 'ListSize', [150 50],...
            'ListString',{'Yes','No', 'Pepi is crap'});
        if KeepPepi==1
            subplot(numplots,1,3)
            [mxNew,mnNew] = SelectPepi(Time, Pepi, pfilt, Vdot, tarousal, stage, evts,1);
            set(h2a,'visible', 'off');
            set(h2b,'visible', 'off');
            plot(Time(mnNew), Pepi(mnNew),'go')
            plot(Time(mxNew), Pepi(mxNew),'g+')
            DelPepi = Pepi(mxNew)-Pepi(mnNew);
        elseif KeepPepi == 2;
            DelPepi = Pepi(peakind(bI))-Pepi(mn(bI,1));
        elseif KeepPepi==3
            DelPepi = nan;
        end
        text(double(Time(peakind(bI))), Pepi(peakind(bI))-2, num2str(round(DelPepi)));
    else
        DelPepi = nan;
    end

    %GG analysis
    if ~isempty(GG)
        indMn = find(TimeGG>=Time(indVI(aind-1)) & TimeGG<=Time(indVE(aind)));
        indMx = find(TimeGG>=Time(indVE(aind)) & TimeGG<=Time(indVI(aind)));
        [mxGG, mxIndGG] = max(GG(indMx));
        [mnGG, mnIndGG] = min(GG(indMn));
        mxIndGG = mxIndGG + indMx(1)-1;
        mnIndGG = mnIndGG + indMn(1)-1;

        GGtonic = mnGG;
        GGpeak =mxGG;
        GGphasic = mxGG-mnGG;

        subplot(numplots,1,7);
        plot(TimeGG(mxIndGG), mxGG,'r+');
        plot(TimeGG(mnIndGG), mnGG, 'ro');
        axis([TimeGG(1) TimeGG(end), min(GG), max(GG)]);
    else
        GGphasic = nan;
        GGtonic = nan;
        GGpeak = nan;
    end
else
    DelPepi = nan;
    GGphasic = nan;
    GGtonic = nan;
    GGpeak = nan;
end
end
