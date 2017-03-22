function [Breaths] = GGanalysis(Fs, Time, Vdot, Pmask, stage, evts, Pepi, TimeGG, GG)
% Integrate the raw flow signal to get volume
vol =cumtrapz(Vdot*Fs);
%detrend volume signal to remove leak
vol = detrend(vol, 'linear');
%Calculate end inspiration and end expiration points, tidal volumes
[indVI,indVE] = peakVol2(Time,vol,0.1);
if length(indVE)>length(indVI)
    extra_indVE = indVE(end);
    indVE(end)=[];
else
    extra_indVE = NaN;
end
TidalVol=vol(indVI)-vol(indVE);

cpap = Pmask(indVE);
% cpap(end) = [];
numplots = 7;

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

ax(7) = subplot(numplots,1,7);
plot(TimeGG, GG,'k');
ylabel('GG')
hold on;

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

% Determine vdote
clear Ttot Ti Te VdotE
for n = 1:length(indVE);
    Ti(n) = Time(indVI(n))-Time(indVE(n));
    if n<length(indVE)
        Ttot(n) = Time(indVE(n+1))-Time(indVE(n));
        Te(n) = Ttot(n)-Ti(n);
        VdotE(n) = (TidalVol(n))/(Ttot(n));
    else
        Ttot(n) = NaN;
        Te(n) = NaN;
        VdotE(n) = NaN;
    end
end
VdotE = VdotE.*60;
a = find(VdotE<0);
VdotE(a) = 0;
VdotE = VdotE';
clear n;
AvgVdotE = nanmean(VdotE);


%Plot VdotE
figure(fh)
ax(6) = subplot(numplots,1,6);
plot(Time(indVI(1:end-1)), VdotE(1:end-1),'k.-')
hold on
plot([Time(indVI(1)), Time(indVI(end))], [AvgVdotE,AvgVdotE],'r')
% text(double(Time(indVI(round(length(indVI)/2))))-1,round(AvgVdotE)+1, strcat('VdotE mean =', num2str(AvgVdotE,3)))
ylabel('VdotE (L/min)')
xlim([Time(1) Time(end)]);

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
    plot(Time(mn(:,1)), Pepi(mn(:,1)),'ro')
    plot(Time(peakind), Pepi(peakind),'r+')
    
    
    KeepPepi = listdlg('PromptString','Pepi quality',...
        'SelectionMode','single', 'ListSize', [150 50],...
        'ListString',{'Keep','Do not use','Repick'});
    if KeepPepi==1
        DelPepi=[];
        DelPepi=pfilt2(peakind)-pfilt2(mn(:,1));
    end
    if KeepPepi==2
        DelPepi = nan;
    end
    if KeepPepi==3
        [mn, peakind] = FixPepi(Time, pfilt2, Vdot, mn, peakind);
        peakind'
        DelPepi = pfilt2(peakind)-pfilt2(mn(:,1));
        figure(fh)
        subplot(numplots,1,3)
        plot(Time(mn(:,1)), pfilt2(mn(:,1)),'g+');
        plot(Time(peakind), pfilt2(peakind),'go' );
    end
    PepiTime = Time(peakind);
    PepiTimeNadir = Time(mn(:,1));
    PepiTimeNadiri = mn(:,1);
else
    DelPepi = nan;
end

%Sync Pepi with breaths from Ventilation analysis
clear DelPepi2 DelPepi2_time
for i=1:length(indVE)
    temp=(PepiTimeNadir-Time(indVE(i)));
    tempi = find(temp>0,1);
    if ~isempty(tempi)
        DelPepi2(i) = DelPepi(tempi);
        DelPepi2_time(i) = PepiTimeNadir(tempi);
        DelPepi2_i(i) = PepiTimeNadiri(tempi);
    else
        DelPepi2(i) = NaN;
        DelPepi2_time(i) = NaN;
        DelPepi2_i(i) = NaN;
    end
end

%GG analysis
clear mnGG mxGG mnIndGG mxIndGG
if ~isempty(GG)
    mnGG(1)=NaN;
    mnIndGG(1)=NaN;
    for i = 1:length(indVI)
        indMx = find(TimeGG>=Time(indVE(i)) & TimeGG<=Time(indVI(i)));
        [mxGG(i), mxIndGG(i)] = max(GG(indMx));
        mxIndGG(i) = mxIndGG(i) + indMx(1)-1;
    end
    for i = 2:length(indVI)
        indMn = find(TimeGG>=Time(indVI(i-1)) & TimeGG<=Time(indVE(i)));
        [mnGG(i), mnIndGG(i)] = min(GG(indMn));
        mnIndGG(i) = mnIndGG(i) + indMn(1)-1;
    end
    
    GGtonic = mnGG;
    GGpeak = mxGG;
    GGphasic = mxGG-mnGG;
    
    subplot(numplots,1,7);
    plot(TimeGG(mxIndGG), mxGG,'r+');
    plot(TimeGG(mnIndGG(~isnan(mnIndGG))), mnGG(~isnan(mnIndGG)), 'ro');
    %     plot([TimeGG(mnIndGG(1)), TimeGG(mnIndGG(end))], [mean(GGtonic), mean(GGtonic)],'g')
    %     plot([TimeGG(mxIndGG(1)), TimeGG(mxIndGG(end))], [mean(GGphasic), mean(GGphasic)],'b')
    %     plot([TimeGG(mxIndGG(1)), TimeGG(mxIndGG(end))], [mean(GGpeak), mean(GGpeak)],'r')
else
    GGphasic = nan;
    GGtonic = nan;
    GGpeak = nan;
end

clear Breaths
Breaths.Time_i(:,1) = indVE(1:end);
Breaths.Time_i(:,2) = indVI(1:end);
Breaths.Time(:,1) = Time(indVE(1:end));
Breaths.Time(:,2) = Time(indVI(1:end));
Breaths.VdotE(:,1) = [VdotE(1:end)];
Breaths.Ti(:,1) = Ti(1:end);
Breaths.Te(:,1) = Te(1:end);
Breaths.Ttot(:,1) = Ttot(1:end);
Breaths.Vt(:,1) = TidalVol(1:end);
Breaths.GGtonic(:,1) = GGtonic(1:end);
Breaths.GGphasic(:,1) = GGphasic(1:end);
Breaths.GGpeak(:,1) = GGpeak(1:end);
Breaths.GGpeak_i(:,1) = mxIndGG(1:end);
Breaths.GGtonic_i(:,1) = mnIndGG(1:end);
Breaths.Cpap(:,1) = cpap(1:end);
Breaths.Pepi(:,1) = -DelPepi2(1:end);
Breaths.PepiTime(:,1) = DelPepi2_time(1:end);
Breaths.Pepi_i(:,1) = DelPepi2_i(1:end);

%Have User select data to exclude
Excludedata = listdlg('PromptString','Exclude any data?',...
    'SelectionMode','single', 'ListSize', [150 50],...
    'ListString',{'No','Yes'});

if Excludedata==2
    title('Select data you would like to exclude')
    while 1
        % plot current volumes with start insp and peak insp
        % First you have the option to delete points
        
        figure(1)
        I=Breaths.Pepi_i(~isnan(Breaths.Pepi_i));
        ax1(1)=subplot(6,1,1); plot(Time,Pepi,'k',Time(I),Pepi(I),'r.');
        I=Breaths.Time_i(~isnan(Breaths.Time_i));
        ax1(2)=subplot(6,1,2); plot(Breaths.PepiTime,Breaths.Pepi,'g.','MarkerSize',22);
        ax1(3)=subplot(6,1,3); plot(Time,vol,'k',Time(I),vol(I),'r.');
        ax1(4)=subplot(6,1,4); plot(Time,Vdot,'k',Time(I),Vdot(I),'r.');
        ax1(5)=subplot(6,1,5); plot(Breaths.Time(:,2),Breaths.VdotE,'r.','MarkerSize',22);
        I=Breaths.GGtonic_i(~isnan(Breaths.GGtonic_i));
        I2=Breaths.GGpeak_i(~isnan(Breaths.GGpeak_i));
        ax1(6)=subplot(6,1,6); plot(TimeGG,GG,'k',TimeGG(I),Breaths.GGtonic(~isnan(Breaths.GGtonic_i)),'b.',TimeGG(I2),Breaths.GGpeak(~isnan(Breaths.GGpeak_i)),'r.');
        linkaxes(ax1,'x');
        [x,y,button]=ginput(1);
        if button==3
            break
            close(fh)
        end
        TimeFromClickToData = Breaths.Time(:,2)-x;
        [~,i]=min(abs(TimeFromClickToData));
        %Breaths.Time_i(i,1) = NaN;
        %Breaths.Time_i(i,2) = NaN;
        Breaths.Time(i,1) = NaN;
        Breaths.Time(i,2) = NaN;
        Breaths.VdotE(i,1) = NaN;
        Breaths.Ti(i,1) = NaN;
        Breaths.Te(i,1) = NaN;
        Breaths.Ttot(i,1) = NaN;
        Breaths.Vt(i,1) = NaN;
        Breaths.GGtonic(i,1) = NaN;
        Breaths.GGphasic(i,1) = NaN; % GGphasic(1:end);
        Breaths.GGpeak(i,1) = NaN; % GGpeak(1:end);
        %Breaths.GGpeak_i(i,1) = NaN; % GGpeak(1:end);
        %Breaths.GGtonic_i(i,1) = NaN; % GGpeak(1:end);
        Breaths.Cpap(i,1) = NaN; % cpap(1:end);
        Breaths.Pepi(i,1) = NaN; % -DelPepi2(1:end);
        Breaths.PepiTime(i,1) = NaN; % DelPepi2_time(1:end);
        %Breaths.Pepi_i(i,1) = NaN;
    end
end
