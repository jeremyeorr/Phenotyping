function [Breaths] = GGanalysis(Fs, Time, Vdot, Pmask, stage, evts, Pepi, TimeGG, GG)
% Integrate the raw flow signal to get volume
vol =cumtrapz(Vdot*Fs);
%detrend volume signal to remove leak
vol = detrend(vol, 'linear');
%Calculate end inspiration and end expiration points, tidal volumes
[TidalVol,indVI,indVE] = peakVol(Time,vol,0.1);
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

%Have User select where to analyze
subplot(numplots,1,1)
title('Select start and end where flow limitation is occuring')
[x,y] = ginput(2); 
indxNew = find(Time>=x(1) & Time<=x(2));

% Integrate the raw flow signal to get volume
volNew =cumtrapz(Vdot(indxNew)*Fs);
%detrend volume signal to remove leak
volNew = detrend(volNew, 'linear');
TimeNew = Time(indxNew);
[TidalVol,indVI,indVE] = peakVol(Time(indxNew),volNew,0.1);

%Rescale figure
axis([x(1) x(2) -1 6])
subplot(numplots,1,2)
axis([x(1) x(2) min(Pmask(indxNew)) max(Pmask(indxNew))])
if ~isempty(Pepi)
    subplot(numplots,1,3)
    axis([x(1) x(2) min(Pepi(indxNew)) max(Pepi(indxNew))])
end
subplot(numplots,1,4)
set(v1temp, 'visible', 'off')
set(hI, 'Visible', 'off')
set(hE, 'Visible', 'off')
plot(TimeNew, volNew,'k'); 
hI = plot(TimeNew(indVI), volNew(indVI),'r+');
hE = plot(TimeNew(indVE), volNew(indVE),'ro');

axis([x(1) x(2) min(volNew) max(volNew)])
subplot(numplots,1,5)
axis([x(1) x(2) min(Vdot(indxNew)) max(Vdot(indxNew))])

% Fix any volume points?
RedoVol = listdlg('PromptString','Fix any volumes?',...
    'SelectionMode','single', 'ListSize', [150 50],...
    'ListString',{'Yes','No'});

if RedoVol==1
    [indVI, indVE] = FixVolumes(TimeNew, volNew, indVI, indVE);
    TidalVol = volNew(indVI)-volNew(indVE);
    set(hI, 'Visible', 'off')
    set(hE, 'Visible', 'off')
    subplot(numplots,1,4)
    plot(TimeNew(indVI), volNew(indVI),'r+');
    plot(TimeNew(indVE), volNew(indVE),'ro');
end
cpap = Pmask(indxNew(indVE));
cpap(end) = [];
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
AvgVdotE = mean(VdotE);


%Plot VdotE
figure(fh)
ax(6) = subplot(numplots,1,6);
plot(TimeNew(indVI(1:end-1)), VdotE,'k.-')
hold on
plot([TimeNew(indVI(1)), TimeNew(indVI(end))], [AvgVdotE,AvgVdotE],'r')
text(double(TimeNew(indVI(round(length(indVI)/2))))-1,round(AvgVdotE)+1, strcat('VdotE mean =', num2str(AvgVdotE,3)))
ylabel('VdotE (L/min)')
xlim([x(1) x(2)]);
% Find peaks/valleys of Pepi signal
if ~isempty(Pepi)
    Pepi = Pepi(indxNew);
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
    subplot(numplots,1,3)
    plot(TimeNew(mn(:,1)), Pepi(mn(:,1)),'ro')
    plot(TimeNew(peakind), Pepi(peakind),'r+')
    DelPepi = Pepi(peakind)-Pepi(mn(:,1));
    PepiTime = Time(peakind);
    
    KeepPepi = listdlg('PromptString','Keep Pepi?',...
                       'SelectionMode','single', 'ListSize', [150 50],...
                       'ListString',{'Yes','Pepi is crap'});
    if KeepPepi==2
        DelPepi = nan;
    else
       for i = 1:length(indVI)
           a = find(abs(TimeNew(indVI(i))-TimeNew(mn(:,1)))<1.5);
           if isempty(a)
               PeakPepi(i,:) = [0,0];
               MinPepi(i,:) = [0,0];
           elseif length(a)==1
               PeakPepi(i,:) = [TimeNew(peakind(a)), Pepi(peakind(a))];
               MinPepi(i,:) = [TimeNew(mn(a,1)), Pepi(mn(a,1))];
           else
              [b,bI] = min(TimeNew(indVI(i))-TimeNew(mn(:,1)));
              PeakPepi(i,:) = [TimeNew(peakind(bI)), Pepi(peakind(bI))];
              MinPepi(i,:) = [TimeNew(mn(bI,1)), Pepi(mn(bI,1))];
           end
       end
    end
    
else
    DelPepi = nan;
end

%GG analysis
if ~isempty(GG)
      for i = 2:length(indVI)
        indMn = find(TimeGG>=TimeNew(indVI(i-1)) & TimeGG<=TimeNew(indVE(i)));
        indMx = find(TimeGG>=TimeNew(indVE(i)) & TimeGG<=TimeNew(indVI(i)));
      [mnGG(i-1), mnIndGG(i-1)] = min(GG(indMn));  
      [mxGG(i-1), mxIndGG(i-1)] = max(GG(indMx));
       mxIndGG(i-1) = mxIndGG(i-1) + indMx(1)-1;
       mnIndGG(i-1) = mnIndGG(i-1) + indMn(1)-1;
    end
      
    GGtonic = mnGG; 
    GGpeak = mxGG;
    GGphasic = mxGG-mnGG; 
    
    subplot(numplots,1,7);
    plot(TimeGG(mxIndGG), mxGG,'r+'); 
    plot(TimeGG(mnIndGG), mnGG, 'ro');
    plot([TimeGG(mnIndGG(1)), TimeGG(mnIndGG(end))], [mean(GGtonic), mean(GGtonic)],'g')
    plot([TimeGG(mxIndGG(1)), TimeGG(mxIndGG(end))], [mean(GGphasic), mean(GGphasic)],'b')
    plot([TimeGG(mxIndGG(1)), TimeGG(mxIndGG(end))], [mean(GGpeak), mean(GGpeak)],'r')
else
    GGphasic = nan;
    GGtonic = nan;
    GGpeak = nan;
end
length(PepiTime)
length(GGtonic)
length(indVI)

Breaths.Time(:,1) = TimeNew(indVE(2:end-1));
Breaths.Time(:,2) = TimeNew(indVI(2:end-1));
Breaths.VdotE(:,1) = VdotE(2:end); 
Breaths.Ti(:,1) = Ti(2:end); 
Breaths.Te(:,1) = Te(2:end); 
Breaths.Ttot(:,1) = Ttot(2:end); 
Breaths.Vt(:,1) = TidalVol(2:end); 
Breaths.GGtonic(:,1) = GGtonic(1:end-1);
Breaths.GGphasic(:,1) = GGphasic(1:end-1);
Breaths.GGpeak(:,1) = GGpeak(1:end-1);
Breaths.Cpap(:,1) = cpap(2:end);
Breaths.Pepi(:,1) = PeakPepi(2:end-1,2)-MinPepi(2:end-1,2);

a = find(Breaths.Pepi(:,1)==0);
Breaths.Pepi(a) = [];
Breaths.VdotE(a) = [];
Breaths.Ti(a) = []; 
Breaths.Te(a) = []; 
Breaths.Ttot(a) = []; 
Breaths.Vt(a) = []; 
Breaths.GGtonic(a) = [];
Breaths.GGphasic(a) = [];
Breaths.GGpeak(a) = [];
Breaths.Cpap(a) = [];
Breaths.Time(a,:) = [];