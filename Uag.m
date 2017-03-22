% code for getting active v0 and the Upper Airway Gain (uag)
function [Vactive, Vdrive, Ug, Cdown, DelPepi, GGphasic, GGtonic, GGpeak, Ti_Vactive, Te_Vactive, Ttot_Vactive, VT_Vactive, Vdown, DelPepipredrop, CPAPdrop] = Uag(Fs, Time, Vdot, Pmask, stage, evts, Pepi, TimeGG, GG, Veupnea, LG, VV0, O2)
% code for getting active v0 and the Upper Airway Gain (uag)
%Lisa M Campana, PhD
%December 2012

%Calc volume by integrating flow signal
vol = cumtrapz(Vdot)*Fs;

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
end
hold on
ylabel('Pepi (cmH2O)')

ax(4) = subplot(numplots,1,4); %Volume
v1temp = plot(Time, vol, 'k');
hold on;
ylabel('Volume (L)')

ax(5) = subplot(numplots,1,5); %Flow
plot(Time, Vdot,'k')
hold on
plot([Time(1), Time(end)], [0,0],'r')
ylabel('Flow (L/s)')

ax(6) = subplot(numplots,1,6);
if numplots==7
    ax(7) = subplot(numplots,1,7);
    plot(TimeGG, GG,'k')
    ylabel('GG')
    hold on
end
linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
% set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])

%Select where pressure is dropped
subplot(numplots,1,1)
title('Select Begining of drop');
[x,y] = ginput(1);
tdown = x(1);
clear x y

subplot(numplots,1,1)
title('Select End of drop');
[x,y] = ginput(1);
tup = x(1);
clear x y;

% Determine Leak
Region1 = find(Time<=tdown);
Region2 = find(Time>tdown & Time<=tup);
Region3 = find(Time>tup);

p1 = polyfit(Time(Region1),vol(Region1),1);
p2 = polyfit(Time(Region2),vol(Region2),1);
p3 = polyfit(Time(Region3),vol(Region3),1);

Vdot(Region1) = Vdot(Region1) - p1(1);
Vdot(Region2) = Vdot(Region2) - p2(1);
Vdot(Region3) = Vdot(Region3) - p3(1);

vol = cumtrapz(Vdot)*Fs; % Andrew's volume
vol2 = detrend(vol, 'linear', [Region2(1), Region3(1)]);
[TidalVol,indVI,indVE] = peakVol(Time,vol,0.05);

subplot(numplots,1,4)
plot(Time, vol, 'r')
plot(Time, vol2,'b')
set(v1temp,'visible', 'off');
hI = plot(Time(indVI), vol(indVI),'r+');
hE = plot(Time(indVE), vol(indVE),'ro');
hI2 = plot(Time(indVI), vol2(indVI),'c+');
hE2 = plot(Time(indVE), vol2(indVE),'co');
axis([Time(Region2(1)) Time(Region2(end)) min(vol(Region2)) max(vol(Region2))])

subplot(numplots,1,5)
axis([Time(Region2(1)) Time(Region2(end)) min(Vdot(Region2)) max(Vdot(Region2))])

%Is there 0 flow?
subplot(numplots,1,1)
title('Check to see if there is zero flow');

[zf] = listdlg('PromptString','Is there zero flow?',...
    'SelectionMode','single', 'ListSize', [150 50],...
    'ListString',{'Yes','No'});
if zf == 1
    Vactive = 0;
    Ti_Vactive = nan;
    Te_Vactive = nan;
    Ttot_Vactive = nan;
    VT_Vactive = 0;
    Ug = 0;
    Vdown = 0;
    Vdrive = 0;
    CPAPdrop=mean(Pmask(Region2));
    Cdown = mean(Pmask(Region1));
    timeVI = Time(indVI)
    breathindx = find(timeVI<tdown,1,'last');
    numbreath = 2;
    subplot(numplots,1,4)
    axis auto
else
    %Ask if we need to redo volume calculation
    CPAPdrop=mean(Pmask(Region2));
    RedoVol = listdlg('PromptString','Fix any volumes?',...
        'SelectionMode','single', 'ListSize', [150 50],...
        'ListString',{'Yes','No'});

    if RedoVol==1
        [indVI, indVE] = FixVolumes(Time, vol, indVI, indVE);
        TidalVol = vol(indVI)-vol(indVE);
        set(hI, 'Visible', 'off')
        set(hE, 'Visible', 'off')
        set(hI2, 'Visible', 'off')
        set(hE2, 'Visible', 'off')
        subplot(numplots,1,4)
        plot(Time(indVI), vol(indVI),'r+');
        plot(Time(indVE), vol(indVE),'ro');
        plot(Time(indVI), vol2(indVI),'c+');
        plot(Time(indVE), vol2(indVE),'co');
    end
    timeVI = Time(indVI);
    timeVE = Time(indVE);

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
    
    %Lisa's tidal volume
    TidalVol2 = vol2(indVI)-vol2(indVE);
    VdotE2 = [];
    for n = 1:length(indVI)-1
        VdotE2(n) = (TidalVol2(n))/(Ttot(n));
    end
    VdotE2 = VdotE2.*60;
    a = find(VdotE2<0);
    VdotE2(a) = 0;
    VdotE2 = VdotE2';
    clear a n;

%     aind = find(timeVI<tdown,1,'last');
%     aind = aind+1; %index is of first breath of the drop
    cpap = Pmask(indVE);
    cpap(end) = [];
    
    subplot(numplots,1,1);
    title('Select breaths: click to the left of first breath and to the right of the last breath you wish to analyze');
    [x2,y2] = ginput(2);
    breaths= find(Time(indVI)>=x2(1,1) & Time(indVI)<=x2(2,1));
    aind=breaths(1);
    Vactive = mean(VdotE(breaths));
    Ti_Vactive = mean(Ti(breaths));
    Te_Vactive = mean(Te(breaths));
    Ttot_Vactive = mean(Ttot(breaths));
    VT_Vactive = mean(TidalVol(breaths));
    
%     if  Time(indVI(aind+2))<tup
%         Vactive = mean(VdotE(aind+1:aind+2)); %breaths 2 and 3 of drop
%         Ti_Vactive = mean(Ti(aind+1:aind+2));
%         Te_Vactive = mean(Te(aind+1:aind+2));
%         Ttot_Vactive = mean(Ttot(aind+1:aind+2));
%         VT_Vactive = mean(TidalVol(aind+1:aind+2));
%     else
%         Vactive = VdotE(aind+1);
%         Ti_Vactive = (Ti(aind+1));
%         Te_Vactive = (Te(aind+1));
%         Ttot_Vactive = (Ttot(aind+1));
%         VT_Vactive = (TidalVol(aind+1));
%     end
    
    breathindx = find(timeVI<tdown,1,'last');
    %breathindx = aind+1; %index is of first breath of the drop
    if breathindx>5
        Vdown = mean(VdotE(breathindx-4:breathindx)); %5 breaths before drop
        Cdown = mean(cpap(breathindx-4:breathindx));
    else
        Vdown = mean(VdotE(1:breathindx));
        Cdown = mean(cpap(1:breathindx));
    end
    if ~isnan(LG)
        Vdrive = (Veupnea-Vdown)*LG;
    else
        Vdrive = NaN;
    end
    Ug = (Vactive-VV0)/Vdrive;


    %Plot VdotE, Vup, Vdown, Veupnea
    subplot(numplots,1,4)
    subplot(numplots,1,4);
    plot(Time(indVI(breaths)), vol(indVI(breaths)),'g+');
    plot(Time(indVE(breaths)), vol(indVE(breaths)),'go');
    plot(Time(indVE(breaths(1)):indVI(breaths(end))), vol(indVE(breaths(1)):indVI(breaths(end))), 'g')
    if breathindx>5
        plot(Time(indVI(breathindx-4:breathindx)), vol(indVI(breathindx-4:breathindx)),'g+');
        plot(Time(indVE(breathindx-4:breathindx)), vol(indVE(breathindx-4:breathindx)),'go');
        plot(Time(indVE(breathindx-4):indVI(breathindx)), vol(indVE(breathindx-4):indVI(breathindx)), 'g')
    else
        plot(Time(indVI(1:breathindx)), vol(indVI(1:breathindx)),'g+');
        plot(Time(indVE(1:breathindx)), vol(indVE(1:breathindx)),'go');
        plot(Time(indVE(1):indVI(breathindx)), vol(indVE(1):indVI(breathindx)), 'g')
    end
    axis([Time(1) Time(end) min(vol) max(vol)])

    subplot(numplots,1,6);
    plot(Time(indVI(1:end-1)), VdotE,'k.-')
    hold on
    ylabel('VdotE (L/min)')
    plot([Time(indVI(1)),Time(indVI(end))], [Veupnea Veupnea],'k--');
    if breathindx>5
        plot([Time(indVI(breathindx-4)),Time(indVI(breathindx))], [Vdown Vdown],'b');
    else
        plot([Time(indVI(1)),Time(indVI(breathindx))], [Vdown Vdown],'b');
    end
    plot([Time(indVI(breaths(1))),Time(indVI(breaths(end)))], [Vactive Vactive],'g');
    axis([Time(1) Time(end) min(VdotE) max(VdotE)])
    text(double(Time(indVI(breaths(1))))-1,mean(VdotE)+1, strcat('Vactive mean =', num2str(Vactive,3)))
end
  
%calculate pepi delta right before arousal
% Find peaks/valleys of Pepi signal
if ~isempty(Pepi)
    %low pass filter pepi signal, cutoff 2hz
    wn = 2/(1/(2*Fs));
    [den,num] = butter(2, wn);
    pfilt = filtfilt(den,num,Pepi);

    cutoff=5;
    wn = cutoff/(1/(2*Fs)); 
    [den,num] = butter(2, wn);
    pfilt2 = filtfilt(den,num,Pepi);
    
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
    plot(Time(mn(:,1)), Pepi(mn(:,1)),'ro')
    plot(Time(peakind), Pepi(peakind),'r+')
    
    %correct Pepi tracing
    KeepPepi = listdlg('PromptString','Pepi quality',...
                       'SelectionMode','single', 'ListSize', [150 50],...
                       'ListString',{'Keep','Do not use','Repick'});
    if KeepPepi==1
        allPepi=[];
        allPepi=pfilt2(peakind)-pfilt2(mn(:,1));
    end
    if KeepPepi==2
        DelPepi = nan;
        allPepi=[];
    end
    if KeepPepi==3
        [mn, peakind] = FixPepi(Time, pfilt2, Vdot, mn, peakind);
        peakind'
        allPepi=[];
        allPepi=pfilt2(peakind)-pfilt2(mn(:,1));
        figure(fh)
        subplot(numplots,1,3)
        plot(Time(mn(:,1)), pfilt2(mn(:,1)),'g+');
        plot(Time(peakind), pfilt2(peakind),'go' );
    end
    
    %Calculating Pepi swing on Vactive breaths
    if isempty(allPepi)
       DelPepipredrop=[];
    else
       DelPepipredrop=mean(allPepi(breathindx-4):allPepi(breathindx)); 
    end
    
    %Calculating Pepi swing on Vactive breaths
    if zf==2 %not zero flow
        for i = 1:length(breaths) 
            [b,bI(i)] = min(abs(Time(mn(:,1))-Time(indVI(breaths(1)-1+i))));
        end
        DelPepi = mean(Pepi(peakind(bI))-Pepi(mn(bI,1)));
        h2a = plot(Time(mn(bI,1)), Pepi(mn(bI,1)),'go');
        h2b = plot(Time(peakind(bI)), Pepi(peakind(bI)),'g+');
    elseif zf==1 %zero flow
        aind = find(Time(mn(:,1))<tdown,1,'last');
        aind = aind+1;
        DelPepi = mean(Pepi(peakind(aind+2:aind+3))-Pepi(mn(aind+2:aind+3,1)));
        h2a = plot(Time(mn(aind+2:aind+3,1)), Pepi(mn(aind+2:aind+3,1)),'go');
        h2b = plot(Time(peakind(aind+2:aind+3)), Pepi(peakind(aind+2:aind+3)),'g+');
    end
    
  

%GG analysis
if ~isempty(GG)

    if zf==2
        for i = 1:length(breaths) 
            indMn = find(TimeGG>=Time(indVI(breaths(1)-2+i)) & TimeGG<=Time(indVE(breaths(1)-2+i+1)));
            indMx = find(TimeGG>=Time(indVE(breaths(1)-2+i+1)) & TimeGG<=Time(indVI(breaths(1)-2+i+1)));
            [mxGG(i), mxIndGG(i)] = max(GG(indMx));
            [mnGG(i), mnIndGG(i)] = min(GG(indMn));
            mxIndGG(i) = mxIndGG(i) + indMx(1)-1;
            mnIndGG(i) = mnIndGG(i) + indMn(1)-1;
        end
    elseif zf==1
        for i = 1:2
            indMn = find(TimeGG>=Time(mn(aind+i-1,1)) & TimeGG<=Time(peakind(aind+i)));
            indMx = find(TimeGG>=Time(peakind(aind+i)) & TimeGG<=Time(mn(aind+i,1)));
            [mxGG(i), mxIndGG(i)] = max(GG(indMx));
            [mnGG(i), mnIndGG(i)] = min(GG(indMn));
            mxIndGG(i) = mxIndGG(i) + indMx(1)-1;
            mnIndGG(i) = mnIndGG(i) + indMn(1)-1;
        end
    end
    GGtonic = mean(mnGG);
    GGpeak =mean(mxGG);
    GGphasic = mean(mxGG-mnGG);

    subplot(numplots,1,7);
    plot(TimeGG(mxIndGG), mxGG,'g+');
    plot(TimeGG(mnIndGG), mnGG, 'go');
    axis([TimeGG(1) TimeGG(end), min(GG), max(GG)]);
    
else
    GGphasic = nan;
    GGtonic = nan;
    GGpeak = nan;
end

end

