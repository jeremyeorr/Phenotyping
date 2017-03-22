
function [v0, Ti_v0, Te_v0, Ttot_v0, VT_v0, CPAPdrop] = Vzero(Fs, Time, Vdot, Pmask, stage, evts)
%Function to calculate passive ventilation
%Lisa M Campana, PhD
%December 2012
thresh=0.01;
%Calc volume by integrating flow signal
vol = cumtrapz(Vdot)*Fs; %note that Fs is not really sampling rate but 1/sampling rate (1/Hz).



% Plot figure
numplots = 5;
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

ax(3) = subplot(numplots,1,3); %Volume
v1temp = plot(Time, vol, 'k');
hold on;
ylabel('Volume (L)')

ax(4) = subplot(numplots,1,4); %Flow
plot(Time, Vdot,'k')
hold on
plot([Time(1), Time(end)], [0,0],'r')
ylabel('Flow (L/s)')

ax(5) = subplot(numplots,1,5);
linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
% set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])

%Select Beginning and end of CPAP drop
subplot(numplots,1,1)
title('Select Beginning of drop');
[x,y] = ginput(1);
tdown = x(1);
clear x y

subplot(numplots,1,1)
title('Select End of drop');
[x,y] = ginput(1);
tup = x(1);
clear x y;

% Detrend volume signal according to each pressure setting
Region1 = find(Time<=tdown);
Region2 = find(Time>tdown & Time<=tup);
Region3 = find(Time>tup);

leakRegion2=mean(Vdot(Region2));
Vdot_predetrend = Vdot;
Vdot=Vdot-leakRegion2;
vol_predetrend = vol;
vol = cumtrapz(Vdot)*Fs;
%Calculated actual CPAP level during the drop
CPAPdrop=mean(Pmask(Region2));

[TidalVol,indVI,indVE] = peakVol(Time,vol,thresh); %lisa had 0.1

figure(200)
plot(Time,vol_predetrend,Time,vol,Time,vol_predetrend-vol)
figure(201)
plot(Time,Vdot_predetrend,Time,Vdot,Time,Vdot_predetrend-Vdot)
figure(fh) %come back to main figure.

%Plot new volumes
subplot(numplots,1,3)
plot(Time, vol, 'r')
hI = plot(Time(indVI), vol(indVI),'r+');
hE = plot(Time(indVE), vol(indVE),'ro');
set(v1temp, 'visible', 'off');
axis([tdown-10 tup+10 min(vol(Region2))-0.1 max(vol(Region2))+0.1])

%Is there zero flow?
subplot(numplots,1,1)
title('Check to see if there is zero flow');
subplot(numplots,1,4)
axis([tdown-10 tup+10 min(Vdot(Region2))-0.1 max(Vdot(Region2))+0.1])
[zf] = listdlg('PromptString','Is there zero flow?',...
    'SelectionMode','single', 'ListSize', [150 50],...
    'ListString',{'Yes','No'});
if zf == 1 %if zero flow then V0 = 0
    v0 = 0;
    Ti_v0 = nan;
    Te_v0 = nan;
    Ttot_v0 = nan;
    VT_v0 = 0;
    aind = [];
    subplot(numplots,1,5)
    plot([tdown tup], [0 0],'k')
    axis([Time(1) Time(end) -1 1])
     subplot(numplots,1,3)
    axis([Time(1) Time(end) min(vol) max(vol)])
     subplot(numplots,1,4)
   axis([Time(1) Time(end) min(Vdot) max(Vdot)])
else
    %Fix any volume points?
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

    
    %aind = find(timeVI<tdown,1,'last');
    %aind = aind+1;
    
%     %Determine breaths to analyze (breath 3 and 4 after CPAP drop) old
%     code by Lisa
%     v0 = mean(VdotE(aind+2:aind+3));
%     Ti_v0 = mean(Ti(aind+2:aind+3));
%     Te_v0 = mean(Te(aind+2:aind+3));
%     Ttot_v0 = mean(Ttot(aind+2:aind+3));
%     VT_v0 = mean(TidalVol(aind+2:aind+3));

    % Choose which breaths to analyze, calculate and plot breaths
    subplot(numplots,1,1);
    title('Select breaths: click to the left of first breath and to the right of the last breath you wish to analyze');
    [x2,y2] = ginput(2);
    breaths= find(Time(indVI)>=x2(1,1) & Time(indVI)<=x2(2,1));
    aind=breaths(1);
    v0 = mean(VdotE(breaths));
    Ti_v0 = mean(Ti(breaths));
    Te_v0 = mean(Te(breaths));
    Ttot_v0 = mean(Ttot(breaths));
    VT_v0 = mean(TidalVol(breaths));
    
    subplot(numplots,1,3);
    plot(Time(indVI(breaths)), vol(indVI(breaths)),'g+');
    plot(Time(indVE(breaths)), vol(indVE(breaths)),'go');
    axis([Time(1) Time(end) min(vol) max(vol)])

    subplot(numplots,1,5);
    plot(Time(indVI(1:end-1)), VdotE,'k.-')
    hold on
    plot(Time(indVI(breaths)), VdotE(breaths),'g*')
    axis([Time(1) Time(end) min(VdotE) max(VdotE)])
    text(double(Time(indVI(aind)))-1,mean(VdotE)+1, strcat('V0 mean =', num2str(v0,3)))
    ylabel('VdotE (L/min)')

    subplot(numplots,1,5)
    axis auto

    hold on
    if ~isnan(v0)
        subplot(numplots,1,4)
        plot([Time(indVI(aind+1)),Time(indVI(aind+1))], [v0 v0],'r');
    end
end



end
