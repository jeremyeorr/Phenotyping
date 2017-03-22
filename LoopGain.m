function [LG, LG2, LG3, LG4, Vup, Vup_1br, Vdown, Vup2, Vup2_1br, Vdown2, Cdown] = LoopGain(Fs, Time, Vdot, Pmask, stage, evts, Pepi, Veupnea)
% Code for getting loop gain
% Must run the Veup script first
%Lisa M Campana, PhD
%December 2012

%Calc volume by integrating flow signal
vol = cumtrapz(Vdot)*Fs;

% Plot figure
fh = figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(6,1,1); %Sleep stage
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

ax(2) = subplot(6,1,2); % Pmask
plot(Time, Pmask,'k')
hold on;
ylabel('Pmask (cmH2O')

ax(3) = subplot(6,1,3);%Pepi (if it exists)
if ~isempty(Pepi)
    plot(Time, Pepi,'k')
end
hold on
ylabel('Pepi (cmH2O)')

ax(4) = subplot(6,1,4); %Volume
v1temp = plot(Time, vol, 'k');
hold on;
ylabel('Volume (L)')

ax(5) = subplot(6,1,5); %Flow
plot(Time, Vdot,'k')
hold on
plot([Time(1), Time(end)], [0,0],'r')
ylabel('Flow (L/s)')

ax(6) = subplot(6,1,6);

linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
% set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])

subplot(6,1,1)
title('Select Start of increase');
[x,y] = ginput(1);
tup = x(1);
clear x y

subplot(6,1,1)
title('Select End of increase');
[x,y] = ginput(1);
tdown = x(1);
clear x y;

% Determine Leak
Region1 = find(Time<=tup);
Region2 = find(Time>tup & Time<=tdown);
Region3 = find(Time>tdown);

p1 = polyfit(Time(Region1),vol(Region1),1);
p2 = polyfit(Time(Region2),vol(Region2),1);
p3 = polyfit(Time(Region3),vol(Region3),1);

Vdot(Region1) = Vdot(Region1) - p1(1);
Vdot(Region2) = Vdot(Region2) - p2(1);
Vdot(Region3) = Vdot(Region3) - p3(1);

vol2 = cumtrapz(Vdot)*Fs; %andrew's volume
vol = detrend(vol, 'linear', [Region2(1), Region3(1)]); %Lisa's volume

[TidalVol,indVI,indVE] = peakVol(Time,vol,0.1);
% [TidalVol2,indVI2,indVE2] = peakVol(Time,vol2,0.1);
subplot(6,1,4)
plot(Time, vol, 'r')
plot(Time, vol2,'b')
set(v1temp,'visible', 'off'); 
hI = plot(Time(indVI), vol(indVI),'r+');
hE = plot(Time(indVE), vol(indVE),'ro');
hI2 = plot(Time(indVI), vol2(indVI),'b+');
hE2 = plot(Time(indVE), vol2(indVE),'bo');
axis([Time(1) Time(end) min([vol;vol2]) max([vol;vol2])])

% Correct volumes
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
    subplot(6,1,4)
    plot(Time(indVI), vol(indVI),'r+');
    plot(Time(indVE), vol(indVE),'ro');
    plot(Time(indVI), vol2(indVI),'b+');
    plot(Time(indVE), vol2(indVE),'bo');
end


% Calc vdote
VdotE = [];
for n = 1:length(indVI)-1
    VdotE(n) = (TidalVol(n))/(Time(indVE(n+1))-Time(indVE(n)));
end
VdotE = VdotE.*60;
a = find(VdotE<0);
VdotE(a) = 0;
VdotE = VdotE';
clear a n;

TidalVol2 = vol2(indVI)-vol2(indVE);
VdotE2 = [];
for n = 1:length(indVI)-1
    VdotE2(n) = (TidalVol2(n))/(Time(indVE(n+1))-Time(indVE(n)));
end
VdotE2 = VdotE2.*60;
a = find(VdotE2<0);
VdotE2(a) = 0;
VdotE2 = VdotE2';
clear a n;

%Calc Vup, Vdown, LG
aind = find(Time(indVI)<tup,1,'last');
aind = aind+1;

%Calculate overshoot 2 different ways
Vup = mean(VdotE(aind:aind+1)); %Calculate LG using an average of breaths 1 and 2 from overshoot (Lisa's VT leak removal)=LG
Vup_1br =VdotE(aind); %Calculate LG using an breaths 1 from overshoot (Lisa's VT leak removal)=LG3
Vdown = mean(VdotE(aind-5:aind-1));

Vup2 = mean(VdotE2(aind:aind+1)); %Calculate LG using an average of breaths 1 and 2 from overshoot (Andrews VT leak removal)=LG2
Vup2_1br = VdotE2(aind); %Calculate LG using an breaths 1 from overshoot (Andrews VT leak removal)=LG4
Vdown2 = mean(VdotE2(aind-5:aind-1));

%Calculate LG using 2 breath average
if Veupnea-Vdown < 0.5 || Vup-Veupnea < 0.5
    LG = NaN;
    LG2 = NaN;
else
    LG = (Vup-Veupnea)/(Veupnea-Vdown);
    LG2 = (Vup2-Veupnea)/(Veupnea-Vdown2);
end

%Calculate LG using first breath only
if Veupnea-Vdown < 0.5 || Vup_1br-Veupnea < 0.5
    LG3 = NaN;
    LG4 = NaN;
else
    LG3 = (Vup_1br-Veupnea)/(Veupnea-Vdown);
    LG4 = (Vup2_1br-Veupnea)/(Veupnea-Vdown2);
end

%Calculate CPAP level from which CPAP was determined
cpap = Pmask(indVE);
Cdown = mean(cpap(aind-5:aind-1));

%Plot VdotE, Vup, Vdown, Veupnea
subplot(6,1,4)
plot(Time(indVI(aind:aind+1)), vol(indVI(aind:aind+1)),'g+');
plot(Time(indVE(aind:aind+1)), vol(indVE(aind:aind+1)),'go');
plot(Time(indVE(aind):indVI(aind+1)), vol(indVE(aind):indVI(aind+1)), 'g')
plot(Time(indVI(aind-5:aind-1)), vol(indVI(aind-5:aind-1)),'b+');
plot(Time(indVE(aind-5:aind-1)), vol(indVE(aind-5:aind-1)),'bo');
plot(Time(indVE(aind-5):indVI(aind-1)), vol(indVE(aind-5):indVI(aind-1)), 'b')

subplot(6,1,6);
plot(Time(indVI(1:end-1)), VdotE,'k.-')
hold on
ylabel('VdotE (L/min)')
plot([Time(indVI(aind-5)),Time(indVI(aind-1))], [Vdown Vdown],'b');
plot([Time(indVI(aind)),Time(indVI(aind+1))], [Vup Vup],'g');
plot([Time(indVI(1)),Time(indVI(end))], [Veupnea Veupnea],'k--');
 text(double(Time(indVI(aind-5)))-1,mean(VdotE)+1.5, strcat('Veup-Vdown =', num2str(Veupnea-Vdown,3)))
  text(double(Time(indVI(aind+1)))-1,mean(VdotE)+1.5, strcat('Vup-Veup =', num2str(Vup-Veupnea,3)))
  if isnan(LG)
       text(double(Time(indVI(aind-1)))-1,Veupnea+1.5, 'LG = nan')
  else
   text(double(Time(indVI(aind-2)))-1,Veupnea+1.5, strcat('LG (2br avg) =', num2str(LG,3)))
  end
end