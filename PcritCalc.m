function [PmaskDrop, VdotMax, VdotE, TimeBreath, Pcrit Vcrit] = PcritCalc(Fs, Time, Vdot, Pmask, stage, evts, Threshold)

thresh=0.01;
% Integrate the raw flow signal to get volume
vol =cumtrapz(Vdot*Fs);
%detrend volume signal to remove leak
%vol = detrend(vol, 'linear');
%Calculate end inspiration and end expiration points, tidal volumes
[TidalVol,indVI,indVE] = peakVol(Time,vol,0.1);

numplots = 4;
 
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

ax(3) = subplot(numplots,1,3); %Volume
plot(Time, vol, 'k')
hold on;
hI = plot(Time(indVI), vol(indVI),'r+');
hE = plot(Time(indVE), vol(indVE),'ro');
ylabel('Volume (L)')

ax(4) = subplot(numplots,1,4); %Flow
plot(Time, Vdot,'k')
hold on
plot([Time(1), Time(end)], [0,0],'r')
ylabel('Flow (L/s)')

linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
% set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])
set(fh,'Toolbar', 'Figure')

%Find pressure drops based on maximum and minimum derivative of Pmask signal
%a = find(Pmask<0.75*max(Pmask)); Lisa original code (I lowered due to
%errors in detecting mask pressure during drops...
a = find(Pmask<Threshold*max(Pmask)); 
b = find(diff(a)>2/Fs);
b(end+1) = length(a);
indstart = 1;
subplot(numplots,1,2)
plot([Time(1) Time(end)], [0.75*max(Pmask), 0.75*max(Pmask)],'g')
for i = 1:length(b)
    indx(i,1) = a(indstart)-5/Fs;
    indx(i,2) = a(b(i))+5/Fs; 
    indstart = b(i)+1;
    subplot(numplots,1,2)
    plot(Time(indx(i,1):indx(i,2)), Pmask(indx(i,1):indx(i,2)), 'r')
end
VdotE = zeros(3,size(indx,1));
VdotMax = zeros(3,size(indx,1));
PmaskDrop = zeros(3,size(indx,1));
TimeBreath = zeros(3,size(indx,1));
subplot(numplots,1,1)
title('Press any key to continue')
pause;
for i = 1:size(indx,1)
    tempTime = Time(indx(i,1):indx(i,2)); 
    tempVdot = Vdot(indx(i,1):indx(i,2));
    tempPmask = Pmask(indx(i,1):indx(i,2));
    tempVol = vol(indx(i,1):indx(i,2));
    
    %reframe figure
    subplot(numplots,1,2)
    axis([tempTime(1) tempTime(end) min(tempPmask) max(tempPmask)])

    detrendeachdrop=1;
    
if detrendeachdrop    
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
Region1 = find(tempTime<=tdown);
Region2 = find(tempTime>tdown & tempTime<=tup);
Region3 = find(tempTime>tup);

%detrend volume signal to remove leak


if 0
tempVol_predetrend=tempVol;
tempVol = detrend(tempVol, 'linear', [Region2(1), Region3(1)]); %overwrite
ydetrendVol = tempVol_predetrend-tempVol;
ydetrendVdot=gradient(ydetrendVol)/Fs;
tempVdot_predetrend=tempVdot;
if 0
    tempVdot = detrend(tempVdot, 'linear', [Region2(1), Region3(1)]); %overwrite
else
    tempVdot=tempVdot-ydetrendVdot;
end
end

leakRegion2=mean(tempVdot(Region2));
tempVdot_predetrend=tempVdot;
tempVdot=tempVdot-leakRegion2;
tempVol_predetrend = tempVol;
tempVol = cumtrapz(tempVdot)*Fs;

figure(200)
plot(tempTime,tempVol_predetrend,tempTime,tempVol,tempTime,tempVol_predetrend-tempVol)
figure(201)
plot(tempTime,tempVdot_predetrend,tempTime,tempVdot,tempTime,tempVdot_predetrend-tempVdot)
figure(fh)
end


    [VtTemp,indVItemp,indVEtemp] = peakVol(tempTime, tempVol,thresh);
    
    
    %continue reframing / updating figure
    subplot(numplots, 1, 3)
    
    axis([tempTime(1) tempTime(end) min(tempVol) max(tempVol)])
    set(hI, 'Visible', 'off')
    set(hE, 'Visible', 'off')
    hItemp = plot(tempTime(indVItemp), tempVol(indVItemp),'r+');
    hEtemp = plot(tempTime(indVEtemp), tempVol(indVEtemp),'ro');
    htempVol = plot(tempTime,tempVol,'r');
    
    
    subplot(numplots, 1, 4)
    htempVdot = plot(tempTime,tempVdot,'r');
    axis([tempTime(1) tempTime(end) min(Vdot) max(tempVdot)])
    
    
    [zf] = listdlg('PromptString','Is there zero flow?',...
        'SelectionMode','single', 'ListSize', [150 50],...
        'ListString',{'Yes','No'});
    if zf==1
        subplot(numplots,1,1);
        title('Select start and end of zero flow');
        [x,y] = ginput(2);
        VdotE(1,i) = 0; 
        VdotMax(1,i) = 0;
        indx0 = find(tempTime>=x(1,1) & tempTime<x(2,1));
        PmaskDrop(1,i) = mean(tempPmask(indx0));
        TimeBreath(1,i) = x(1,1);
    elseif zf==2
        UseBreaths = listdlg('PromptString','Use Any Breaths?',...
            'SelectionMode','single', 'ListSize', [150 50],...
            'ListString',{'Yes','No'});
        if UseBreaths==1
            RedoVol = listdlg('PromptString','Fix any volumes?',...
                'SelectionMode','single', 'ListSize', [150 50],...
                'ListString',{'Yes','No'});
            if RedoVol==1
                [indVItemp,indVEtemp] = FixVolumes(tempTime, tempVol, indVItemp, indVEtemp);
                set(hItemp, 'Visible', 'off')
                set(hEtemp, 'Visible', 'off')
                subplot(numplots,1,3)
                plot(tempTime(indVItemp), tempVol(indVItemp),'r+');
                plot(tempTime(indVEtemp), tempVol(indVEtemp),'ro');
            end
            subplot(numplots,1,1); 
            title('Select breaths: click to the left of first breath and to the right of the last breath you wish to analyze');
            [x2,y2] = ginput(2);
            breaths = find(tempTime(indVItemp)>=x2(1,1) & tempTime(indVItemp)<=x2(2,1));
            for n = 1:length(breaths)
                VdotE(n,i) = 60*(tempVol(indVItemp(breaths(n)))-tempVol(indVEtemp(breaths(n))))/(tempTime(indVEtemp(breaths(n)+1))-tempTime(indVEtemp(breaths(n))));
                [VdotMax(n,i),posVdotmax(n)] =max(tempVdot(indVEtemp(breaths(n)):indVItemp(breaths(n))));  
                posVdotmax(n) = posVdotmax(n)+ indVEtemp(breaths(n));
                PmaskDrop(n,i) = tempPmask(posVdotmax(n));
                TimeBreath(n,i) = tempTime(indVEtemp(breaths(n)));
            end
           subplot(numplots,1,3)
           plot(tempTime(indVItemp(breaths)), tempVol(indVItemp(breaths)),'g+')
           plot(tempTime(indVEtemp(breaths)), tempVol(indVEtemp(breaths)),'go')
           text(double(tempTime(1)), max(tempVol), strcat('VdotMax = ', num2str(mean(VdotMax(1:length(breaths),i)))));
           subplot(numplots,1,2);
           plot(tempTime(posVdotmax(1:length(breaths))), PmaskDrop(1:length(breaths),i),'g+');
           subplot(numplots,1,4);
           plot(tempTime(posVdotmax(1:length(breaths))), VdotMax(1:length(breaths),i),'g+');
           
        else
            VdotE(1,i) = nan;
            PmaskDrop(1,i) = nan;
            VdotMax(1,i) = nan;
            TimeBreath(1,i) = nan;
        end
    else
        VdotE = [];
        PmaskDrop = [];
        VdotMax = [];
        TimeBreath = [];
        close all
        return
    end
end
subplot(numplots,1,3)
axis([Time(1) Time(end), min(vol), max(vol)])
subplot(numplots,1,1)
title('Press any key to continue')
pause

fh2 = figure; 
b = find(VdotMax>0);
%b = find(VdotE>0);
cpap = PmaskDrop(b);
Vmax = VdotMax(b);
Ve = VdotE(b); 
 
fit1 = polyfit(cpap,Vmax,1);
Pcrit = -fit1(2)/fit1(1) 
pc = floor(Pcrit);
xtemp = [min(0,pc):max(max(PmaskDrop))+1]; 
ytemp = xtemp.*fit1(1)+fit1(2);

plot(cpap,Vmax,'b.','MarkerSize',22)
hold on; 
plot(xtemp, ytemp,'b--')
plot([xtemp(1) xtemp(end)], [0 0],'k')
ylabel('Vdot Max (L/s)')
xlabel('Pmask (cmH2O)')
title({['File time: ', num2str(round(Time(1)))]; ['Pcrit: ', num2str(Pcrit)]})
saveas(fh2, ['Pcrit at Time', num2str(round(Time(1)))]);


fh3 = figure;
fit2 = polyfit(cpap, Ve, 1);
Vcrit=fit2(2)
xtemp2 = [0:max(cpap)+1];
ytemp2 = xtemp2.*fit2(1)+fit2(2);
    
plot(cpap,Ve,'g.','MarkerSize',22)
hold on; 
plot(xtemp2, ytemp2,'g--')
plot([xtemp2(1) xtemp2(end)], [0 0],'k')
ylabel('VdotE (L/min)')
xlabel('Pmask (cmH2O)')
title({['File time: ', num2str(round(Time(1)))]; ['Vcrit: ', num2str(Vcrit)]})
saveas(fh3, ['Vcrit at Time', num2str(round(Time(1)))]);

pause
end

