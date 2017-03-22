function [mxInd,mnInd] = SelectPepi(Time, Pepi, pfilt, Vdot, t, stage, evts, num); 
%Lisa M Campana, PhD
%December 2012

%Inputs: Time, Pepi, pfilt, Vdot vectors of raw data. t = time picked by
%user. stage & evts = staging, num = number of breaths you wish to identify
%by hand
%Outputs: new max and min indexs

%Plot Figure
fh = figure;
ax(1) = subplot(3,1,1); %Sleep stage
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
plot([t t], [0 6], 'r')
ylabel('Sleep Stage')
xlim([Time(1) Time(end)])
ylim([0 6])
title('Select max Pepi then min Pepi')

ax(2) = subplot(3,1,2);
hold on; 
plot(Time, Pepi,'b'); hold on 
plot(Time, pfilt,'k', 'linewidth', 1)
plot([t t], [min(Pepi) max(Pepi)],'r')

ax(3) = subplot(3,1,3);
hold on ;
plot(Time, Vdot,'k')
plot([t t], [min(Vdot) max(Vdot)],'r')
plot([Time(1) Time(end)], [0 0],'k')

linkaxes(ax,'x')
scrsz = get(0,'ScreenSize');
set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])

subplot(3,1,2)

%Pick max then min
for i = 1:num
    [x,y] = ginput(1);
    mxInd(i) = find(Time>=x,1,'first');
    mx(i) = Pepi(mxInd(i)); 
    plot(Time(mxInd), mx,'r+')

    [x,y] = ginput(1);
    mnInd(i) = find(Time>=x,1,'first');
    mn(i) = Pepi(mnInd(i));
    plot(Time(mnInd), mn,'ro')

end

close(fh)
