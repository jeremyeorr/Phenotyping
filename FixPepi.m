function [mn, peakind] = FixPepi(Time, Pepi, Vdot, mn, peakind)
%Brad Edwards, PhD
%Inputs: Time vector, Pepi vector, and the indices of start inspiration
%(peakind) and peak pressure (mn)
%Outputs: new indcies of start inspiration and and peak pressure on user input

peakind=peakind';
% Filter Pepi signal
dt=(Time(end)-Time(1))/(length(Time)-1);
backup_mn=mn;

fh = figure;
while 1
% plot current volumes with start insp and peak insp
% First you have the option to delete points
subplot(2,1,1), plot(Time, Pepi); hold on
hmin = plot(Time(mn(:,1)), Pepi(mn(:,1)),'r.','MarkerSize',20);
hmax = plot(Time(peakind), Pepi(peakind),'g.','MarkerSize',20);
title('Left click deletes BAD Pepi Nadir data point (also deletes paired insp value), right click to finish')
hold off;
subplot(2,1,2), plot(Time, Vdot,'b',[Time(1), Time(end)], [0,0],'r');
scrsz = get(0,'ScreenSize');
set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);

   [x,y,button]=ginput(1);
   if button==3
       break
       close(fh)
   end
   TimeFromClickToPeakPepi = Time(peakind)-x;
   [~,i]=min(abs(TimeFromClickToPeakPepi));
   mn(i,:)=[];
   peakind(i)=[];
    
end

while 1
% This code moves the insp value to where you want it
subplot(2,1,1), plot(Time, Pepi); hold on
hmin = plot(Time(mn(:,1)), Pepi(mn(:,1)),'r.','MarkerSize',20);
hmax = plot(Time(peakind), Pepi(peakind),'g.','MarkerSize',20);
title('Left click to move peak Pepi to x location, right click to finish')
hold off;
subplot(2,1,2), plot(Time, Vdot,'b',[Time(1), Time(end)], [0,0],'r');
scrsz = get(0,'ScreenSize');
set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);

   [x,y,button]=ginput(1)
   if button==3
       break
       close(fh)
   end
   TimeFromClickToPeakPepi = Time(peakind)-x;
   [y,i]=min(abs(TimeFromClickToPeakPepi));
   oldpeakind=peakind(i)';
   newpeakind=find(Time>x,1);
   peakind(i)=newpeakind;
    
end


