function [scale,offset] = CalibrateGG(Fs, Time, GG, w, p)
if nargin<5
    p = 0.1; % window over which to remove DC offset (must be multiple of Fs)
    if nargin<4
        w = 0.2; % window over which to smooth GG (must be multiple of Fs)
    end
end

%DC remove, calculates average value of GG between t-w to t+w
if rem(w,Fs)~=0
    txt = ['Window is not divisible by ',num2str(Fs), ' Enter new w:'];
    w = inputdlg(txt,'',1, {'0.1'}); 
    w = str2num(cell2mat(w));
end
%create a moving window of length 2*w/Fs +1
MovWin = ones((2*w/Fs)+1,1)./((2*w/Fs)+1);
MovAvg = conv(GG, MovWin);
GGremove = GG - MovAvg((w/Fs+1):end-w/Fs); 

%Rectify GG
GGrectify = abs(GGremove);                                                                         

%Smooth GG with period p
SmoothWin = ones((2*p/Fs)+1,1)./((2*p/Fs)+1);
SmoothAvg = conv(GGrectify, SmoothWin);
GGsmooth = SmoothAvg((p/Fs+1):end-p/Fs);

figure;
hold all
ax(1) = subplot(4,1,1);
plot(Time, GG, 'k'); 
ylabel('GG Volts');
title('Press enter to continue')
ax(2) = subplot(4,1,2);
plot(Time, GGremove, 'k'); 
ylabel('GG DCremove');

ax(3) = subplot(4,1,3);
plot(Time, GGrectify, 'k'); 
ylabel('GG rectify');

ax(4) = subplot(4,1,4);
plot(Time, GGsmooth, 'k'); 
ylabel('GG smooth');

linkaxes(ax,'x');
axis([Time(1) Time(end) min(GGsmooth) max(GGsmooth)+0.1])

pause;
close all

fh = figure;
plot(Time, GGsmooth,'k')
hold on; 
set(fh,'Toolbar', 'Figure')
scrsz = get(0,'ScreenSize');
x = [];

set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
pbh = uicontrol(fh,'Style','pushbutton','String','Select Electrical Zero',...
                'Position',[20 750 150 50], 'UserData', x, 'Callback', @SelectEE);
pbh2 = uicontrol(fh,'Style','pushbutton','String','No Electrical Zero',...
                'Position',[20 650 150 50]);            

[GG100, indGG100] = max(GGsmooth); 
plot(Time(indGG100), GG100,'r+')

% set(pbh, 'Callback', 'uiresume; [x] = ginput(2); set(pbh, ''UserData'', x)');  
set(pbh2, 'Callback', 'uiresume;');
set(pbh, 'Callback', {@SelectEE, x});
uiwait
x = get(pbh, 'UserData');
if isempty(x)
    GG0 = 0; 
else
    a = find(Time>=x(1) & Time<=x(2));
    GG0 = mean(GGsmooth(a));
    plot([x(1) x(2)], [GG0 GG0],'r');
end

scale = 100/(GG100-GG0);
offset = 100-scale*GG100;
GGscale = GGsmooth*scale+offset; 
figure; 
plot(Time, GGscale,'k')
title('Scaled GG')
end

function SelectEE(hobj, eventdata, handles)
uiresume;
x = get(hobj, 'UserData');
x = ginput(2); 
set(hobj, 'UserData', x);
end
