function [scale,offset] = CalibrateGG(Fs, Time, GGsmooth)
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
