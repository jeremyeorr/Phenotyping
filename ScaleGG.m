function [t, g, GGscaled] = ScaleGG(Fs, GGsmooth, scale, offset)

clear g t
str = {'1000','2000','5000','10000', '20000', '50000', '100000'};
initalGain = inputdlg('Calibration gain (Hz): ', 'Gain', 1, {'5000'});
g(1) = str2double(cell2mat(initalGain)); 
t(1) = 0; 
 c = 1;
 d = 1; 
 while d == 1;
     c = c+1;
     temp = inputdlg({'Time of Gain Change(s)'}, 'time', 1, {'0'});
     t(c) = str2double(cell2mat(temp));
     if isempty(temp)
         d=0;
     else
        [opt] = listdlg('PromptString','Choose the new Gain (Hz)',...
                'SelectionMode','single',...
                'ListString',str);
        g(c) = str2double(cell2mat(str(opt)));
     end
 end
 
t(c) = length(GGsmooth)*Fs;
GGscaled = zeros(length(GGsmooth),1);
for i = 1:length(t)-1
    GainAmp = g(1)/g(i);
    GGscaled(t(i)/Fs+1:t(i+1)/Fs) = GGsmooth(t(i)/Fs+1:t(i+1)/Fs).*GainAmp*scale +offset;
end
plot(GGsmooth(1:100:end),'k'); hold on;
plot(GGscaled(1:100:end),'r')