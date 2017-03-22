function [t, g, GGscaled] = ScaleGG(Fs, GG, scale, offset, w,p)
if nargin<6
    p = 0.2; % window over which to remove DC offset (must be multiple of Fs)
    if nargin<5
        w = 0.1; % window over which to smooth GG (must be multiple of Fs)
    end
end

%DC remove, calculates average value of GG between t-w to t+w
if rem(w,Fs)~=0
    txt = ['Window is not divisible by ',num2str(Fs), ' Enter new w:'];
    w = inputdlg(txt,'',1, {'0.1'}); 
    w = str2num(cell2mat(w));
end

%create a moving window of length 2*w/Fs +1, then rectify
MovWin = ones((2*w/Fs)+1,1)./((2*w/Fs)+1);
MovAvg = conv(GG, MovWin);
GGrectify= abs(GG - MovAvg((w/Fs+1):end-w/Fs)); 
clear MovAvg
clear GG
                                                                        

%Smooth GG with period p
SmoothWin = ones((2*p/Fs)+1,1)./((2*p/Fs)+1);
GGsmooth = conv(GGrectify, SmoothWin);
GGsmooth = GGsmooth((p/Fs+1):end-p/Fs);
clear GGrectify

clear g t
str = {'5000','10000', '20000', '50000'};
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