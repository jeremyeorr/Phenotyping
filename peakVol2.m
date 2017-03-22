function [indVI,indVE] = peakVol2(time,V,thresh)
%Lisa Campana PhD
%peak detection algorithm for volume signal

if nargin<3
    thresh = 0.1;
end
[indVE,indVI] = peakdet(-V,thresh);
indVI(:,2)=[];
indVE(:,2)=[];
if indVE(1)>indVI(1)
    indVE(1)=[];
end





