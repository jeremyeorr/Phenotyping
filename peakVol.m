function [TidalVol,indVI,indVE] = peakVol(time,V,thresh)
%Lisa Campana PhD
%peak detection algorithm for volume signal

if nargin<3
    thresh = 0.1;
end
[max_list, min_list] = peakdet(V,thresh);

[indVI, indVE] = CheckVind(max_list, min_list);

TidalVol = V(indVI)-V(indVE);  




