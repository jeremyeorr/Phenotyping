clear all
close all

x = [.008:.008:120];
Region1 = [1:5000];
Region2 = [5001:10000];
Region3 = [10001:15000];

z1 = 0.1;
z2 = 0.3;
z3 = 0.1;

er = rand(1,15000);

y(Region1) = 1.*sin(x(Region1))+z1+er(Region1);
y(Region2) = 2.*sin(x(Region2))+z2+er(Region2);
y(Region3) = 1.*sin(x(Region3))+z3+er(Region3);
plot(x,y)

v = cumtrapz(y)/125;
v1 = detrend(v, 'linear', [Region2(1), Region3(1)]);

p1 = polyfit(x(Region1),v(Region1),1);
p2 = polyfit(x(Region2),v(Region2),1);
p3 = polyfit(x(Region3),v(Region3),1);

y2(Region1) = y(Region1)-p1(1);
y2(Region2) = y(Region2)-p2(1);
y2(Region3) = y(Region3)-p3(1);

y3(Region1) = y(Region1)-mean(y(Region1));
y3(Region2) = y(Region2)-mean(y(Region2));
y3(Region3) = y(Region3)-mean(y(Region3));

v2 = cumtrapz(y2)/125; 
v3 = cumtrapz(y3)/125;
figure
plot(x,v,'k'); hold on;
plot(x,v1,'r--')
plot(x,v2,'g')
plot(x,v3,'b--')
[TidalVol,indVI,indVE] = peakVol(x,v1,0.1);

[TidalVol2,indVI2,indVE2] = peakVol(x,v2,0.1);

[TidalVol3,indVI3,indVE3] = peakVol(x,v3,0.1);
TidalVol = TidalVol';
TidalVol2 = TidalVol2';
TidalVol3 = TidalVol3';
d = TidalVol-TidalVol2;
d2 = TidalVol2-TidalVol3;