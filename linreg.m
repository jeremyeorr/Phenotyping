function [fit]=linreg(x,y,sp,xlab,ylab)

fit = polyfit(x,y,1);
xtemp1 = [min(x):.1:max(x)];
ytemp1 = xtemp1.*fit(1) + fit(2);

%calculating r-squared
yfit = polyval(fit,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
fit(1,3) = 1 - SSresid/SStotal;

subplot(3,2,sp);
hold on;
cc=rand(1,3);
plot(x,y,'r.','MarkerSize',22,'color',cc);
plot(xtemp1,ytemp1,'color',cc);
title(['slope = ',num2str(fit(1), '%0.2f'), ' intercept = ', num2str(fit(2), '%0.2f'), ' r-sq = ', num2str(fit(3), '%0.2f')])
xlabel(xlab)
ylabel(ylab)