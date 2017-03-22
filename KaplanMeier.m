%% Determine Arousal Threshold from Kaplan-Meier plot
function [fit]=KaplanMeier(ArTH Pepi ) 

for n = 1:R
    if Rate(n) == 2
        ArTh(n) = NaN;
        PAT(n) = NaN;
    end
end
    
Ma = max(VDRIVE,[],2)';   %maximum drive on each of the runs
ArTh(isnan(ArTh)) = 0;  %make NaNs zero

%if max drive is > AT, then use AT instead
for n = 1:R
    if ArTh(n) ~= 0
        Ma(n) = ArTh(n);
    end
end
c = (Ma>ArTh);   %get Boolean values when the maximum is > 0
[f,x,flo,fup] = ecdf(Ma,'censoring',c);

% Determine arousal threshold by fitting a line
if max(f) > 0.5
    i2 = find(f>=0.5,1,'first');
    i1 = i2-1;
    ff = [f(i1) f(i2)];
    xx = [x(i1) x(i2)];
    pp = polyfit(xx,ff,1);
    TA = (0.5 - pp(2))/pp(1);
else
    TA = max(x);
end; clear i2 i1 ff xx;

% Create figure
figure1 = figure;
 
% Create axes
axes1 = axes('Parent',figure1);
ylim([0 1]);
box('on');
hold('all');
 
% Create stairs
stairs1 = stairs(x,f,'LineWidth',2);
 
% Create stairs
stairs2 = stairs(...
  x,flo,...
  'Color',[1 0 0],...
  'LineWidth',2,...
  'LineStyle',':');
 
% Create stairs
stairs3 = stairs(...
  x,fup,...
  'Color',[1 0 0],...
  'LineWidth',2,...
  'LineStyle',':');
 
% Create scatter
scatter1 = scatter(...
  TA,0.5,...
  'MarkerEdgeColor',[0 0 1],...
  'MarkerFaceColor',[0 0 1]);
 
% Create xlabel
xlabel('Ventilatory Drive (L/min)','FontSize',14);
 
% Create ylabel
ylabel('Probability of Arousal','FontSize',14);
 
% Create title
title('Arousal Threshold','FontSize',14);

gtext(['Arousal Threshold = ' num2str(round(TA*100)/100) ' L/sec']);

clear figure1 stairs1 stairs2 stairs3 scatter1 axes1

saveas(gcf,'ArThr');
print -dmeta

pause(2);
close;