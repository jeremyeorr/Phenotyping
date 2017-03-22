function [ta, uag]=SummarizeData(subid, night, MR, Veupnea, Varousal, Vactive, V0, Ug, LG, Cpap, ArPepi)
% 
nremahi = inputdlg('Enter NREM AHI: ');
remahi = inputdlg('Enter REM AHI: ');

% %only use Vactive when CPAP is lowest
% if length(Cpap)>3
%     tempCpap = round(Cpap); 
%     [minCpap] = min(tempCpap); 
%     minInd = find(tempCpap==minCpap); 
%     if length(minInd)>3
%         lowCdown = Cpap(minInd);
%         Vactlow = mean(Vactive(minInd));
%     else
%         temp = find(tempCpap>minCpap);
%         minInd = find(tempCpap==minCpap | tempCpap==min(tempCpap(temp)));
%         lowCdown = Cpap(minInd);
%     end
% else
%     minInd = [1:length(Cpap)];
%     lowCdown = Cpap;
% end
% VactiveLow = nanmean(Vactive(minInd));
% UgLow = nanmean(Ug(minInd)); 

% loop gain line
x = Veupnea:.01:50;
b = Veupnea - (1/LG)*Veupnea;
y = (1/LG)*x + b;

% Arousal threshold line
ta = LG*(Varousal-b);

% upper airway gain line
uag = (Vactive-V0)/(ta-Veupnea);
b1 = V0 - uag*Veupnea;
y1 = uag*x + b1;

% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'FontSize',12);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 20]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 10]);
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(Veupnea,Veupnea,Veupnea,Varousal,Veupnea,Vactive,Veupnea,V0,'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',15,...
    'Marker','.',...
    'LineStyle','none',...
    'Color',[0 0 0],...
    'Parent',axes1);
set(plot1(2),'DisplayName','Varousal vs Veupnea');
set(plot1(3),'DisplayName','Vactive vs Veupnea');
set(plot1(4),'DisplayName','V0 vs Veupnea');

% Create xlabel
xlabel('Ventilatory drive (L/min)','FontSize',14);

% Create ylabel
ylabel('Ventilation (L/min)','FontSize',14);

% Create title

title({['Subject: ', subid, ', ', night]; ['NREM AHI ', cell2mat(nremahi)];['REM AHI ', cell2mat(remahi)]},'FontSize',14);

% Create multiple lines using matrix input to plot
plot2 = plot(x,y,x,y1,'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'LineWidth',1,...
    'Color',[0 0 0]);
set(plot2(1),'DisplayName','y vs x');
set(plot2(2),'DisplayName','y1 vs x');

% Create plot
plot([ta ta],[-10 50],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'LineWidth',1,...
    'DisplayName','[-10 50] vs [ta ta]',...
    'Color',[0 0 0]);

% Create multiple lines using matrix input to plot
plot(ta,Varousal,ta,Vactive,'Visible','off','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',15,...
    'Marker','.',...
    'LineStyle','none',...
    'DisplayName','Vactive vs ta',...
    'Color',[0 0 0]);

% Create multiple lines using matrix input to plot
plot3 = plot([Veupnea ta-.4],[Varousal Varousal],[Veupnea ta-.4],[Vactive Vactive],'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'LineStyle',':',...
    'Color',[0 0 0]);
set(plot3(1),'DisplayName','[Varousal Varousal] vs [Veupnea ta-.1]');
set(plot3(2),'DisplayName','[Vactive Vactive] vs [Veupnea ta-.1]');

% Create multiple lines using matrix input to plot
plot4 = plot(ta-.4,Varousal,ta-.4,Vactive,'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',5,...
    'Marker','>',...
    'Color',[0 0 0]);
set(plot4(1),'DisplayName','Varousal vs ta-.5');
set(plot4(2),'DisplayName','Vactive vs ta-.5');

% Create text
text('Parent',axes1,'String',['V_e_u_p_n_e_a = ' num2str(round(Veupnea*10)/10)],...
    'Position',[Veupnea-4.1 Veupnea 0],...
    'Color',[0 0 0]);

% Create text
text('Parent',axes1,'String',['V_a_r_o_u_s_a_l = ' num2str(round(Varousal*10)/10)],...
    'Position',[Veupnea-4.1 Varousal 0],...
    'Color',[0 0 0]);

% Create text
text('Parent',axes1,'String',['active V_0 = ' num2str(round(Vactive*10)/10)],...
    'Position',[Veupnea-4.4 Vactive 0],...
    'Color',[0 0 0]);

if V0 < 1
    text('Parent',axes1,'String',['passive V_0 = ' num2str(round(V0*10)/10)],...
        'Position',[Veupnea-4.9 V0+.4 0],...
        'Color',[0 0 0]);
else
    text('Parent',axes1,'String',['passive V_0 = ' num2str(round(V0*10)/10)],...
        'Position',[Veupnea-4.9 V0 0],...
        'Color',[0 0 0]);
end

% Create text
text('Parent',axes1,'String',['ArThr (L/min) = ' num2str(round(ta*10)/10)],...
    'Position',[ta+.2 9 0],...
    'Color',[0 0 0]);
text('Parent',axes1,'String',['ArThr (cmH2O) = ' num2str(round(ArPepi))],...
    'Position',[ta+.2 8 0],...
    'Color',[0 0 0]);
text('Parent',axes1,'String',['MR (%max/cmH2O) = ' num2str(sprintf('%.2f',MR))],...
    'Position',[ta+.2 7 0],...
    'Color',[0 0 0]);
% Create text
text('Parent',axes1,'String',['LG = ' num2str(round(LG*10)/10)],...
    'Position',[ta+2 Veupnea 0],...
    'Color',[0 0 0]);

% Create text
if V0 < 1
text('Parent',axes1,'String',['UAG = ' num2str(round(uag*10)/10)],...
    'Position',[ta+2 V0+0.3 0],...
    'Color',[0 0 0]);
else
text('Parent',axes1,'String',['UAG = ' num2str(round(uag*10)/10)],...
    'Position',[ta+2 V0 0],...
    'Color',[0 0 0]);
end

saveas(gcf,'Model');

