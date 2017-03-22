function pcrit_model_run
% close all
global Pmaskdata VdotMaxdata VdotEdata
if 0 %make up some data
Pmaskdata  =[0 3 4 4.5 5 5 5 6 6 ] 
VdotMaxdata=[0 0.3 0 0.2 0.4 0.4 0.5 0.6 0.6]
VdotEdata=VdotMaxdata*10;
end
fh=figure;
set(fh,'color',[1 1 1]);
ax1(1)=subplot(2,2,1); 
pcrit_model_run_once(0,1)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('Vpeak','fontname','arial narrow','fontsize',12);
ax1(2)=subplot(2,2,2); 
pcrit_model_run_once(1,1)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('Vpeak','fontname','arial narrow','fontsize',12);
ax2(1)=subplot(2,2,3); 
pcrit_model_run_once(0,0)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
ax2(2)=subplot(2,2,4); 
pcrit_model_run_once(1,0)
set(ax1,'fontname','arial narrow','fontsize',12)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
linkaxes(ax1);linkaxes(ax2);

function [VcritSEM,Vcritvalue,PcritSEM,Pcritvalue,PVSlope]=pcrit_model_run_once(exclude_zeroflowdata,usepeakflownotVE)
global Pmaskdata VdotMaxdata VdotEdata

Pmask=Pmaskdata;
if usepeakflownotVE
Vflow=VdotMaxdata;
else
Vflow=VdotEdata;    
end

if exclude_zeroflowdata
Pmask(Vflow==0)=[];
Vflow(Vflow==0)=[];
end

plot(Pmask,Vflow,'.','markersize',16); hold('on'); box('off');

lsqoptions=optimset('display','off','maxiter',500,'tolx',10E-3);
Pcrit_lowerlimit=-20;
Pcrit_upperlimit=+20;
lower=[0 Pcrit_lowerlimit]
upper=[max(Vflow)*10 Pcrit_upperlimit]
parameters=[1 Pcrit_lowerlimit];
modeloption=1; %not used just yet...
pcrit_model(parameters,Pmask,1)
for i=1:5
    [parameters,~,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
end
%have a good start point - but can we do better:
lower=[parameters(1)/3 parameters(2)-5]
upper=[parameters(1)*2 parameters(2)+5]
for i=1:50
    parameters_start=randn*(upper-lower)+lower;
    [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters_start,Pmask,Vflow,lower,upper,lsqoptions);
end
parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[~,i]=min(Fres); parameters=parametersi(i,:);

%have a good better point - but can we do even better again:
lower=[parameters(1)/1.5 parameters(2)-1]
upper=[parameters(1)*1.5 parameters(2)+1]
for i=1:100
    parameters_start=randn*(upper-lower)+lower;
    [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters_start,Pmask,Vflow,lower,upper,lsqoptions);
end
parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[~,i]=min(Fres); parameters=parametersi(i,:);

Pmaskline=min(Pmask):0.01:max(Pmask);
%[modelVdotMax] = pcrit_model(parameters,Pmaskline,modeloption);
%confidence intervals
[modelVdotMax,delta] = nlpredci(@(parameters,Pmaskline) pcrit_model(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN);

upperSEM=modelVdotMax+delta'/1.96;
lowerSEM=modelVdotMax-delta'/1.96;

filly=[upperSEM fliplr(lowerSEM)]
fillx=[Pmaskline fliplr(Pmaskline)]
fill(fillx,filly,[0.8 0.2 0.2],'linestyle','none','facealpha',0.5); 
plot(Pmask,Vflow,'.','markersize',16); hold('on'); box('off');
plot(Pmaskline,modelVdotMax,'r'); 

CI_parameters = NLPARCI(parameters,RESIDUAL,'jacobian',JACOBIAN);
Pcritvalue=parameters(2);
PcritSEM=(CI_parameters(2,2)-Pcritvalue)/1.96
%ci = nlparci(beta,resid,'jacobian',J)

title(['Pcrit=' num2str(Pcritvalue,4) setstr(177) num2str(PcritSEM,3)])

PVSlope=parameters(1);
Vcritvalue=-PVSlope*Pcritvalue

Vcrit_lowerlimit=Vcritvalue-1;
Vcrit_upperlimit=Vcritvalue+1
lower=[PVSlope/2 Vcrit_lowerlimit]
upper=[PVSlope*2 Vcrit_upperlimit]
parameters2=[PVSlope Vcritvalue];
modeloption=2; %not used just yet...
pcrit_model(parameters2,Pmask,2)
for i=1:50
    [parameters2,~,RESIDUAL2,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN2] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
end
[modelVdotMax2,delta2] = nlpredci(@(parameters2,Pmaskline) pcrit_model(parameters2,Pmaskline,modeloption),Pmaskline,parameters2,RESIDUAL2,'Jacobian',JACOBIAN2);

if 0 %check that models are identical
filly=[modelVdotMax2+delta2' fliplr(modelVdotMax2-delta2')]
fillx=[Pmaskline fliplr(Pmaskline)]
fill(fillx,filly,[0.2 0.8 0.2],'linestyle','none','facealpha',0.5); 
plot(Pmaskline,modelVdotMax2,'b:'); 
end

plot(0,Vcritvalue,'bo'); 



CI_parameters2 = NLPARCI(parameters2,RESIDUAL2,'jacobian',JACOBIAN2);

VcritSEM=abs((CI_parameters2(2,2)-parameters2(2))/1.96)

if Pcritvalue>0
I=find(abs(upperSEM)>0);
x=[0 Pmaskline(I)];
yupper=[Vcritvalue+VcritSEM upperSEM(I)];
ylower=[Vcritvalue-VcritSEM lowerSEM(I)];
Pmaskline2=0:0.01:Pmaskline(I);
upperSEMinterp=interp1(x,yupper,Pmaskline2,'spline')
lowerSEMinterp=interp1(x,ylower,Pmaskline2,'spline')

plot([0 Pmaskline(I)],[Vcritvalue modelVdotMax(I)],'b:'); 
plot(Pmaskline2,lowerSEMinterp,'b:'); 
plot(Pmaskline2,upperSEMinterp,'b:'); 
end

title(['Pcrit=' num2str(parameters(2),4) setstr(177) num2str(PcritSEM,3) '; ' 'Vcrit=' num2str(parameters2(2),4) setstr(177) num2str(VcritSEM,3) ] )
xlim([min([Pmask,0])-0.2,max(Pmask)])

herr=errorbar(0,Vcritvalue,VcritSEM,'-')


function y = pcrit_model(x,xdata,modeloption)

if modeloption==1
y=x(1)*(xdata-x(2)); %x(2) is Pcrit
y(y<0)=0;
end
if modeloption==2
y=x(1)*(xdata)+x(2); %x(2) is Vcrit
y(y<0)=0;
end
