function [indVI, indVE] = FixVolumes(Time, vol, indVI, indVE)
%Lisa M Campana, PhD
%December 2012
%Inputs: Time vector, volume vector, and the indices of end inspiration
%(indVI) and end expiration (indVE)
%Outputs: new indcies of end inspiration and expiration based on user input

% plot current volumes with end I and end E
fh = figure;
plot(Time, vol); hold on
hI = plot(Time(indVI), vol(indVI),'k+');
hE = plot(Time(indVE), vol(indVE),'ko');
scrsz = get(0,'ScreenSize');
set(fh,'OuterPosition',[1 1 scrsz(3) scrsz(4)])

%Loop until user right clicks twice
but = 1;
while but==1;
    title('Select around location where you would like to remove or add, right click twice when done');
    [x,y,but] = ginput(2);
    %added by SS Dec3, 2013 (denoted by lines with %%):
    %if x(1)<Time(indVI(1)) %%
        %aI=1;              %% 
        %aE=1;              %%
    %else                   %%
    aI = find(Time(indVI)>=x(1) & Time(indVI)<=x(2));
    aE = find(Time(indVE)>=x(1) & Time(indVE)<=x(2));
    %end                    %%
    if but==1
        if isempty(aI) & isempty(aE) %if no max or min exist, then assume you want to add a point
            if x(1)<Time(indVI(1))
                if indVI(1)>indVE(1) %leftmostpointisamaximum
                    bI=1; bE=2;
                else
                    bI=2; bE=1;
                end
            else
            bI = find(Time(indVI)<=x(1),1,'last'); %find whether last point is min or maximum Volume
            bE = find(Time(indVE)<=x(1),1,'last');
%             if bE>bI
%                 leftmostpointisamaximum=1
%             end
            end
            btemp = find(Time>=x(1) & Time<=x(2));
            if bE>bI %last point before is a minimum, so now look for a maximum
                [mx, indmx] = max(vol(btemp)); %find local maximum between points you selected
                tempVI = zeros(length(indVI)+1,1);
                tempVI(1:bI) = indVI(1:bI);
                tempVI(bI+1) = indmx+btemp(1)-1;
                tempVI(bI+2:end) = indVI(bI+1:end);
                plot(Time(tempVI(bI+1)), vol(tempVI(bI+1)),'r+');
                %ss workaround:
                tempVI=sort(tempVI);
                tempVE = indVE;
            else %last point before is a maximum, so now look for a min
                [mn, indmn] = min(vol(btemp)); %find local maximum between points you selected
                tempVE = zeros(length(indVE)+1,1);
                tempVE(1:bE) = indVE(1:bE);
                tempVE(bE+1) = indmn+btemp(1)-1;
                tempVE(bE+2:end) = indVE(bE+1:end);
                plot(Time(tempVE(bE+1)), vol(tempVE(bE+1)),'ro');
                %ss workaround:
                tempVE=sort(tempVE);
                tempVI = indVI;
            end
            clear bI bE btemp
        else
            indVI(aI) = [];
            indVE(aE) = [];
            set(hI, 'xdata', Time(indVI), 'ydata', vol(indVI));
            set(hE, 'xdata', Time(indVE), 'ydata', vol(indVE));
            tempVI = indVI;
            tempVE = indVE;
        end
        clear indVI indVE
        indVI = tempVI;
        indVE = tempVE;
    end
end
if length(indVI)~=length(indVE)
    errordlg('# of end expiration points does not equal # of end inspiration!!');
end
close(fh)   
end