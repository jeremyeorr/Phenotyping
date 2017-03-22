 a = find(Time.Vdot>=Data.Traits2EvtTimeStart(indxGG(i),1) & Time.Vdot<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        a2 = find(Time.Epochs>=Data.Traits2EvtTimeStart(indxGG(i),1)-30 & Time.Epochs<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        a3 = find(Data.EvtTimeEnd>=Data.Traits2EvtTimeStart(indxGG(i),1) & Data.EvtTimeStart<=Data.Traits2EvtTimeEnd(indxGG(i),1));
        if isempty(Data.Pepi)
            peptemp = [];
        else
            peptemp= Data.Pepi(a);
        end
        tempGG = GGscaled(ceil(Data.Traits2EvtTimeStart(indxGG(i),1)/Fs.GG): floor(Data.Traits2EvtTimeEnd(indxGG(i),1)/Fs.GG));
        timeGG = [ceil(Data.Traits2EvtTimeStart(indxGG(i),1)/Fs.GG)*Fs.GG: Fs.GG: floor(Data.Traits2EvtTimeEnd(indxGG(i),1)/Fs.GG)*Fs.GG]';
        [Breaths(i)] = GGanalysis(Fs.Vdot, Time.Vdot(a), Data.Vdot(a), Data.Pmask(a), [Time.Epochs(a2), Data.Epochs(a2)], [Data.EvtTimeStart(a3), Data.EvtTimeEnd(a3)], peptemp, timeGG, tempGG);
        subplot(7,1,1)
        title({['Number ', num2str(i), ' of ', num2str(length(indxGG))]; 'Press Enter to continue'})
        pause;
        close all