function [ind_vmax, ind_vmin] = CheckVind(max_list, min_list)
%%make sure order is min-max-...-min-max (end expire-end inspire-)
if max_list(1,1)<min_list(1,1)  %(First point is end inspire point)
    ind_vmax = max_list(2:end,1);
    if length(max_list(:,1))==length(min_list(:,1)); %(max-min is order)
        ind_vmin = min_list(1:end-1,1);
    elseif length(max_list(:,1))>length(min_list(:,1)) %(max-min-max is order)
        ind_vmin = min_list(1:end,1);
    end
else
    if length(max_list(:,1))==length(min_list(:,1));
        ind_vmax = max_list(1:end,1);
        ind_vmin = min_list(1:end,1);
    else
        ind_vmin = min_list(1:end-1,1);
        ind_vmax = max_list(1:end,1);
    end
end