function [comp_ind] = findICA(comp, thr_keep, thr_leave)
nICs = size(comp,1);
comp_ind =[];
for iIC=1:nICs
    [p, ind ]= max(comp(iIC, :));
    if (ind==1)
       comp_ind = [comp_ind iIC];
    elseif comp(iIC,1 )>thr_keep
        if (p<=thr_leave)
            comp_ind = [comp_ind iIC];
        end
        
    end
end

fprintf('We retained %d ICs of %d . \n', length(comp_ind), nICs);

end

