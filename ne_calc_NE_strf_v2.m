function exp_site_nedata = ne_calc_NE_strf_v2(exp_site_nedata, stim_mat, taxis, faxis)
% calculating strf for NEtrain and NEtrain_2018
if isfield(exp_site_nedata.nedata,'NEstrf')
    return
end

fprintf('Calculating NEtrain STRF\n')
nedata = exp_site_nedata.nedata;
% stack all NEtrains together for quick_calc_sta
NEtrain = [nedata.NEtrain.all;nedata.NEtrain.posi; nedata.NEtrain.neg; nedata.sta_NEtrain_2018];
nNE = size(nedata.Patterns,2);

nlags = nedata.nlags;
% downsample stimulus
if nedata.df < 10
    stimulus = stim_mat(:, 4:4:end);
else
    stimulus = stim_mat(:, 10:10:end);
    if nedata.df > 10
        diff = size(stimulus, 2) - size(NEtrain, 2);
        if diff> 0
            stimulus(:,end-diff+1:end) = [];
        end
    end
end

NE_stamat = quick_calc_sta(stimulus, NEtrain, 'nlags', nlags, 'chunks', 1, 'suppressprint', 1);
exp_site_nedata.nedata.NEstrf.all = NE_stamat(1:nNE, :);
exp_site_nedata.nedata.NEstrf.posi = NE_stamat(nNE+1:nNE*2, :);
exp_site_nedata.nedata.NEstrf.neg = NE_stamat(nNE*2+1:nNE*3, :);
exp_site_nedata.nedata.NEstrf_2018 = NE_stamat(nNE*3+1:nNE*4, :);


exp_site_nedata.nedata.strf_faxis = faxis;
exp_site_nedata.nedata.strf_taxis = taxis(1:nlags);
