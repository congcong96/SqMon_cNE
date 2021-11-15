function ne_batch_save_sta_NEtrain(nefiles, NEthreshalpha, sponopt)

if ~exist('NEthreshalpha','var')
    NEthreshalpha = 99.5;
end

if ~exist('sponopt','var')
    sponopt = 0;
end

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    fprintf('Processing %d', nefiles{i})
    nedata = exp_site_nedata.nedata;
    NEact = nedata.Activities;
    
    if sponopt
        sponlen = nedata.sponlen;
        NEact = NEact(:, sponlen + 1:end);
    end
    
    NEmembers = nedata.NEmembers;
    spktrain = nedata.sta_spktrain;
    
    NEtidx = nedata.NEthresh_alpha_2018 == NEthreshalpha;
    
    NEthresh = nedata.NEthresh_2018(NEtidx,:);
    
    NEtrain = ne_upsample_NEact_using_member_neuron_activity(NEact, NEmembers, spktrain, NEthresh);
    exp_site_nedata.nedata.sta_NEtrain_2018 = NEtrain;
    
    save(nefiles{i}, 'exp_site_nedata', '-append');

end
    