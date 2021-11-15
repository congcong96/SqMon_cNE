function sm_ne_batch_save_NEtrain(nefiles)


for ii = 1:length(nefiles)
    load(nefiles{ii}, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    fprintf('Processing %s\n', nefiles{ii})
    if isfield(nedata, 'sta_spktrain')
        spktrain = nedata.sta_spktrain;
    else
        spktrain = nedata.spktrain;
    end
    
    NEmembers = nedata.NEmembers;
    NEmembers_posi = cell(size(NEmembers));
    NEmembers_neg = cell(size(NEmembers));
    Patterns = nedata.Patterns;
    nNE = length(NEmembers);
    NEact_posi = zeros(nNE, size(nedata.spktrain,2));
    NEact_neg = zeros(nNE, size(nedata.spktrain,2));
    spkNEtrain = nedata.spktrain;
    spkNEtrain(spkNEtrain>0) = 1;
    for jj = 1:nNE
        MemberWeights = Patterns(NEmembers{jj},jj);
        idxPosi = NEmembers{jj}(MemberWeights > 0);
        idxNeg = NEmembers{jj}(MemberWeights < 0);
        NEmembers_posi{jj} = idxPosi;
        NEmembers_neg{jj} = idxNeg;
        if length(idxPosi) > 1
            act_tmp = sum(spkNEtrain(idxPosi,:));
            act_tmp(act_tmp <= 1) = 0;
%             exclude_tmp = sum(spkNEtrain(idxNeg,:), 1);
%             act_tmp(exclude_tmp > 0) = 0;
            NEact_posi(jj,:) = act_tmp;
        end
        if length(idxNeg) > 1
            act_tmp = sum(spkNEtrain(idxNeg,:));
            act_tmp(act_tmp <= 1) = 0;
%             exclude_tmp = sum(spkNEtrain(idxPosi,:),1);
%             act_tmp(exclude_tmp > 0) = 0;
            NEact_neg(jj,:) = act_tmp;
        end
    end
    actThresh = zeros(1, size(NEact_posi,1));
    NEtrain.posi = ne_upsample_NEact_using_member_neuron_activity(NEact_posi, NEmembers_posi, spktrain, actThresh);
    NEtrain.neg = ne_upsample_NEact_using_member_neuron_activity(NEact_neg, NEmembers_neg, spktrain, actThresh);
    NEtrain.all = NEtrain.posi + NEtrain.neg;
    exp_site_nedata.nedata.NEtrain = NEtrain;
    save(nefiles{ii}, 'exp_site_nedata', '-append');
end
