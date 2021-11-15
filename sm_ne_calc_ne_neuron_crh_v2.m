function exp_site_nedata = sm_ne_calc_ne_neuron_crh_v2(exp_site_nedata, mtffile)

% Inputs:
% exp_site_nedata:varaible saved in ensemble .mat files
% stimfile: string of stimfile where stim_mat is stored
if isfield(exp_site_nedata.nedata, 'NEcrh')
    return
end
nedata = exp_site_nedata.nedata;
fprintf('Calculating NEtrain CRH\n')
% 3 NE trians
NEtrain = nedata.NEtrain.all;
NEtrain_posi = nedata.NEtrain.posi;
NEtrain_neg = nedata.NEtrain.neg;
NEtrain_2018 = nedata.sta_NEtrain_2018;
nNE = size(NEtrain,1);
% calculate CRH
alltrain = [NEtrain; NEtrain_posi; NEtrain_neg; NEtrain_2018];
[mtfhist, tmf, smf] = sm_calculate_CRH(alltrain, mtffile);
% save results
exp_site_nedata.nedata.NEcrh.all = mtfhist(1: nNE,:);
exp_site_nedata.nedata.NEcrh.posi = mtfhist(nNE + 1:2 * nNE,:);
exp_site_nedata.nedata.NEcrh.neg = mtfhist(2 * nNE + 1:3 * nNE,:);
exp_site_nedata.nedata.NEcrh_2018 = mtfhist(3 * nNE + 1:4 * nNE,:);

% calculate neuron crh
if ~isfield(exp_site_nedata.nedata, 'neuroncrh')
    if isfield(nedata, 'sta_spktrain')
        spktrain = nedata.sta_spktrain;
    else
        spktrain = nedata.spktrain;
    end
    lendiff = size(NEtrain,2) - size(spktrain,2);
    if lendiff ~= 0
        spktrain = spktrain(:, 1:end+lendiff);
    end
    [mtfhist, ~, ~] = sm_calculate_CRH(spktrain, mtffile);
    exp_site_nedata.nedata.neuroncrh = mtfhist;
end
exp_site_nedata.nedata.crh_tmf = tmf;
exp_site_nedata.nedata.crh_smf = smf;



