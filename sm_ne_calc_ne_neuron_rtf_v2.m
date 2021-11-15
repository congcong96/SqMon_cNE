function exp_site_nedata = sm_ne_calc_ne_neuron_rtf_v2(exp_site_nedata)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
if isfield(exp_site_nedata.nedata, 'NErtf')
    return
end

nedata = exp_site_nedata.nedata;

fprintf('Calculating NEtrain RTF\n')

nf = nedata.nf;
nlags = nedata.nlags;
taxis = nedata.strf_taxis;
faxis = nedata.strf_faxis;
stimstr.paramfile = fullfile('/data/congcong/SqMoPhys_Josh/stim', 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt10_DFf5_param.mat');
load(stimstr.paramfile, 'MaxFM','MaxRD')

% ger rtf of individual neurons
if ~isfield(exp_site_nedata.nedata, 'neuronstrf')
    
    neuronrtf = strf2rtf(nedata.stamat, nf, nlags, taxis, faxis, MaxFM, MaxRD);
    exp_site_nedata.nedata.neuronrtf = neuronrtf;
end
% get rtf of NE sta
NErtf_all = strf2rtf(nedata.NEstrf.all, nf, nlags, taxis, faxis, MaxFM, MaxRD);
exp_site_nedata.nedata.NErtf.all = NErtf_all;

% get rtf of NE_posi sta
NErtf_posi = strf2rtf(nedata.NEstrf.posi, nf, nlags, taxis, faxis, MaxFM, MaxRD);
exp_site_nedata.nedata.NErtf.posi = NErtf_posi;

% get rtf of NE_neg sta
NErtf_neg = strf2rtf(nedata.NEstrf.neg, nf, nlags, taxis, faxis, MaxFM, MaxRD);
exp_site_nedata.nedata.NErtf.neg = NErtf_neg;

% get rtf of NE_2018
[NErtf_2018, tmf, xmf] = strf2rtf(nedata.NEstrf_2018, nf, nlags, taxis, faxis, MaxFM, MaxRD);
exp_site_nedata.nedata.NErtf_2018 = NErtf_2018;
   
exp_site_nedata.nedata.rtf_tmf = tmf;
exp_site_nedata.nedata.rtf_smf = xmf;
end

function [rtf, tmf, xmf] = strf2rtf(strf, nf, nlags, taxis, faxis, MaxFM, MaxRD)
tmodbins = 16;
smodbins = 20;
rtf = zeros(size(strf,1), (smodbins + 1) * (2 * tmodbins + 1));

for ii = 1:size(strf, 1)
    [tmf, xmf, temprtf] = sm_mtf_sta2rtf(reshape(strf(ii,:), nf, nlags), taxis, faxis, MaxFM, MaxRD, tmodbins, smodbins);
    rtf(ii,:) = temprtf(:)';
end
end
