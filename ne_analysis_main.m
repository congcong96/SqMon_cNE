%% configure paths
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

data_path = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh';
save_path = '/data/congcong/SqMoPhys_Josh/cNE_analysis/all_bins';
stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
stimfile = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');

%%  ---- PART1:cNE calculation ----------------------------------
cd(data_path)
files = dir('*-dmr*curated-thresh.mat');
files = {files.name};
%badfiles = ne_batch_save_halfms_binned_spktrain(files, stimfile);
binsize = [150, 100, 80, 50, 30, 20, 10, 8, 5, 2]; 
ne_rn_data_processing2(files, 2*binsize, 1, save_path)

%% ---- PART2: determine member neurons of cNE--------------------------------
% IC weights v.s. correlation of cNE train
cd(save_path)
nefiles = dir('*20dft.mat');
thresh = [0 25 37 50 75]/100;
for ii = 1:length(nefiles)
    fprintf('Processing %s\n', nefiles(ii).name)
    %sm_ne_corr_NEtrain_spktrain(nefiles(ii).name, 1)
    sm_ne_NEmember_thresh(nefiles(ii).name, thresh)
end

thresh = 0.37; % use 37% as threshold
for ii = 1:length(nefiles)
    fprintf('Processing %s\n', nefiles(ii).name)
    sm_ne_NEmember_thresh_set(nefiles(ii).name, thresh)
end

% to check IC weigth distribution and plot relationship of ICweights and
% correlation matrix, run 'sm_ne_plot_ICweight_distribution'

%% ---- PART2.2: get NEtrain ---------------
cd('/data/congcong/SqMoPhys_Josh/cNE_analysis')
nefiles = {nefiles.name};
sm_ne_batch_save_NEtrain(nefiles)
ne_batch_save_sta_NEtrain(nefiles)

%% ---- PART3: calculate cNE/neuron strf, rtf and crh--------------------------------
cd('/data/congcong/SqMoPhys_Josh/cNE_analysis')
nefiles = dir('*20dft.mat');
if ~exist('stim_mat', 'var')
    load(stimfile, 'stim_mat')
end
mtfstimfile = regexprep(stimfile, '_stim', '_mtf');
stimparamfile = regexprep(stimfile, '_stim', '_param');
load(stimparamfile, 'taxis', 'faxis')
for ii = 1:length(nefiles)
    fprintf('\n\nProcessing %s\n', nefiles(ii).name)
    load(nefiles(ii).name, 'exp_site_nedata')
    
    if exp_site_nedata.df < 10
        taxis_tmp = taxis(1:4:end);
    else
        taxis_tmp = taxis(1:10:end);
    end
    taxis_tmp = taxis_tmp(1:exp_site_nedata.nedata.nlags);
    
    % calculate strf
    exp_site_nedata = ne_calc_NE_strf_v2(exp_site_nedata, stim_mat, taxis_tmp, faxis); 
    % calculate rtf
    exp_site_nedata = sm_ne_calc_ne_neuron_rtf_v2(exp_site_nedata);
    % calculate crh
    exp_site_nedata = sm_ne_calc_ne_neuron_crh_v2(exp_site_nedata, mtfstimfile);
    
    % save result
    fprintf('Saving data in exp_site_nedata\n')
    save(nefiles(ii).name, 'exp_site_nedata', '-append');
    clear exp_site_nedata
end

%% ---- Scripts for plots--------------------------------
% stem plot of ICweights + CCG + strf and crh (2018, posi, neg)
sm_ne_plot_2018_posi_neg_ccg
% calculate and plot moron's I, unit location, and cell type
sm_ne_plot_MoronI
% get cell types and detemine it's relationship to negative members
sm_ne_plot_CellType
% get how many negative members are there in each NE
sm_plot_NE_type_distribution
% get CCG of individual neurons
sm_plot_posi_neg_CCG


