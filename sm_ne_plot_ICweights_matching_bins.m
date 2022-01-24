save_path = '/data/congcong/SqMoPhys_Josh/cNE_analysis/all_bins';
binsize = [150, 100, 80, 50, 30, 20, 10, 8, 5, 2]; % in ms
reference = 10;

%% get matched ICs ----------------------------------------------------

% pearson's correlation
figurefolder = sprintf('/data/congcong/SqMoPhys_Josh/figure/cNE/binsize_stability_%dms_pearson', reference);
if ~exist(figurefolder, 'dir')
    mkdir(figurefolder)
end
sm_ne_ICweights_matching_bins(save_path, 2*binsize, 2*reference, figurefolder, 'Pearson')
% spearman's correlation
figurefolder = sprintf('/data/congcong/SqMoPhys_Josh/figure/cNE/binsize_stability_%dms_spearman', reference);
if ~exist(figurefolder, 'dir')
    mkdir(figurefolder)
end
sm_ne_ICweights_matching_bins(save_path, 2*binsize, 2*reference, figurefolder, 'Spearman')

%% distribution of significantly matched ICs----------------------------
reference = 150;
cd(save_path)
df = binsize * 2;
% get total number of cNEs
nNE = zeros(1, length(binsize));
for ii = 1:length(df)
    files = dir(sprintf('*-%ddft.mat', df(ii)));
    for jj = 1:length(files)
        load(files(jj).name,'exp_site_nedata')
        nNE(ii) = nNE(ii) + size(exp_site_nedata.nedata.Patterns, 2);
    end
end

% get number of matched NEs
matchedNE = zeros(1, length(binsize));
ref_idx = find(df == reference * 2);
files = dir(sprintf('*-%ddft.mat', reference * 2));
for ii = 1:length(files)
    load(files(ii).name, 'IC_matched_pearson')
    IC_matched = IC_matched_pearson;
    for jj = 1:length(IC_matched)
        p = IC_matched(jj).p;
        df_match = IC_matched(jj).df;
        df_match_ref_idx = find(df_match == reference*2);
        idx = find(p(df_match_ref_idx :end) >= 0.05);
        if ~isempty(idx)
            sig_df_idx = ref_idx : (ref_idx+idx);
        else
            sig_df_idx = ref_idx : (ref_idx+length(p)-df_match_ref_idx+1);
        end
        matchedNE(sig_df_idx) = matchedNE(sig_df_idx) + 1;
        
        idx = find(flip(p(1:df_match_ref_idx-1)) >= 0.05);
        if ~isempty(idx)
            sig_df_idx = ref_idx-idx+1:ref_idx-1;
        else
            sig_df_idx = ref_idx-df_match_ref_idx+1:ref_idx-1;
        end
        matchedNE(sig_df_idx) = matchedNE(sig_df_idx) + 1;

        
    end
end
%
h = bar([matchedNE', (nNE-matchedNE)'], 'stacked');
h(1).FaceColor = 0.4*[1 1 1];
h(2).FaceColor = 0.9*[1 1 1];
ylabel('Number of cNEs')
xticks(1:10)
xticklabels(df/2)
xlabel('Binsize (ms)')
title("Matching significance by pearson's correlation")


%% distribution of significantly matched ICs (core and belt)----------------------
probrecord = '/data/congcong/SqMoPhys_Josh/ProbeRecord.mat';
load(probrecord, 'ProbeRecord')
reference = 150;
cd(save_path)
df = binsize * 2;
% get total number of cNEs
nNE = zeros(1, length(binsize));
for ii = 1:length(df)
    files = dir(sprintf('*-%ddft.mat', df(ii)));
    for jj = 1:length(files)
        load(files(jj).name,'exp_site_nedata')
        nNE(ii) = nNE(ii) + size(exp_site_nedata.nedata.Patterns, 2);
    end
end

% get number of matched NEs
matchedNE = zeros(1, length(binsize));
ref_idx = find(df == reference * 2);
files = dir(sprintf('*-%ddft.mat', reference * 2));
for ii = 1:length(files)
    load(files(ii).name, 'IC_matched_pearson')
    IC_matched = IC_matched_pearson;
    for jj = 1:length(IC_matched)
        p = IC_matched(jj).p;
        df_match = IC_matched(jj).df;
        df_match_ref_idx = find(df_match == reference*2);
        idx = find(p(df_match_ref_idx :end) >= 0.05);
        if ~isempty(idx)
            sig_df_idx = ref_idx : (ref_idx+idx);
        else
            sig_df_idx = ref_idx : (ref_idx+length(p)-df_match_ref_idx+1);
        end
        matchedNE(sig_df_idx) = matchedNE(sig_df_idx) + 1;
        
        idx = find(flip(p(1:df_match_ref_idx-1)) >= 0.05);
        if ~isempty(idx)
            sig_df_idx = ref_idx-idx+1:ref_idx-1;
        else
            sig_df_idx = ref_idx-df_match_ref_idx+1:ref_idx-1;
        end
        matchedNE(sig_df_idx) = matchedNE(sig_df_idx) + 1;

        
    end
end
%
h = bar([matchedNE', (nNE-matchedNE)'], 'stacked');
h(1).FaceColor = 0.4*[1 1 1];
h(2).FaceColor = 0.9*[1 1 1];
ylabel('Number of cNEs')
xticks(1:10)
xticklabels(df/2)
xlabel('Binsize (ms)')
title("Matching significance by pearson's correlation")