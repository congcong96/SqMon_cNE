function ne_rn_data_processing2(files, DF, roundopt, savefolder)

% files: cell array of names of files to be processed (files with spk)
% DF: vector of DFs (for time) to be processed
% roundopt: if 1, floors downsampled vector to remove the last incomplete
% bin.

% to use this code, binned spiketrains must first be calculated. see 
% ne_batch_save_halfms_binned_spktrain.m

if ~exist('roundopt','var')
    roundopt = 0;
end

% load stimulus
stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
uniquestim = {'contra-10'};
if isempty(uniquestim{1})
    uniquestim = {'rn1-15'};
end
if length(uniquestim) == 1
    parts = regexp(uniquestim{1}, '-', 'split');
    
    stimfilename = sprintf('%s-*%smin_DFt1_DFf5_stim.mat',parts{1}, parts{2});

    stimfile = fullfile(stimfolder,stimfilename);
    stimfile = gfn(stimfile,1);
    load(stimfile{1});
end    


for ii = 1:length(files)
    load(files{ii});

    if exist('evoked_spktrain','var')
        spktrain = evoked_spktrain;
        edges = evoked_edges;
        clear('evoked_spktrain', 'evoked_edges');
    end

    if ~exist('spktrain','var')
        error('Run ''ne_batch_save_halfms_binned_spktrain.m'' first!')
    end
    
    if isempty(spktrain)
        warning('There are no units in %s. Skipping...', files{ii})
        continue
    end
    
    if isfield(spk.spk, 'spiketimes_dmr')
        spktimes = {spk.spk.spiketimes_dmr};
        [spk.spk.spiketimes] = spktimes{:};
    end
    
    
    
    %base1 = regexp(files{i},'\<site.*db\>','match','once');
    base2 = regexp(files{ii},'(?<=min-).*(?=.mat)','match','once');
    base = [spk.exp '-dmr-' spk.stimlength '-' base2];
        
    if length(uniquestim) ~= 1
        
        stim = regexp(base,'rn\d{1,2}','match','once');
        stimlength = regexp(base, '\d{1,3}(?=(min))', 'match', 'once');

        stimfilename = sprintf('%s-*%smin_DFt1_DFf5_matrix.mat',stim, stimlength);
        stimfile = fullfile(stimfolder,stimfilename);
        stimfile = gfn(stimfile,1);
        load(stimfile{1});
        
    end
    
    position = {spk.spk.position};
    
    for j = 1:length(DF)
        
        DFt = DF(j);
        outfile = sprintf('%s-ne-%ddft.mat',base,DFt);
        outfile = fullfile(savefolder,outfile);
        
        outfile = regexprep(outfile, 'sponrn1', 'rn1'); % temp
        
        if exist(outfile,'file')
            continue
        else
            %% calculate NE
            nedata = ne_calc_neuronal_ensemble_sta(spktrain, stim_mat, 96000, DFt, edges, position, roundopt);
            if isempty(nedata)
                continue
            end
            save(outfile,'nedata')
            exp_site_nedata = ne_create_exp_site_nedata_file(outfile);
            
            %% get confidence interval for IC weights
%             CI = ne_calc_ICA_threshold(exp_site_nedata,'circular', 100, 'stdev', 1.5); %threshold currently at 1.5 stdev
%             exp_site_nedata.nedata.CI = CI;
            
            %% get membr neuron ids
%             NEmembers = ne_identify_NEmembers(exp_site_nedata.nedata.Patterns, exp_site_nedata.nedata.CI);
%             exp_site_nedata.nedata.NEmembers_2018 = NEmembers;
%             
%             %
%             if isfield(exp_site_nedata.nedata, 'NEthresh')
%                 fprintf('\ncNE activity thresh for %s already calculated\n', outfile);
%             else
%                 [thresh, alpha] = ne_calc_NE_act_thresholds(exp_site_nedata,'circular', 50, 99:0.1:99.9);
%                 exp_site_nedata.nedata.NEthresh_2018 = thresh;
%                 exp_site_nedata.nedata.NEthresh_alpha_2018 = alpha;
%             end
            
            close all
            save(outfile,'exp_site_nedata');      
        end
    end
end
%%
% cd('/data/congcong/SqMoPhys_Josh/cNE_analysis/dmr_std4_thresh')
% nefiles = gfn('*dmr*dft.mat');
% 
% ne_batch_calc_num_NE_events(nefiles); % number of elements in nedata.Activities cross NEthresh

%%
% cd('/data/congcong/SqMoPhys_Josh/cNE_analysis/dmr_std4_thresh')

%%
% stim = fullfile('/data/congcong/SqMoPhys_Josh/stim', 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');
% ne_batch_calc_NE_stamat_with_sta_NEtrain(nefiles, stim);

%%
%ne_batch_calc_ne_neuron_rtf(nefiles, 0);
% 
% %%
% ne_batch_calc_STA_statistics(nefiles);
% % 
% % %%
% ne_batch_calc_significant_sta(nefiles);
% % 
% % %%
% ne_batch_classify_significant_STAs(nefiles);


end


