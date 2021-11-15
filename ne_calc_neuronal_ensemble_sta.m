function [nedata] = ne_calc_neuronal_ensemble_sta(spktrain, stim_mat, FsDVD, DF, edges, position, roundopt, plotopt)
% ca_calc_cell_assembly_sta Estimate STAs from PCA/ICA cell assembly analysis
% 
%     [cadata] = ca_calc_cell_assembly_sta(spk, stimulus, trigger, FsDVD, totalDF)
%     --------------------------------------------------------------------
%     spk : struct array of spike times. Obtaines from saveSpikeSort.m
%
%     stimulus : ripple stimulus in matrix format
%
%     trigger : trigger for ripple stimulus
%
%     FsDVD : sampling of DVD Audio system
%
%     DF : downsampling factor (1 is 0.5 ms, 2 is 1 ms, etc.)
%
%     cadata : struct holding calculations. Has the form:
%         cadata.spktrain = binned spike train, one row per neuron
%         cadata.fsdvd = sampling rate of sound stimulus
%         cadata.df = total downsampling factor of ripple noise envelope
%         cadata.position = recording depth. position(i) -> spktrain(i,:)
%         cadata.Patterns = cell assemblies
%         cadata.Activities = time course of cell assembly activity; one row per assembly
%         cadata.nf = number of frequencies in STA
%         cadata.nlags = number of time lags in STA
%         cadata.stamat = matrix of spike train STAs.
%         cadata.ca_stamat = matrix of cell assembly STAs.
%
%   updated 7/20/16 JS

if ~exist('roundopt','var')
    roundopt = 0;
end

if ~exist('plotopt','var')
    plotopt = 0;
end

spkmat = downsample_spiketrain(spktrain, DF);
stimulus = stim_mat(:,1:DF:end);

if roundopt == 1 && mod(size(spktrain,2), DF) ~= 0
    spkmat = spkmat(:,1:end-1);
    stimulus = stimulus(:,1:end-1);
end

edges = edges(1:DF:end);

if roundopt == 0
    maxdiff = max(unique(diff(edges)));
    edges = [edges edges(end)+maxdiff];
end


% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spkmat);
if isempty(Patterns)
    nedata = [];
    return
end

if plotopt
    ca_plot_cell_assemblies_data(spkmat, Patterns, Activities, position);
end

% Get stimulus-response matrices for later STA calculations
fprintf('\nForming Stimulus observation matrix\n');
if DF <= 5 
    nlags = 40;
else
    nlags = 20;
end

if DF > 10
    sta_spktrain = downsample_spiketrain(spktrain, 10);
    sta_stimulus = stim_mat(:,1:10:end);
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(sta_stimulus, sta_spktrain, nlags);
else
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(stimulus, spkmat, nlags);
end


% Calculate STAs in one block for each neuron.
%   stim = #trials X #dims, #dims = nfreq * nlags
%   resp = # neurons X #trials
%   resp * stim = #neurons X #dims
stamat = resp * stim; % one sta per row
[nf, ~] = size(stimulus);
nlags = size(stamat,2) / nf; % Get # time lags from sta and stimulus


% Calculate STAs in one block for each cell assembly.
%   stim = #trials X #dims, #dims = nfreq * nlags
%   Activities = # cell assemblies X #trials
%   Activities * stim = #neurons X #dims
if DF > 10
    mult = round(DF / 10);
    sta_act = zeros(size(Activities,1), size(Activities,2) * mult);
    for i = 1:size(Activities,1)
        sta_act(i,:) = interp(Activities(i,:), mult);
    end
    if size(sta_act,2) == size(sta_spktrain, 2)
        ne_stamat = sta_act(:, nlags:end) * stim;
    elseif size(sta_act,2) > size(sta_spktrain, 2)
        difference = size(sta_act,2) - size(sta_spktrain,2);
        ne_stamat = sta_act(:, nlags:end-difference) * stim;
    else
        difference = size(sta_spktrain,2) - size(sta_act,2);
        ne_stamat = sta_act(:, nlags:end) * stim(1:end-difference,:);
    end
else
    ne_stamat = Activities(:,nlags:end) * stim;
end

if plotopt
    % Plot the cell assemblies:
    ca_plot_cell_assembly_stamat(ne_stamat, nf, nlags, 'CA');

    % Plot the single units:
    ca_plot_cell_assembly_stamat(stamat, nf, nlags, 'STA');    
end

% Assign data to struct for output argument
nedata.spktrain = spkmat;
nedata.fsdvd = FsDVD;
nedata.df = DF;
nedata.position = position;
nedata.Patterns = Patterns;
nedata.Activities = Activities;
nedata.nf = nf;
nedata.nlags = nlags;
nedata.stamat = stamat;
nedata.ca_stamat = ne_stamat;
nedata.edges = edges;

if DF > 10
    nedata.sta_spktrain = sta_spktrain;
else
    nedata.sta_spktrain = spkmat;
end



return;




% Code for significance testing, if we want it.
% nreps = 4;
% pval = 0.05;
% for i = 1:size(stamat,1)
%     [sta_sig, siglevel, rand_dist] = ...
%         ca_sig_sta_from_stim_obs_resp(stamat(i,:), resp(i,:), stim, nreps, pval);
%     figure;
%     subplot(1,2,1);
%     imagesc(reshape(stamat(i,:),nf,nlags));
%     subplot(1,2,2);
%     imagesc(reshape(sta_sig,nf,nlags));
% end % (for i)




















