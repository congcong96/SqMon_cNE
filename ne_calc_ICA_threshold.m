function [CI] = ne_calc_ICA_threshold(exp_site_nedata,varargin)

% Get thresholds for cell assembly data to determine which units are
% 'considered part of an assembly'. Makes use of a random permutation test 
% to get the full range of ICA values by overriding the Marcenko-Pastur
% test and taking the first N PCs, where N = number of assemblies obtained
% from running the Marcenko-Pastur on the original dataset.
%
% Description of inputs:
%
%   spktrain:   Matrix of spike counts per time bin. Rows are single units;
%               columns are time bins.
%   perm:       If 'random', spike trains are shuffled randomly, preserving
%               only number of spikes in each spike train. If 'circular',
%               spike trains are circularly shifted, preserving the number
%               of spikes and all ISIs except for those at the edges.
%   num_iter:   Number of iterations for random permutation test. Default
%               set to 100.
%   prctile:    Confidence interval. For example, 95% returns 2.5 and 97.5
%               percentiles. Default set to [].
%   stdev:      Number of standard deviations from the mean. For example, 2
%               returns 2.5 and 97.5 percentile. Default set to 1.5. 
%   opts:       Struct array containing options for assembly_patterns.m.
%               Optional.
%
%   Default number of iterations for fast_ICA for original data is 500.
%   However, this is 100 for the random permutations for speed. This 
%   can be changed in opts. 
%
%   Updated 3/17/17 by js, added std for alpha
%   Updated 10/20/17 by js, moved the bulk of the code for shuffling to
%   ne_calc_shuffled_spktrain_ICs.m such that other codes can use the
%   shuffling procedure.

p = inputParser;
addOptional(p, 'perm', 'circular', @(x) strcmp(x, 'circular') || strcmp(x, 'random'))
addOptional(p, 'num_iter', 100, @(x) x >= 1)
addParameter(p, 'prctile', [], @(x) isempty(x) || (x >= 0 && x <= 100))
addParameter(p, 'stdev', 1.5, @(x) isempty(x) || (x >=0 && x < 10))
addParameter(p, 'opts', [], @(x) isempty(x) || isstruct(x))
% addParameter(p, 'multiple', [], @iscell)
parse(p, varargin{:})

perm = p.Results.perm;
num_iter = p.Results.num_iter;
prctile = p.Results.prctile;
stdev = p.Results.stdev;
opts = p.Results.opts;

if ~isempty(prctile) && ~isempty(stdev)
    error('Choose either prctile or stdev.')
elseif isempty(prctile) && isempty(stdev)
    error('Enter either prctile or stdev.')
end

if isfield(exp_site_nedata, 'spktrain')
    spktrain = exp_site_nedata.spktrain;
else
    spktrain = exp_site_nedata.nedata.spktrain;
end

ICshufflemat = ne_calc_shuffled_spktrain_ICs(spktrain,perm,num_iter, 'opts', opts);

%get lower and upper percentiles based on alpha
if ~isempty(prctile)
    lowerprc = (100 - prctile)/2;
    upperprc = 100 - lowerprc;
    CI = prctile(ICshufflemat(:),[lowerprc upperprc]);
elseif ~isempty(stdev)
    miu = mean(ICshufflemat(:));
    sigma = std(ICshufflemat(:));
    CI = [miu - sigma*stdev, miu + sigma*stdev];
end



    
