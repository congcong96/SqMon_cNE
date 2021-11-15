function ICshufflemat = ne_calc_shuffled_spktrain_ICs(spktrain,varargin)

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
%   stdev:      Number of standard deviations from the mean. For example,
%               1.5 returns 2.5 and 97.5 percentile. Default set to 1.5. 
%   opts:       Struct array containing options for assembly_patterns.m.
%               Optional.
%
%   Default number of iterations for fast_ICA for original data is 500.
%   However, this is 100 for the random permutations for speed. This 
%   can be changed in opts. 
%
%   Updated 3/17/17 by js, added std for alpha

p = inputParser;
addOptional(p, 'perm', 'circular', @(x) strcmp(x, 'circular') || strcmp(x, 'random'))
addOptional(p, 'num_iter', 100, @(x) x >= 1)
addParameter(p, 'opts', [], @(x) isempty(x) || isstruct(x))
parse(p, varargin{:})

perm = p.Results.perm;
num_iter = p.Results.num_iter;
opts = p.Results.opts;


     
% if isstruct(spktrain)
%     if isfield(spktrain, 'spktrain')
%         spktrain = spktrain.spktrain;
%     else
%         spktrain = spktrain.nedata.spktrain;
%     end
% else
%     spktrain = spktrain;
% end

%process original data
if isempty(opts)
    patterns = assembly_patterns(spktrain);
else
    patterns = assembly_patterns(spktrain,opts);
end

num_assemblies = size(patterns,2);

%for all the ICA values from permutation.
ICshufflemat = zeros(size(patterns,1),size(patterns,2),num_iter);

for i = 1:num_iter
    
    fprintf('%d of %d iterations...\n',i,num_iter);
    spktrainperm = zeros(size(spktrain));
    [M, N] = size(spktrainperm);
    
    switch perm
        case 'random'
   
            for j = 1:M
                X = rand(N,1);
                [~,idx] = sort(X);
                spktrainperm(j,:) = spktrain(j,idx);      
            end

        case 'circular'
        
            shiftidx = random('unid',N,[M 1]);
    %         shiftidx = cell2mat(arrayfun(@(x) [x:N,1:x-1],shiftidx,'UniformOutput',0));            

            for j = 1:M
                spktrainperm(j,:) = circshift(spktrain(j,:),[0 shiftidx(j)]);
            end
    end
            
        
    %run each iteration through the modified assembly_patterns code
    if isempty(opts)
        perm_patterns = assembly_patterns_permutation(spktrainperm, ...
            num_assemblies);
    else
         perm_patterns = assembly_patterns_permutation(spktrainperm, ...
            num_assemblies,opts);
    end
    
    ICshufflemat(:,:,i) = perm_patterns;
    
end


    
