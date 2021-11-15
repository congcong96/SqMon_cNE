function weightsmat = create_rf_spatial_autocorr_weights_matrix(nrows, ncols, diagopt)

% Create spatial autocorrelations weights matrix. Mainly for use with
% morans_i.m.

% Inputs:
%
%   nrows: Number of rows for receptive field.
%   ncols: Number of columns for receptive field
%   diagopt: Determines if adjacent diagonal pixels are considered adjacent.
%            Default is True.


if ~exist('diagopt','var')
    diagopt = 1;
end

% define searches
eq{1} = @(i,j) [i-1 j];
eq{2} = @(i,j) [i+1 j];
eq{3} = @(i,j) [i j-1];
eq{4} = @(i,j) [i j+1];

if diagopt %include diagonals as contiguous pixels
    eq{5} = @(i,j) [i-1 j-1];
    eq{6} = @(i,j) [i+1 j-1];
    eq{7} = @(i,j) [i-1 j+1];
    eq{8} = @(i,j) [i+1 j+1];
end


numelements = nrows * ncols;
weightsmat = zeros(numelements);
refmat = reshape(1:numelements, nrows, ncols);

n = 1;

for j = 1:ncols
    
    for i = 1:nrows

        % apply all equations to find 2d-indices of contiguous pixels
        rcidx = cellfun(@(x) x(i,j), eq, 'UniformOutput', 0);
        % remove all 2d-indices smaller than 1 and larger than nrows/ncols
        rcidx(cellfun(@(x) x(1) < 1 | x(1) > nrows | x(2) < 1 | x(2) > ncols, rcidx)) = [];
        % get 1d-indices from 2d-indices
        idx = cellfun(@(x) refmat(x(1), x(2)), rcidx);
                
        weightsmat(sub2ind(size(weightsmat), n*ones(length(idx),1), idx(:))) = 1;
        
        n = n+1;
    end
end
            
        
    
    
end

