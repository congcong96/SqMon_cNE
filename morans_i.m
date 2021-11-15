function I = morans_i(matrix, weightsmat) 

% Runs Moran's I. 

% Inputs:
%   matrix: Receptive field (or something else) on which Moran's I
%   needs to be calculated
%
%   weightsmat: Spatial autocorrelation weights mat. If none is provided,
%               then it is assumed that the matrix input is a receptive
%               field and the default spatial autocorrelation matrix
%               (where diagonals are considered adjacent) is generated.



% Calculate weights matrix if not given as input
if ~exist('weightsmat','var')
    [nrows, ncols] = size(matrix);
    weightsmat = create_rf_spatial_autocorr_weights_matrix(nrows, ncols, 1);
end

% Calculate N and W
N = numel(matrix);
W = sum(weightsmat(:));

% Calculate denominator
xbar = mean(matrix(:));
dvec = matrix(:) - xbar;
deno = sum(dvec.^2);

% Calculate numerator
nummat = zeros(N);
[i,j] = find(weightsmat);

for k = 1:length(i)
    
    nummat(i(k),j(k)) = dvec(i(k)) * dvec(j(k)); 

end

nume = sum(nummat(:));

I = (N/W) * (nume/deno); 