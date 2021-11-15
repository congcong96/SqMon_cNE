function NEmembers = ne_identify_NEmembers(patterns, CI)

% Identifies assembly members based on thresholds calculated in
% ca_calc_ICA_threshold.m.
%   
%   cadata: must include cadata.CI

%   Updated 7/20/16 by js
%   Updated 10/9/17 by js. Inputs are now ICweights (patterns) and CI
%   (from ne_calc_ICA_threshold.m)

NEmembers = cell(size(patterns,2),1);

for i = 1:size(patterns,2)
    
    idx1 = find(patterns(:,i) <= CI(1));
    idx2 = find(patterns(:,i) >= CI(2));
    NEmembers{i} = sort([idx1; idx2]);

end

