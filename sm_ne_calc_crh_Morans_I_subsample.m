function I = sm_ne_calc_crh_Morans_I_subsample(spktrain, mtffile, weightsmat, nevent, nrepeat)

% INPUTS:
%   spktrain: single spike train to resample and calculate Moran's I from
%   mtffile: string of the path of stimulus file to calculate crh
%   weightsmat: 1xn cell contains weightsmat to calculate Moran's I
%   nevent: number of events in subsampled spike trains
%   nrepeat: number of repeat to subsample

% get the the time bin index for every spike
spk_idx = find(spktrain > 0);
for ii = 1:10
    idx_tmp = find(spktrain > ii);
    if isempty(idx_tmp)
        break
    else
        spk_idx = [spk_idx, idx_tmp];
    end
end

% get subsampled spike trains
spktrain_subsample = zeros(nrepeat, size(spktrain,2));
for ii = 1:nrepeat
    idx_subsample = spk_idx(randsample(length(spk_idx), nevent));
    [c, idx_subsample] = histcounts(idx_subsample, [unique(idx_subsample), Inf]);
    spktrain_subsample(ii, idx_subsample(1:end-1)) = c;
end

%calculate crh for resampled spike trains
crh = sm_calculate_CRH(spktrain_subsample, mtffile);

% calculate Moran' I for resampled crh
I = zeros(length(weightsmat), nrepeat);
for ii = 1:nrepeat
    for jj = 1:length(weightsmat)
        I(jj, ii) = morans_i(crh(ii,:), weightsmat{jj});
    end
end
