function sm_ne_NEmember_thresh(nefile, thresh, flag_plot)

if nargin == 2
    flag_plot = 0;
end
load(nefile, 'exp_site_nedata')
nedata = exp_site_nedata.nedata;
Patterns = nedata.Patterns;
nNE = size(Patterns, 2);
NEmembers = struct([]);
for ii = 1:length(thresh)
    prc = thresh(ii);
    NEmembers(ii).thresh = prc;
    member = cell(1,nNE);
    for jj = 1:nNE
        pattern = Patterns(:,jj);
        m = max(abs(pattern));
        member{jj} = find(abs(pattern) > m * prc);
    end
    NEmembers(ii).member = member;
end
exp_site_nedata.nedata.NEmember = NEmembers;
save(nefile, 'exp_site_nedata', '-append');

% plot the membership of each threshold
if flag_plot
    savefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/NEmember_thresh/max_prc';
    for ii = 1:size(Patterns, 2)
        savename = sprintf('%s-cNE%d-%dms.jpg', exp_site_nedata.exp, ii, nedata.df/2);
        
        figure
        figuresetup2savepdf(15, 20)
        [pattern, idx] = sort(Patterns(:,ii), 'descend' );
        nplot = length(thresh);
        for jj = 1:nplot
            subplot(nplot,1,jj)
            stem(pattern, 'k','linewidth', 2)
            hold on
            prc = NEmembers(jj).thresh;
            title(sprintf('Threshold at %.0f%% of the maximum', prc*100))
            member_tmp = NEmembers(jj).member{ii};
            if ~isempty(member_tmp)
                member = zeros(size(member_tmp));
                for kk = 1:length(member_tmp)
                    member(kk) = find(idx == member_tmp(kk));
                end
                idx_posi = find(pattern(member) > 0);
                idx_neg = find(pattern(member) < 0);
            end
            if ~isempty(idx_posi)
                stem(member(idx_posi), pattern(member(idx_posi)), 'r', 'linewidth', 2)
            end
            if ~isempty(idx_neg)
                stem(member(idx_neg), pattern(member(idx_neg)), 'b', 'linewidth', 2)
            end
            thr = prc * max(abs(pattern));
            xlim([0 length(pattern)+1])
            plot(xlim, thr*[1 1], 'k--')
            plot(xlim, -thr*[1 1], 'k--')
        end
        saveas(gcf, fullfile(savefolder, savename))
    end
end


