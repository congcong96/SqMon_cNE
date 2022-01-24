function sm_ne_ICweights_matching_bins(datafolder, df, reference, figurefolder, corrType)


cd(datafolder)
if nargin == 4
    corrType = 'Pearson';
end
df = sort(df, 'descend');
ref_idx = find(df == reference);
% split timebins into two parts: 
% df1 the part larger than reference timebin
% df2 the part shorter than the reference
df1 = flip(df(1:ref_idx));
df2 = df(ref_idx:end);
files = dir(sprintf('*-%ddft.mat', reference));
for ii = 1:length(files)    
    load(files(ii).name, 'exp_site_nedata')
    % get mached IC and corresponding correlations and p-values
    % for the two parts of time bins split at the reference
    [Patterns1, corrVal1, pVal1, df_match1, nNE]= match_ICweights(files(ii).name, df1, corrType);
    [Patterns2, corrVal2, pVal2, df_match2, ~]= match_ICweights(files(ii).name, df2, corrType);
    
    % combine results from the two parts
    Patterns = [flip(Patterns1), Patterns2(2:end)];
    Patterns = cell2mat(Patterns);
    nNeuron =  size(Patterns,1);
    Patterns = reshape(Patterns, nNeuron, nNE, size(Patterns,2)/nNE);
    
    corrVal = [fliplr(corrVal1) corrVal2];
    pVal = [fliplr(pVal1) pVal2];
    
    IC_matched = struct('ICs', [], 'corr', [], 'p', [], 'df', []);
    for jj = 1:nNE
        ICs = squeeze(Patterns(:, jj, :));
        c = corrVal(jj,:);
        p = pVal(jj,:);
        df_match = [flip(df_match1{jj}) reference df_match2{jj}];
        ICs(:,ICs(1,:)==0) = [] ;
        c(c == 0) = [];
        p(p == 0) = [];
       
        IC_matched(jj).ICs = ICs;
        IC_matched(jj).corr = c;
        IC_matched(jj).p = p;
        IC_matched(jj).df = df_match;
    end
    if strcmp(corrType, 'Pearson')
        IC_matched_pearson = IC_matched;
        save(files(ii).name, 'IC_matched_pearson', '-append')
    elseif strcmp(corrType, 'Spearman')
        IC_matched_spearman = IC_matched;
        save(files(ii).name, 'IC_matched_spearman', '-append')
    end
    % draw ICs and correlation matrix
    for jj = 1:nNE
        figure
        figuresetup2savepdf(30, 20)
        patterns = IC_matched(jj).ICs;
        thresh = 1/sqrt(size(patterns,1));
        df_match = IC_matched(jj).df;
        idx_start = find(df == df_match(1));
        idx_end = find(df == df_match(end));
        for kk = 1:size(patterns,2)
            subplot(3,10, idx_start - 1 + [kk, 10+kk])
            % stemplot
            stem(patterns(:,kk), 'k', 'filled')
            member = find(abs(patterns(:,kk)) > thresh);
            hold on
            stem(member, patterns(member,kk), 'r', 'filled')
            % fill the subthreshold region
            fill([1, nNeuron, nNeuron, 1], [-thresh, -thresh, thresh, thresh], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            
            view(90, 90)
            box off
            set(gca, 'XDir','reverse')
            xlim([1, nNeuron])
            ylabel(sprintf('%dms', df(kk)/2))
            if kk > 1
                set(gca, 'xcolor','none')
            else
                xlabel('Neuron #')
            end
        end
        
        subplot(3,10, 21:30)
        x = idx_start:idx_end-1;
        plot(x, IC_matched(jj).corr, '-ok', 'MarkerFaceColor','w')
        hold on
        p = IC_matched(jj).p;
        sigidx = find(p < 0.05);
        plot(x(sigidx), IC_matched(jj).corr(sigidx), 'ok','MarkerFaceColor','k')
        ylim([0 1])
        xlim([0 10])
        box off
        ylabel(sprintf("%s's coefficient", corrType))
        xticks(1:9)
        xticklabels({'150-100', '100-80', '80-50', '50-30', '30-20', '20-10', '10-8', '8-5', '5-2'})
        xlabel('compared binwidth (ms)')
        
        sgtitle(sprintf('%s  cNE#%d',exp_site_nedata.exp, jj), 'Interpreter', 'none')
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('%s-matchedIC-cNE%d', exp_site_nedata.exp, jj)));
        close
    end
end
end

function Patterns = pattern_direction(Patterns)
for ii = 1:size(Patterns, 2)
    [~, idx] = max(abs(Patterns(:,ii)));
    if Patterns(idx,ii) < 0
        Patterns(:,ii) = -Patterns(:,ii);
    end
end
end

function [Patterns, corrVal, pVal, df_match, nNE]= match_ICweights(filename, df, corrType)

fileID = filename(1:13);
load(filename, 'exp_site_nedata')
Pattern1 = exp_site_nedata.nedata.Patterns;
Pattern1 = pattern_direction(Pattern1);
nNE = size(Pattern1, 2);
Patterns = cell(1, length(df));
Patterns{1} = Pattern1;
corrVal = zeros(nNE, length(df)-1);
pVal = zeros(nNE, length(df)-1);

for jj = 2:length(df)
    filename = dir(sprintf('%s-*-%ddft.mat', fileID, df(jj)));
    if ~isempty(filename)
        Patterns{jj} = zeros(size(Pattern1));
        load(filename.name, 'exp_site_nedata')
        Pattern2 = exp_site_nedata.nedata.Patterns;
        Pattern2 = pattern_direction(Pattern2);
        [CorrMat, p] = corr(Pattern1, Pattern2, 'Type', corrType);
        CorrMat = abs(CorrMat);
        % match ICs
        while any(CorrMat(:) > 0)
            [M, I] = max(CorrMat, [], 'all', 'linear');
            [row, col] = ind2sub(size(CorrMat), I);
            CorrMat(row,:) = 0;
            CorrMat(:,col) = 0;
            corrVal(row, jj-1) = M;
            pVal(row, jj-1) = p(row, col);
            Patterns{jj}(:,row) = Pattern2(:,col);
        end
    else
        break
    end
    Pattern1 = Patterns{jj};
end

df_match = cell(1, size(corrVal,1));
for ii = 1:length(df_match)
    idx_stop = find(corrVal(ii,:) == 0);
    if idx_stop
        df_match{ii} = df(2:idx_stop);
    else
        df_match{ii} = df(2:end);
    end
end
end