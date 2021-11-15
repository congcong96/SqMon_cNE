%% set path
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

datapath = '/data/congcong/SqMoPhys_Josh/cNE_analysis';
spkfolder = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh';
cd(datapath)
nefiles = dir('*-20dft.mat');
%% plot CCG of posi and neg neurons
ccgfolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/ccg';
for ii = 3%1:length(nefiles)
    
    % load data
    load(nefiles(ii).name, 'exp_site_nedata')
    exp = exp_site_nedata.exp;
    nedata = exp_site_nedata.nedata;
    spkfile = dir(fullfile(spkfolder, [exp, '*']));
    load(fullfile(spkfile.folder, spkfile.name), 'spktrain')
    spktrain = spktrain(:,1:end-1);
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    nfigure = 1;
    for jj = 2%1:nNE
        member = NEmembers{jj};
        pattern = nedata.Patterns(:, jj);
        neg_idx = member(pattern(member) < 0);
        if isempty(neg_idx)
            continue
        end
        posi_idx = member(pattern(member) > 0);
        posi_train = spktrain(posi_idx, :)';
        for k1 = 1:length(neg_idx)
            neg_train = spktrain(neg_idx(k1), :)'; 
            neg_train = reshape(neg_train, [2, 599966]);
            neg_train = sum(neg_train)';
            for k2 = 1:length(posi_idx)
                nplot = (k1 - 1) * length(posi_idx) + k2;
                nplot = mod(nplot, 16);
                if nplot == 1
                    figure
                    figuresetup2savepdf(30, 30)
                elseif nplot == 0
                    nplot = 16;
                end
                posi_train = spktrain(posi_idx(k2), :)';
                posi_train = reshape(posi_train, [2, 599966]);
                posi_train = sum(posi_train)';
                xc = xcorr(posi_train, neg_train, 200);
                c = corr(posi_train, neg_train);
                subplot(4,4, nplot)
                bar(-200:1:200, xc, 1, 'FaceColor', 'k')
                hold on
                title(sprintf('#%d-#%d CC: %.3f', neg_idx(k1), posi_idx(k2), c))
                if max(xc) < 1
                    ylim([-1 1])
                end
                plot([-200 200], mean(xc)* [1 1], 'b')
                plot([-200 200], (mean(xc) + 3*std(xc))* [1 1], 'g')
                plot([-200 200], (mean(xc) - 3*std(xc))* [1 1], 'g')
                if nplot == 16
                    printPDFandPSC(gcf, fullfile(ccgfolder, sprintf('ccg_%s_cNE%d_%d', exp, jj, nfigure)))
                    nfigure = nfigure + 1;
                    close
                end
            end
        end
        if nplot ~= 16
            printPDFandPSC(gcf, fullfile(ccgfolder, sprintf('ccg_%s_cNE%d_%d', exp, jj,nfigure)))
            nfigure = nfigure + 1;
            close
        end
    end
end