%% configure paths
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

data_path = '/data/congcong/SqMoPhys_Josh/cNE_analysis';
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/strf_rtf_crh_ccg/NE';

stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
stimfile = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');
mtffile = regexprep(stimfile, '_stim', '_mtf');

%% plot stem plot, CCG, STRF and CRH
nefiles = dir('*20dft.mat');

for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns,2);
    spktrain = nedata.spktrain;
    NEmember = nedata.NEmembers;
    for jj = 1:nNE
        pattern = nedata.Patterns(:,jj);
        member = NEmember{jj};
        if ~any(pattern(member) < 0)
            continue
        end
        % get positive and negative member neuron spike trains
        posi_member = member(pattern(member) > 0);
        neg_member = member(pattern(member) < 0);
        spktrain_neg = spktrain(neg_member,:);
        spktrain_posi = spktrain(posi_member,:);
        
        % get positive and negative NE activities
        NEtrain_posi = nedata.NEtrain.posi(jj,:);
        NEtrain_posi = reshape(NEtrain_posi, [2, length(NEtrain_posi)/2]);
        NEtrain_posi = sum(NEtrain_posi);
        if sum(pattern(member) < 0) > 1
            NEtrain_neg = nedata.NEtrain.neg(jj,:);
            negNE = 1;
        else
            NEtrain_neg = spktrain_neg;
            negNE = 0;
        end
        
        % plot figure setup
        figure
        figuresetup2savepdf(15, 20)
        
        % stem plot
        subplot(4, 3,[1 2 3])
        stem(pattern, 'k','linewidth', 2)
        hold on
        stem(posi_member, pattern(posi_member), 'r', 'linewidth', 2)
        stem(neg_member, pattern(neg_member), 'b', 'linewidth', 2)
        xlim([0 length(pattern)+1])
        ylabel('ICweight')
        
        % ccg
        % all positive neurons vs. all negative neurons
        subplot(4,3,4)
        ccg_neuron = xcorr(sum(spktrain_neg,1), sum(spktrain_posi,1), 20);
        bar(-200:10:200, ccg_neuron, 1, 'k')
        xlim([-200 200])
        if max(ylim) < 1
            ylim([-1 1])
        end
        title('nCell vs. pCell')
        
        % all negative neurons vs. positive NEtrain
        subplot(4,3,5)
        ccg_nneuron_pNE = xcorr(sum(spktrain_neg,1), NEtrain_posi, 20);
        bar(-200:10:200, ccg_nneuron_pNE, 1, 'k')
        xlim([-200 200])
        if max(ylim) < 1
            ylim([-1 1])
        end
        title('nCell vs. pNE')
        
        % negative NEtrain vs. positive NEtrain
        if negNE % if exist multiple negative members
            subplot(4,3,6)
            ccg_nNE_pNE = xcorr(NEtrain_neg, NEtrain_posi, 20);
            bar(-200:10:200, ccg_nNE_pNE, 1, 'k')
            xlim([-200 200])
            if max(ylim) < 1
                ylim([-1 1])
            end
        end
        title('nNE vs. pNE')
        
        % strf
        taxis = nedata.strf_taxis;
        faxis = nedata.strf_faxis;
        subplot(4,3,7)
        plot_strf_raw(nedata.NEstrf_2018(jj,:), faxis, taxis)
        title('Paper2018')
        subplot(4,3,8)
        plot_strf_raw(nedata.NEstrf.posi(jj,:), faxis, taxis)
        title('NEposi')
        subplot(4,3,9)
        if negNE
            plot_strf_raw(nedata.NEstrf.neg(jj,:), faxis, taxis)
        else
            plot_strf_raw(nedata.stamat(neg_member,:), faxis, taxis)
        end
        title('NEneg')
        
        % crh
        taxis = nedata.rtf_tmf;
        faxis = nedata.rtf_smf;
        subplot(4, 3, 10)
        plot_CRH(nedata.NEcrh_2018(jj,:), taxis, faxis)
        title(sprintf('%dspk', sum(nedata.sta_NEtrain_2018(jj,:))))
        subplot(4, 3, 11)
        plot_CRH(nedata.NEcrh.posi(jj,:), taxis, faxis)
        title(sprintf('%dspk',  sum(nedata.NEtrain.posi(jj,:))))
        subplot(4, 3, 12)
        if negNE
            plot_CRH(nedata.NEcrh.neg(jj,:), taxis, faxis)
            title(sprintf('%dspk', sum(nedata.NEtrain.neg(jj,:))))
        else
            plot_CRH(nedata.neuroncrh(neg_member,:), taxis, faxis)
            title(sprintf('%dspk', sum(nedata.spktrain(neg_member,:))))
        end
       
        sgtitle(sprintf('%s cNE#%d',exp_site_nedata.exp, jj),'interpreter','none')
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('%s-thresh%d-cNE%d-%dms', exp_site_nedata.exp,nedata.MemberThr*100, jj,exp_site_nedata.df/2)));
        close
    end
end
%% correlation for example ensembles
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/strf_rtf_crh_ccg/CCG';
load('191113_145808-dmr-10min-poly21x48-fs20000-pydict-spk-curated-thresh-ne-20dft.mat')
nedata = exp_site_nedata.nedata;
nNE = size(nedata.Patterns,2);
spktrain = nedata.spktrain;
NEmember = nedata.NEmembers;
% index of the cNE
jj = 1;
pattern = nedata.Patterns(:,jj);
member = NEmember{jj};
% get positive and negative member neuron spike trains
posi_member = member(pattern(member) > 0);
neg_member = member(pattern(member) < 0);
spktrain_neg = spktrain(neg_member,:);
spktrain_posi = spktrain(posi_member,:);
nplot = 1;
nfig = 1;
for ii = 1:size(spktrain_neg,1)
    for jj = 1:size(spktrain_posi,1)
        ccg = xcorr(spktrain_neg(ii,:), spktrain_posi(jj,:), 20);
        c = corr(spktrain_neg(ii,:)', spktrain_posi(jj,:)');
        if nplot == 1
            figure
            figuresetup2savepdf(20, 20)
        elseif nplot == 10
            nplot = 1;
        end
        subplot(3, 3, nplot)
        bar(-200:10:200, ccg, 1, 'k')
        xlim([-200 200])
        if max(ylim) < 1
            ylim([-1 1])
        end
        title(sprintf('#%d:#%d (%.3f)', neg_member(ii), posi_member(jj), c))
        nplot = nplot+1;
        if nplot == 9
            printPDFandPSC(gcf, fullfile(figurefolder, sprintf('%s-thresh%d-cNE%d-%dms-%d', exp_site_nedata.exp,nedata.MemberThr*100, jj,exp_site_nedata.df/2, nfig)));
            nfig = nfig+1;
            close
        end
    end
end
if nplot ~= 9
    printPDFandPSC(gcf, fullfile(figurefolder, sprintf('%s-thresh%d-cNE%d-%dms-%d', exp_site_nedata.exp,nedata.MemberThr*100, jj,exp_site_nedata.df/2, nfig)));
    close
end
%% plot STRF and CRH for individual neurons
nefiles = dir('*20dft.mat');
% load stimulus
if ~exist('stim_mat', 'var')
    load(stimfile, 'stim_mat')
end
stimulus = stim_mat(:, 10:10:end);

for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns,2);
    spktrain = nedata.spktrain;
    NEmember = nedata.NEmembers;
    stamat = nedata.stamat;
    neuroncrh = nedata.neuroncrh;
    Thr = nedata.MemberThr;
    for jj = 1:nNE
        pattern = nedata.Patterns(:,jj);
        member = NEmember{jj};
        
        if ~any(pattern(member) < 0)
            continue
        end
        
        % get coincidence train of neurons and cNEs
        NEtrain = nedata.NEtrain.all(jj,:);
        NEtrain(NEtrain > 0) = 1;
        NEtrain = reshape(NEtrain, [2, length(NEtrain)/2]);
        NEtrain = sum(NEtrain);
        
        spktrain_coincide = zeros(size(spktrain,1), size(spktrain,2)*2);
        for kk = 1:size(spktrain,1)
            spktrain_coincide_tmp = zeros(size(NEtrain));
            spktrain_coincide_tmp(NEtrain > 0) = spktrain(kk, NEtrain > 0);
            spktrain_coincide(kk,:) = ne_upsample_NEact_using_member_neuron_activity(spktrain_coincide_tmp, {member}, nedata.sta_spktrain, 0);;
        end
        
        % calculate strf
        lendiff = size(stimulus, 2) - size(spktrain_coincide, 2);
        stimulus(:,end-lendiff+1:end) = [];
        
        stamat_coincide = quick_calc_sta(stimulus, spktrain_coincide, 'nlags', 20, 'chunks', 1, 'suppressprint', 1);
        % calculate crh
        [mtfhist, tmf, smf] = sm_calculate_CRH(spktrain_coincide, mtffile);
        
        % plot
        [~, neuronidx] = sort(pattern, 'descend');% plot in the descending ICweight
        nfigure = 1;
        for kk = 1:size(spktrain,1)
            nplot = mod(kk,4);
            if nplot == 1
                figure
                figuresetup2savepdf(20, 20)
            elseif nplot == 0
                nplot = 4;
            end
            idx = neuronidx(kk);
            % strf
            subplot(4,4,nplot*4-3)
            hold on
            taxis = nedata.strf_taxis;
            faxis = nedata.strf_faxis;
            plot_strf_raw(stamat(idx,:), faxis, taxis)
            title(sprintf('Neuron#%d spk%d', idx, sum(spktrain(idx,:))));
            % strf_coincide
            subplot(4,4,nplot*4-2)
            hold on
            plot_strf_raw(stamat_coincide(idx,:), faxis, taxis)
            title(sprintf('spk%d', sum(spktrain_coincide(idx,:))));
            % crh
            subplot(4,4,nplot*4-1)
            hold on
            taxis = nedata.rtf_tmf;
            faxis = nedata.rtf_smf;
            plot_CRH(neuroncrh(idx,:), taxis, faxis)
            % crh_coincide
            subplot(4,4,nplot*4)
            hold on
            taxis = nedata.rtf_tmf;
            faxis = nedata.rtf_smf;
            plot_CRH(mtfhist(idx,:), taxis, faxis)
            
            if any(member == idx)
                
                if pattern(idx) > 0
                    c = 'r';
                else
                    c = 'b';
                end
                % strf
                subplot(4,4,nplot*4-3)
                plot([0, 20], [0 0], c, 'linewidth', 3)
                plot([0, 20], [69 69], c, 'linewidth', 3)
                plot([0, 0], [0 69], c, 'linewidth', 3)
                plot([20, 20], [0 69], c, 'linewidth', 3)
                xlim([0 20])
                ylim([0 69])
                % strf_coincide
                subplot(4,4,nplot*4-2)
                plot([0, 20], [0 0], c, 'linewidth', 3)
                plot([0, 20], [69 69], c, 'linewidth', 3)
                plot([0, 0], [0 69], c, 'linewidth', 3)
                plot([20, 20], [0 69], c, 'linewidth', 3)
                xlim([0 20])
                ylim([0 69])
                % crh
                subplot(4,4,nplot*4-1)
                plot([-64, 64], [0 0], c, 'linewidth', 3)
                plot([-64, 64], [4 4], c, 'linewidth', 3)
                plot([-64, -64], [0 4], c, 'linewidth', 3)
                plot([64, 64], [0 4], c, 'linewidth', 3)
                xlim([-64 64])
                ylim([0 4])
                % crh_coincide
                subplot(4,4,nplot*4)
                plot([-64, 64], [0 0], c, 'linewidth', 3)
                plot([-64, 64], [4 4], c, 'linewidth', 3)
                plot([-64, -64], [0 4], c, 'linewidth', 3)
                plot([64, 64], [0 4], c, 'linewidth', 3)
                xlim([-64 64])
                ylim([0 4])
            end
            if nplot == 4
                printPDFandPSC(gcf, fullfile(figurefolder, sprintf('%s-thresh%d-cNE%d-%dms-neuron-%d', exp_site_nedata.exp,Thr*100, jj,exp_site_nedata.df/2, nfigure)));
                nfigure = nfigure + 1;
                close
            end
        end
        if nplot ~= 4
            printPDFandPSC(gcf, fullfile(figurefolder, sprintf('%s-thresh%d-cNE%d-%dms-neuron-%d', exp_site_nedata.exp,Thr*100, jj,exp_site_nedata.df/2, nfigure)));
            close
        end
    end
end