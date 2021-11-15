%% set path
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

datapath = '/data/congcong/SqMoPhys_Josh/cNE_analysis';
cd(datapath)
nefiles = dir('*-20dft.mat');

stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
stimfile = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');
mtfstimfile = regexprep(stimfile, '_stim', '_mtf');
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/MoranI';
%% calculate Moran's I
nrepeat = 100; % number of repeat for calculating subsampled Moran's I
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata');
    nedata = exp_site_nedata.nedata;
    
    % get weigth matrix for calculating Moron's i
    if ~exist('weightsmat', 'var')
        nrows = length(nedata.crh_smf);
        ncols = length(nedata.crh_tmf);
        weightsmat = create_rf_spatial_autocorr_weights_matrix(nrows, ncols, 0);
        weightsmat_diag = create_rf_spatial_autocorr_weights_matrix(nrows, ncols, 1);
    end
    
    fprintf('Processing %s\n', nefiles(ii).name)
    
    % calculate Moran's I for individual neurons
    nNeuron  = size(nedata.neuroncrh,1);
    I_neuron = zeros(1, nNeuron);
    I_neuron_diag = zeros(1, nNeuron);
    for jj = 1:nNeuron
        I_neuron(jj) = morans_i(nedata.neuroncrh(jj,:), weightsmat);
        I_neuron_diag(jj) = morans_i(nedata.neuroncrh(jj,:), weightsmat_diag);
    end
    nedata.MoranI.neuron = I_neuron;
    nedata.MoranI_diag.neuron = I_neuron_diag;
    
    % calculate Moran's I for positive NEtrain
    nNE = size(nedata.Patterns, 2);
    I_pNE = zeros(1, nNE);
    I_pNE_diag = zeros(1, nNE);
    for jj = 1:nNE
        I_pNE(jj) = morans_i(nedata.NEcrh.posi(jj,:), weightsmat);
        I_pNE_diag(jj) = morans_i(nedata.NEcrh.posi(jj,:), weightsmat_diag);
    end
    nedata.MoranI.pNE = I_pNE;
    nedata.MoranI_diag.pNE = I_pNE_diag;
    
    % zeros matrices to hold subsampled moran's I value
    I_pNE_subsample_neuron = zeros(nNE, nrepeat);
    I_neuron_subsample = cell(1, nNE);
    I_pNE_subsample_neuron_diag = zeros(nNE, nrepeat);
    I_neuron_subsample_diag = cell(1, nNE);
    
    I_pNE_nonneg = zeros(1, nNE);
    I_pNE_subsample_nonneg = zeros(nNE, nrepeat);
    I_pNE_nonneg_diag = zeros(1, nNE);
    I_pNE_subsample_nonneg_diag = zeros(nNE, nrepeat);
    for jj = 1:nNE
        
        member = nedata.NEmembers{jj};
        pattern = nedata.Patterns(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        NEtrain_posi = nedata.NEtrain.posi(jj,:);
        spktrain_posi = nedata.spktrain(posi_idx,:);
        spktrain_neg = nedata.spktrain(neg_idx,:);

        % calculate Moran's I with subsampling
        nevent = min([sum(NEtrain_posi); sum(spktrain_posi,2)]);
        if nevent == 0
            continue
        end
        % get NE Moran's I after subsampling
        I = sm_ne_calc_crh_Morans_I_subsample(NEtrain_posi, mtfstimfile, {weightsmat, weightsmat_diag}, nevent, nrepeat);
        I_pNE_subsample_neuron(jj,:) = I(1,:);
        I_pNE_subsample_neuron_diag(jj,:) = I(2,:);
        % get positve member's  Moran's I after subsampling
        I_neuron_subsample_tmp = zeros(length(posi_idx), nrepeat);
        I_neuron_subsample_diag_tmp = zeros(length(posi_idx), nrepeat);
        for kk = 1:length(posi_idx)
            I = sm_ne_calc_crh_Morans_I_subsample(spktrain_posi(kk,:), mtfstimfile, {weightsmat, weightsmat_diag}, nevent, nrepeat);
            I_neuron_subsample_tmp(kk,:) = I(1,:);
            I_neuron_subsample_diag_tmp(kk,:) = I(2,:);
        end
        I_neuron_subsample{jj} = I_neuron_subsample_tmp;
        I_neuron_subsample_diag{jj} = I_neuron_subsample_diag_tmp;
        
        % get NE train excluding coincidence with negative spikes
        NEtrain_posi_original = sum(reshape(NEtrain_posi, [2, length(NEtrain_posi)/2]));
        NEtrain_posi_nonneg_idx = find( NEtrain_posi_original & ~sum(spktrain_neg,1));
        NEtrain_posi_nonneg_idx = [NEtrain_posi_nonneg_idx*2, NEtrain_posi_nonneg_idx*2 - 1];
        NEtrain_posi_nonneg = zeros(size(NEtrain_posi));
        NEtrain_posi_nonneg(NEtrain_posi_nonneg_idx) = NEtrain_posi(NEtrain_posi_nonneg_idx);
        % calculate crh for NE train excluding coincidence with negative spikes
        if sum(NEtrain_posi_nonneg) < sum(NEtrain_posi)
            crh = sm_calculate_CRH(NEtrain_posi_nonneg, mtfstimfile);
            I_pNE_nonneg(jj) = morans_i(crh, weightsmat);
            I_pNE_nonneg_diag(jj) = morans_i(crh, weightsmat_diag);
            % get Monran's I for subsampled NEtrain
            I = sm_ne_calc_crh_Morans_I_subsample(NEtrain_posi, mtfstimfile, {weightsmat, weightsmat_diag}, sum(NEtrain_posi_nonneg), nrepeat);
            I_pNE_subsample_nonneg(jj,:) = I(1,:);
            I_pNE_subsample_nonneg_diag(jj,:) = I(2,:);
        else
            I_pNE_nonneg(jj) = I_pNE(jj);
            I_pNE_nonneg_diag(jj) = I_pNE_diag(jj);
            I_pNE_subsample_nonneg(jj,:) = I_pNE_nonneg(jj);
            I_pNE_subsample_nonneg_diag(jj,:) = I_pNE_nonneg_diag(jj);
        end
    end
    
    nedata.MoranI.pNeuron_subsample = I_neuron_subsample;
    nedata.MoranI_diag.pNeuron_subsample = I_neuron_subsample_diag;
    nedata.MoranI.pNE_subsample_pNeuron = I_pNE_subsample_neuron;
    nedata.MoranI_diag.pNE_subsample_pNeuron = I_pNE_subsample_nonneg;
    
    nedata.MoranI.pNE_nonneg = I_pNE_nonneg;
    nedata.MoranI.pNE_subsample_nonneg = I_pNE_subsample_nonneg;
    nedata.MoranI_diag.pNE_nonneg = I_pNE_nonneg_diag;
    nedata.MoranI_diag.pNE_subsample_nonneg = I_pNE_subsample_nonneg_diag;
    
    exp_site_nedata.nedata = nedata;
    save(nefiles(ii).name, 'exp_site_nedata', '-append')
end






%% scatter plot 1 - neuron v.s. cNE (subsample, mean)
I_NE_subsample = cell(size(nefiles));
I_neuron_subsample = cell(size(nefiles));
file_number = cell(size(nefiles));
NEnumber = cell(size(nefiles));
NEtype = cell(size(nefiles));
MI_method = 'diag';
average_method = 'median';
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    if strcmp(MI_method, 'diag') && strcmp(average_method, 'mean')
        I_NE_subsample{ii} = mean(nedata.MoranI_diag.pNE_subsample_pNeuron, 2);
        I_neuron_subsample{ii} = cellfun(@(x) mean(mean(x, 2)), nedata.MoranI_diag.pNeuron_subsample);
    elseif strcmp(MI_method, 'diag') && strcmp(average_method, 'median')
        I_NE_subsample{ii} = median(nedata.MoranI_diag.pNE_subsample_pNeuron, 2);
        I_neuron_subsample{ii} = cellfun(@(x) median(median(x, 2)), nedata.MoranI_diag.pNeuron_subsample);
    elseif strcmp(MI_method, 'border') && strcmp(average_method, 'mean')
        I_NE_subsample{ii} = mean(nedata.MoranI.pNE_subsample_pNeuron, 2);
        I_neuron_subsample{ii} = cellfun(@(x) mean(mean(x, 2)), nedata.MoranI.pNeuron_subsample);
    elseif strcmp(MI_method, 'border') && strcmp(average_method, 'median')
        I_NE_subsample{ii} = median(nedata.MoranI.pNE_subsample_pNeuron, 2);
        I_neuron_subsample{ii} = cellfun(@(x) median(median(x, 2)), nedata.MoranI.pNeuron_subsample);
    end
    file_number{ii} = ii * ones(size(I_NE_subsample{ii}));
    NEnumber{ii} = 1:length(I_NE_subsample{ii});
    NEtype{ii} = zeros(size(I_NE_subsample{ii}));
    for jj = 1:length(NEtype{ii})
        member = nedata.NEmembers{jj};
        pattern = nedata.Patterns(:,jj);
        if sum(pattern(member) < 0) == 1
            NEtype{ii}(jj) = 1;
        elseif sum(pattern(member) < 0) > 1
            NEtype{ii}(jj) = 2;
        end
            
    end
end

I_NE_subsample = cell2mat(I_NE_subsample);
I_neuron_subsample = cell2mat(I_neuron_subsample')';
file_number = cell2mat(file_number);
NE_number = cell2mat(NEnumber')';
NEtype = cell2mat(NEtype);
idx = ~isnan(I_neuron_subsample);
I_NE_subsample = I_NE_subsample(idx);
I_neuron_subsample = I_neuron_subsample(idx);
file_number = file_number(idx);
NE_number = NE_number(idx);
NEtype = NEtype(idx);

for ii = 1:length(I_NE_subsample)
    if isnan(I_neuron_subsample(ii))
        continue
    end
    figure
    scatter(I_neuron_subsample, I_NE_subsample, 'r');
    hold on
    xlabel(sprintf("%s Moran's I of positive member neurons", average_method))
    ylabel("Moran's I of positive NEtrain")
    scatter(I_neuron_subsample(NEtype > 0), I_NE_subsample(NEtype > 0), 'b');
    plot([-0.1 .9], [-0.1 .9], 'k--')
    xlim([-.1, .9])
    ylim([-.1, .9])
    % circle example NE
    scatter(I_neuron_subsample(ii), I_NE_subsample(ii), 100, 'k');
    legend({'cNE with only positive members', 'cNE with negative members'})
    filename = nefiles(file_number(ii)).name;
    saveas(gcf, fullfile(figurefolder, sprintf('MoranI_pNE_vs_pMember_%s_%s_%s_cNE%d', average_method, MI_method, filename(1:13), NE_number(ii))))
    saveas(gcf, fullfile(figurefolder, sprintf('MoranI_pNE_vs_pMember_%s_%s_%s_cNE%d.jpg', average_method, MI_method, filename(1:13), NE_number(ii))))
    close all
end
%% histogram 
NEneuron_diff = I_NE_subsample -  I_neuron_subsample;
edges = -0.2:0.05:0.8;
posi = histcounts(I_NE_subsample(NEtype == 0) - I_neuron_subsample(NEtype == 0), edges);
posi_neg = histcounts(I_NE_subsample(NEtype > 0) - I_neuron_subsample(NEtype > 0), edges);
centers = (edges(2:end) + edges(1:end-1))/2;
b = bar(centers, [posi', posi_neg'], 'stacked');
b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
hold on
scatter(mean(NEneuron_diff), 29, 'vk')
scatter(median(NEneuron_diff), 29, 'vk', 'filled')
scatter(mean(NEneuron_diff(NEtype > 0)), 28, 'vb')
scatter(median(NEneuron_diff(NEtype > 0)), 28, 'vb', 'filled')
scatter(mean(NEneuron_diff(NEtype == 0)), 27, 'vr')
scatter(median(NEneuron_diff(NEtype == 0)), 27, 'vr', 'filled')
legend({'cNE with only positive members', 'cNE with negative members'})
xlabel("pNE Moran's I - median Moran's I of positive member neurons")
ylabel("number of cNEs")
signtest(NEneuron_diff)
ranksum(NEneuron_diff(NEtype > 0), NEneuron_diff(NEtype == 0))
signrank(I_NE_subsample(NEtype == 0) ,  I_neuron_subsample(NEtype == 0))
%%  plot all cNE crh - MI
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/MoranI';
nefiles = dir('*20dft.mat');
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns,2);
    NEmember = nedata.NEmembers;
    for jj = 1:nNE
        pattern = nedata.Patterns(:,jj);
        member = NEmember{jj};

        posi_member = member(pattern(member) > 0);
        neg_member = member(pattern(member) < 0);
        
        % plot figure setup
        figure
        figuresetup2savepdf(15, 10)
        
        % stem plot
        subplot(2, 3,[1 2 3])
        stem(pattern, 'k','linewidth', 2)
        hold on
        stem(posi_member, pattern(posi_member), 'r', 'linewidth', 2)
        stem(neg_member, pattern(neg_member), 'b', 'linewidth', 2)
        xlim([0 length(pattern)+1])
        ylabel('ICweight')
        
        % crh
        taxis = nedata.rtf_tmf;
        faxis = nedata.rtf_smf;
        subplot(2, 3, 4)
        plot_CRH(nedata.NEcrh.posi(jj,:), taxis, faxis)
        MI = median(nedata.MoranI_diag.pNE_subsample_pNeuron(jj,:));
        title(sprintf('MI: %.3f',  MI))
       
        sgtitle(sprintf('%s cNE#%d',exp_site_nedata.exp, jj),'interpreter','none')
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('MoranI_%s-thresh%d-cNE%d-%dms', exp_site_nedata.exp,nedata.MemberThr*100, jj,exp_site_nedata.df/2)));
        close
    end
end

%%  plot all member CRH - MI
nefiles = dir('*20dft.mat');

for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns,2);
    NEmember = nedata.NEmembers;
    neuroncrh = nedata.neuroncrh;
    Thr = nedata.MemberThr;
    for jj = 1:nNE
        pattern = nedata.Patterns(:,jj);
        member = NEmember{jj};
        pmember = member(pattern(member) > 0);
        if isempty(nedata.MoranI_diag.pNeuron_subsample{jj})
            continue
        end
        % plot
        nfigure = 1;
        for kk = 1:length(pmember)
            nplot = mod(kk,12);
            if nplot == 1
                figure
                figuresetup2savepdf(15, 20)
            elseif nplot == 0
                nplot = 12;
            end
            idx = pmember(kk);
            % crh
            subplot(4,3,nplot)
            taxis = nedata.rtf_tmf;
            faxis = nedata.rtf_smf;
            plot_CRH(neuroncrh(idx,:), taxis, faxis)
            
            MI = median(nedata.MoranI_diag.pNeuron_subsample{jj}(kk,:));
            
            title(sprintf('#%d  MI: %.3f', idx, MI))
           
            if nplot == 12
                printPDFandPSC(gcf, fullfile(figurefolder, sprintf('MoranI_%s-thresh%d-cNE%d-%dms-neuron-%d', exp_site_nedata.exp,Thr*100, jj,exp_site_nedata.df/2, nfigure)));
                nfigure = nfigure + 1;
                close
            end
        end
        if nplot ~= 12
            printPDFandPSC(gcf, fullfile(figurefolder, sprintf('MoranI_%s-thresh%d-cNE%d-%dms-neuron-%d', exp_site_nedata.exp,Thr*100, jj,exp_site_nedata.df/2, nfigure)));
            close
        end
    end
end








%% scatter plot 2 - pNE v.s. pNE without neg spikes
I_NE_subsample = cell(size(nefiles));
I_nonnegNE = cell(size(nefiles));
file_number = cell(size(nefiles));
NEtype = cell(size(nefiles));
NEnumber = cell(size(nefiles));
MI_method = 'diag';
average_method = 'median';
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    if strcmp(MI_method, 'diag')
        I_nonnegNE{ii} = nedata.MoranI_diag.pNE_nonneg;
        if strcmp(average_method, 'mean')
            I_NE_subsample{ii} = mean(nedata.MoranI_diag.pNE_subsample_nonneg, 2);
        elseif strcmp(average_method, 'median')
            I_NE_subsample{ii} = median(nedata.MoranI_diag.pNE_subsample_nonneg, 2);
        end
    elseif strcmp(MI_method, 'border')
        I_nonnegNE{ii} = nedata.MoranI.pNE_nonneg;
        if strcmp(average_method, 'mean')
            I_NE_subsample{ii} = mean(nedata.MoranI.pNE_subsample_nonneg, 2);
        elseif strcmp(MI_method, 'border') && strcmp(average_method, 'median')
            I_NE_subsample{ii} = median(nedata.MoranI.pNE_subsample_nonneg, 2);
        end
    end
    
    file_number{ii} = ii * ones(size(I_NE_subsample{ii}));
    NEnumber{ii} = 1:length(I_NE_subsample{ii});
    NEtype{ii} = zeros(size(I_NE_subsample{ii}));
    for jj = 1:length(NEtype{ii})
        member = nedata.NEmembers{jj};
        pattern = nedata.Patterns(:,jj);
        if sum(pattern(member) < 0) == 1
            NEtype{ii}(jj) = 1;
        elseif sum(pattern(member) < 0) > 1
            NEtype{ii}(jj) = 2;
        end    
    end
end

I_NE_subsample = cell2mat(I_NE_subsample);
I_nonnegNE = cell2mat(I_nonnegNE')';
NE_number = cell2mat(NEnumber')';
file_number = cell2mat(file_number);
NEtype = cell2mat(NEtype);
idx = find(I_NE_subsample == 0);
I_NE_subsample(idx) = [];
I_nonnegNE(idx) = [];
NE_number(idx) = [];
file_number(idx) = [];
NEtype(idx) = [];
idx = find(NEtype == 0);
I_NE_subsample(idx) = [];
I_nonnegNE(idx) = [];
NE_number(idx) = [];
file_number(idx) = [];
NEtype(idx) = [];
%%
for ii = 1:length(I_NE_subsample)
    figure
    scatter(I_NE_subsample(NEtype > 0), I_nonnegNE(NEtype > 0), 'b');
    hold on
    plot([-0.1 .9], [-0.1 .9], 'k--')
    xlim([-.1, .9])
    ylim([-.1, .9])
    xlabel(" Moran's I of positive NEtrain")
    ylabel("Moran's I of positive NEtrain wihtout negative spikes")
    % circle example NE
    scatter(I_NE_subsample(ii), I_nonnegNE(ii), 100, 'k');
    filename = nefiles(file_number(ii)).name;
    saveas(gcf, fullfile(figurefolder, sprintf('MoranI_pNE_vs_nonnegNE_%s_%s_%s_cNE%d', average_method, MI_method, filename(1:13), NE_number(ii))))
    saveas(gcf, fullfile(figurefolder, sprintf('MoranI_pNE_vs_nonnegNE_%s_%s_%s_cNE%d.jpg', average_method, MI_method, filename(1:13), NE_number(ii))))
    close all
end

%% histogram 
NEneuron_diff = I_nonnegNE(NEtype > 0) -  I_NE_subsample(NEtype > 0);
histogram(NEneuron_diff, -0.04:0.005:0.04, 'FaceColor', 'b');
hold on
scatter(mean(NEneuron_diff), 23, 'vb')
scatter(median(NEneuron_diff), 23, 'vb', 'filled')
xlabel("pNE Moran's I - Moran's I of positive NEtrain without negative spikes")
ylabel("number of cNEs")
signtest(NEneuron_diff)
%%  plot all cNE crh - MI
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/MoranI';
nefiles = dir('*20dft.mat');
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns,2);
    NEmember = nedata.NEmembers;
    for jj = 1:nNE
        pattern = nedata.Patterns(:,jj);
        member = NEmember{jj};

        posi_member = member(pattern(member) > 0);
        neg_member = member(pattern(member) < 0);
        if isempty(neg_member)
            continue
        end
        % plot figure setup
        figure
        figuresetup2savepdf(15, 10)
        
        % stem plot
        subplot(2, 3,[1 2 3])
        stem(pattern, 'k','linewidth', 2)
        hold on
        stem(posi_member, pattern(posi_member), 'r', 'linewidth', 2)
        stem(neg_member, pattern(neg_member), 'b', 'linewidth', 2)
        xlim([0 length(pattern)+1])
        ylabel('ICweight')
        
        % posi crh
        taxis = nedata.rtf_tmf;
        faxis = nedata.rtf_smf;
        subplot(2, 3, 4)
        plot_CRH(nedata.NEcrh.posi(jj,:), taxis, faxis)
        MI = median(nedata.MoranI_diag.pNE_subsample_nonneg(jj,:));
        title(sprintf('posi MI: %.3f',  MI))
        % neg crh
        subplot(2, 3, 5)
        if length(neg_member) == 1
            NEtrain_neg = nedata.spktrain(neg_member,:);
            crh = nedata.neuroncrh(neg_member,:);
        else
            NEtrain_neg = nedata.NEtrain.neg(jj,:);
            crh = nedata.NEcrh.neg(jj,:);
        end
        plot_CRH(crh, taxis, faxis)
        title('neg')
        
        % posi - neg crh
        % get NE train excluding coincidence with negative spikes
        subplot(2, 3, 6)
        spktrain_neg = nedata.spktrain(neg_member,:);
        NEtrain_posi = nedata.NEtrain.posi(jj,:);
        NEtrain_posi_original = sum(reshape(NEtrain_posi, [2, length(NEtrain_posi)/2]));
        NEtrain_posi_nonneg_idx = find( NEtrain_posi_original & ~sum(spktrain_neg,1));
        NEtrain_posi_nonneg_idx = [NEtrain_posi_nonneg_idx*2, NEtrain_posi_nonneg_idx*2 - 1];
        NEtrain_posi_nonneg = zeros(size(NEtrain_posi));
        NEtrain_posi_nonneg(NEtrain_posi_nonneg_idx) = NEtrain_posi(NEtrain_posi_nonneg_idx);
        % calculate crh for NE train excluding coincidence with negative spikes
        crh = sm_calculate_CRH(NEtrain_posi_nonneg, mtfstimfile);
        plot_CRH(crh, taxis, faxis)
        MI = nedata.MoranI_diag.pNE_nonneg(jj);
        title(sprintf('MI: %.3f',  MI))
       
        sgtitle(sprintf('%s cNE#%d',exp_site_nedata.exp, jj),'interpreter','none')
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('MoranI_nonneg_%s-thresh%d-cNE%d-%dms', exp_site_nedata.exp,nedata.MemberThr*100, jj,exp_site_nedata.df/2)));
        close
    end
end