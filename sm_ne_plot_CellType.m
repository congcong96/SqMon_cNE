%% set path
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

datapath = '/data/congcong/SqMoPhys_Josh/cNE_analysis';
spkfolder = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh';
cd(datapath)
nefiles = dir('*-20dft.mat');

%% find the spike waveforms and determine trough to peak delay
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/singleunit/waveform';
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    exp = exp_site_nedata.exp;
    spkfile = dir(fullfile(spkfolder, [exp, '*']));
    load(fullfile(spkfile.folder, spkfile.name), 'waveform_STRF_CRH', 'spk')
    fs = spk.fs;
    probe = spk.probe;
    if contains(probe, 'Tetrode')
        if strcmp(spk.probe, 'TetrodeB1x64')
            prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/TetrodeB1x64.csv';
        elseif strcmp(spk.probe, 'Tetrode1x64')
            prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/Tetrode1x64.csv';
        end
    elseif strcmp(spk.probe, 'poly21x48')
        prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/poly21x48.csv';
    end
    prblayout = readmatrix(prbfile);
    waveform = struct('IntanChan', [], 'position', [], 'average', [], 'std', [], ...
        't_trough', [], 't_peak', [], 'tpd', []);
    for jj = 1:length(waveform_STRF_CRH)
        chan = waveform_STRF_CRH(jj).chan;
        position = abs(spk.spk(jj).position);
        waveform(jj).IntanChan = chan;
        waveform(jj).position = position;
        average = waveform_STRF_CRH(jj).avgwaveform(chan+1,:);
        average_interp = interp1(0:0.05:2,average, 0:0.01:2, 'spline');
        waveform(jj).average = average_interp;
        waveform(jj).std = waveform_STRF_CRH(jj).stdwaveform(chan+1,:);
        waveform(jj).n0 = waveform_STRF_CRH(jj).n0contra;
        
        % get trough peak difference
        [m1,trough_idx] = min(waveform(jj).average);
        t_trough = (trough_idx-1)/fs*1000/5; 
        waveform(jj).t_trough = t_trough;
        [m2,peak_idx] = max(waveform(jj).average(trough_idx:end));
        t_peak = (peak_idx+trough_idx-2)/fs*1000/5; 
        waveform(jj).t_peak = t_peak;
        tpd = waveform(jj).t_peak - waveform(jj).t_trough;
        waveform(jj).tpd = tpd;
        % plot waveform for each neurons
        figure
        figuresetup2savepdf(15, 15)
        curve1 = waveform(jj).average(1:5:end) + waveform(jj).std;
        curve2 = waveform(jj).average(1:5:end) - waveform(jj).std;
        time = (0:(length(waveform(jj).average)-1)/5)/fs*1000;
        if waveform(jj).tpd <= 0.35
            c1 = [0.4, 0.8, 1];
            c2 = 'b';
        else
            c1 = [1, 0.8, 0.8];
            c2 = 'r'; 
        end
        a = fill([time, fliplr(time)], [curve1, fliplr(curve2)], 'r');
        a.FaceColor = c1;
        a.EdgeColor = [1 1 1];
        hold on
        plot((0:(length(waveform(jj).average)-1))/fs*1000/5, waveform(jj).average, c2, 'linewidth', 5)
        xlim([0 2])
        y = [m1 m2];
        Ylim = ylim;
        plot(t_trough*[1 1], y, 'k--', 'LineWidth', 2)
        plot(t_peak*[1 1], y, 'k--', 'LineWidth', 2)
        annotation('doublearrow', [0.13+t_trough/2*0.77, 0.13+t_peak/2*0.77], (0.11+(m1-Ylim(1))/diff(ylim))*0.89*[1, 1], 'LineWidth', 3)
        text((t_peak+t_trough)/2- tpd*.25, m1 - 10, sprintf('%.2fms', tpd), 'fontsize', 13)
        box off
        saveas(gcf, fullfile(figurefolder, sprintf('waveform_%s_%d.jpg', spk.exp, jj)))
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('waveform_%s_%d', spk.exp, jj)))
        
        close
    end
    save(nefiles(ii).name, 'waveform', '-append')
end

%% plot spike waveforms for each recording
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/singleunit/waveform';
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'waveform', 'exp_site_nedata')
    
    nfigure = 1;
    ymin = floor(min(min([waveform.average]))/10)*10;
    ymax = ceil(max(max([waveform.average]))/10)*10;

    for jj = 1:length(waveform)
        
        nplot = mod(jj, 16);
        if nplot == 1
            figure
            figuresetup2savepdf(30, 30)
        elseif nplot == 0
            nplot = 16;
        end
        average = waveform(jj).average;
        std = waveform(jj).std;
        n0 = waveform(jj).n0;
        % get trough peak difference
        t_trough = waveform(jj).t_trough; 
        t_peak = waveform(jj).t_peak;
        tpd = t_peak - t_trough;
        
        % plot waveform 
        curve1 = average(1:5:end) + std;
        curve2 = average(1:5:end) - std;
        time = (0:(length(average)-1)/5)/fs*1000;
        if tpd <= 0.35
            c1 = [0.4, 0.8, 1];
            c2 = 'b';
        else
            c1 = [1, 0.8, 0.8];
            c2 = 'r'; 
        end
        subplot(4,4,nplot)
        a = fill([time, fliplr(time)], [curve1, fliplr(curve2)], 'r');
        a.FaceColor = c1;
        a.EdgeColor = [1 1 1];
        hold on
        plot((0:(length(waveform(jj).average)-1))/fs*1000/5, waveform(jj).average, c2, 'linewidth', 2)
        xlim([0 2])
        ylim(1.2*[ymin, ymax]);
        Ylim = ylim;
        axis off
        title(sprintf('%.2fms %dspk', tpd, n0 ))
        if nplot == 16
            printPDFandPSC(gcf, fullfile(figurefolder, sprintf('waveform_all_%s_%d', exp_site_nedata.exp,nfigure)))
            nfigure = nfigure + 1;
            close
        end   
        
    end
    if nplot ~= 16
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('waveform_all_%s_%d', exp_site_nedata.exp,nfigure)))
        nfigure = nfigure + 1;
    end
    close
end

%% Threshold set at 0.35ms: spike width vs. firing rate
% 1 for narrow spike and 0 for broad spiking
tpd = cell(1, length(nefiles));
fr = cell(1, length(nefiles));
dur = 60*10;
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'waveform')
    tpd{ii} = [waveform.tpd];
    fr{ii} = [waveform.n0]/dur;
end
tpd = cell2mat(tpd);
fr = cell2mat(fr);
scatter(tpd(tpd <= 0.35), fr(tpd <= 0.35), 'b')
hold on
scatter(tpd(tpd > 0.35), fr(tpd > 0.35), 'r')
xlabel('trough to peak delay (ms)')
ylabel('firing rate (Hz)')
%% Threshold set at 0.35ms and label neuron types
% 1 for narrow spike and 0 for broad spiking
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/celltype';
tpd = cell(1, length(nefiles));
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata','waveform')
    tpd{ii} = [waveform.tpd];
    celltype = [waveform.tpd] <= 0.35;
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns, 2);
    NEmember = nedata.NEmembers;
    for jj = 1:nNE
        % stem plot
        pattern = nedata.Patterns(:,jj);
        member = NEmember{jj};
        
        posi_member = member(pattern(member) > 0);
        neg_member = member(pattern(member) < 0);
        
        figure
        figuresetup2savepdf(15, 5)
        stem(pattern, 'k','linewidth', 2)
        hold on
        stem(posi_member, pattern(posi_member), 'r', 'linewidth', 2)
        stem(neg_member, pattern(neg_member), 'b', 'linewidth', 2)
        xlim([0 length(pattern)+1])
        ylabel('ICweight')
        ylim([min([0, min(pattern)-0.2]), max(pattern) + 0.2])
        
        idx = find(celltype & pattern' > 0);
        scatter(idx, pattern(idx) + 0.1, '*b')
        idx = find(celltype & pattern' < 0);
        scatter(idx, pattern(idx) - 0.1, '*b')
        title(sprintf('%s    cNE#%d', exp_site_nedata.exp, jj), 'interpreter', 'none')
        saveas(gcf, fullfile(figurefolder, sprintf('celltype_%s_%d.jpg', exp_site_nedata.exp, jj)))
        printPDFandPSC(gcf, fullfile(figurefolder, sprintf('celltype_%s_%d', exp_site_nedata.exp, jj)));
        
        close
    end
end
%% histogram of distribution of waveform length
tpd = cell2mat(tpd);
edges = 0:0.1:1;
counts = histcounts(tpd, edges);
centers = (edges(1:end-1) + edges(2:end))/2;
bar(centers(1:3), counts(1:3), 1, 'b')
hold on
bar(centers(4)-0.025, counts(4), 0.05, 'b')
bar(centers(5:end), counts(5:end), 1, 'r')
bar(centers(4)+0.025, counts(4), 0.05, 'r')
plot(0.35*[1 1], ylim, 'k--')
xlim([-0.1 1.1])
xlabel('trough to peak delay (ms)')
ylabel('number of neurons')
%% histogram of the number of cNE a neuron belongs to
neuron_nNE = cell(1,length(nefiles));
celltype = cell(1,length(nefiles));
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata', 'waveform')
    nedata = exp_site_nedata.nedata;
    celltype{ii} = [waveform.tpd] <= 0.35;
    member = nedata.NEmembers;
    neuron_nNE{ii} = zeros(1, size(nedata.spktrain,1));
    member = cell2mat(member');
    [n, idx] = groupcounts(member);
    neuron_nNE{ii}(idx) = n;
end
%
celltype = cell2mat(celltype);
neuron_nNE = cell2mat(neuron_nNE);
c_ns = histcounts(neuron_nNE(celltype), 0:1:5, 'normalization', 'probability');
c_bs = histcounts(neuron_nNE(~celltype), 0:1:5, 'normalization', 'probability');
bar((0:4)-0.25, c_ns, 0.5, 'FaceColor', 'b')
hold on
bar((0:4)+0.25, c_bs, 0.5, 'FaceColor', 'r')   
legend({'narrow spiking', 'broad spiking'})
xlabel('number of cNE a neuron belongs to')
ylabel('proportion')

%% bar plot of proportion of inhibitory neurons in positive members and negative members
membertype = cell(1,length(nefiles));
celltype = cell(1,length(nefiles));
for ii = 2:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata', 'waveform')
    nedata = exp_site_nedata.nedata;
    celltype_all = double([waveform.tpd] <= 0.35);
    nNE = size(nedata.Patterns,2);
    NEmember = nedata.NEmembers;
    membertype_tmp = cell(1, nNE);
    celltype_tmp = cell(1, nNE);
    for jj = 1:nNE
        member = NEmember{jj};
        pattern = nedata.Patterns(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        
        if isempty(neg_idx) || length(posi_idx) == 1
            continue
        end
        
        membertype_tmp{jj} = [ones(size(posi_idx)); zeros(size(neg_idx))];
        celltype_tmp{jj} = [celltype_all(posi_idx), celltype_all(neg_idx)];
    end
    
    membertype{ii} = cell2mat(membertype_tmp');
    celltype{ii} = cell2mat(celltype_tmp);
end

celltype = cell2mat(celltype);
membertype = cell2mat(membertype');
%
c_pNE = histcounts(celltype(membertype==0), 0:2, 'normalization', 'probability');
c_nNE = histcounts(celltype(membertype==1), 0:2, 'normalization', 'probability');
c_all = histcounts(celltype, 0:2, 'normalization', 'probability');
bar( [0 2]-0.5, c_pNE, 0.25, 'FaceColor', 'r')
hold on
bar([0 2], c_nNE, 0.25, 'FaceColor', 'b')   
bar([0 2]+0.5, c_all, 0.25, 'FaceColor', 'k')   
legend({'positive NEmember', 'Negative NEmember', 'all neurons'})
xticks([0 2])
xticklabels({'broad spiking', 'narrow spiking'})
ylabel('proportion')

%% histogram of proportion of inhibitory neurons in positive members and negative members
proportion_ne = cell(1,length(nefiles));
proportion_posi = cell(1,length(nefiles));
proportion_neg = cell(1,length(nefiles));
proportion_all = cell(1,length(nefiles));
netype = cell(1,length(nefiles));
file_idx = cell(1,length(nefiles));
ne_idx = cell(1,length(nefiles));
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata', 'waveform')
    nedata = exp_site_nedata.nedata;
    celltype_all = double([waveform.tpd] <= 0.35);
    nNE = size(nedata.Patterns,2);
    proportion_all{ii} = sum(celltype_all)/length(celltype_all) * ones(1, nNE);
    NEmember = nedata.NEmembers;
    proportion_ne_tmp = zeros(1, nNE);
    proportion_posi_tmp = zeros(1, nNE);
    proportion_neg_tmp = zeros(1, nNE);
    netype_tmp = zeros(1, nNE);
    file_idx{ii} = ones(1, nNE) * ii;
    ne_idx{ii} = 1:nNE;
    for jj = 1:nNE
        member = NEmember{jj};
        pattern = nedata.Patterns(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        
        if isempty(neg_idx) 
            netype_tmp(jj) = 1;% type1:posi only
        elseif length(posi_idx) == 1
            netype_tmp(jj) = 3;% type3:only 1 pois, error
        else
            netype_tmp(jj) = 2;% type2:posi_neg
        end
        
        proportion_ne_tmp(jj) = sum(celltype_all(member))/length(celltype_all(member));
        proportion_posi_tmp(jj) = sum(celltype_all(posi_idx))/length(celltype_all(posi_idx));
        proportion_neg_tmp(jj) = sum(celltype_all(neg_idx))/length(celltype_all(neg_idx));

    end
    
    proportion_ne{ii} = proportion_ne_tmp;
    proportion_posi{ii} = proportion_posi_tmp;
    proportion_neg{ii} = proportion_neg_tmp;
    netype{ii} = netype_tmp;
    
end
% delete error cNEs
proportion_all = cell2mat(proportion_all)';
proportion_ne = cell2mat(proportion_ne)';
proportion_posi = cell2mat(proportion_posi)';
proportion_neg = cell2mat(proportion_neg)';
file_idx = cell2mat(file_idx)';
ne_idx = cell2mat(ne_idx)';
netype = cell2mat(netype)';
idx_error = find(netype == 3);
proportion_all(idx_error) = [];
proportion_ne(idx_error) = [];
proportion_posi(idx_error) = [];
proportion_neg(idx_error) = [];
file_idx(idx_error) = [];
ne_idx(idx_error) = [];
netype(idx_error) = [];
%%
subplot(311)
histogram(proportion_all, 0:0.05:1, 'FaceColor','k', 'FaceAlpha', 0.8)
ylabel('counts')
title('all neurons in the recording')
subplot(312)
histogram(proportion_ne, 0:0.05:1, 'FaceColor','m', 'FaceAlpha', 0.8)
ylabel('counts')
title('cNE members')
subplot(313)
histogram(proportion_posi, 0:0.05:1, 'FaceColor','r', 'FaceAlpha', 0.8)
ylabel('counts')
title('positive cNE members')
xlabel('proportion of NS neurons')
%% scatter plot
plot([0 1], [0 1], 'k--')
proportion_ne2 = proportion_ne(netype == 2);
proportion_posi2 = proportion_posi(netype == 2);
proportion_neg2 = proportion_neg(netype == 2);
hold on
scatter(proportion_ne2, proportion_posi2, 'r','filled')
xlabel('proportions of NS neurons in cNE members')
ylabel('proportions of NS neurons in positive cNE members')
signrank(proportion_ne2, proportion_posi2)

figure
plot([0 1], [0 1], 'k--')
hold on
scatter(proportion_ne2, proportion_neg2, 'b','filled')
xlabel('proportions of NS neurons in cNE members')
ylabel('proportions of NS neurons in negative cNE members')
signrank(proportion_ne2, proportion_neg2)
%% firing rate of NS neurons (hist, NS proportion > 0.2 or not)
proportion_all = zeros(1,length(nefiles));
fr = cell(1,length(nefiles));
file_idx = cell(1,length(nefiles));
ne_idx = cell(1,length(nefiles));
nNE = zeros(1,length(nefiles));
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'waveform', 'exp_site_nedata')
    celltype_all = [waveform.tpd] <= 0.35;
    proportion_all(ii) = sum(celltype_all)/length(celltype_all);
    fr{ii} = [waveform(celltype_all).n0]/600;
    nNE(ii) = size(exp_site_nedata.nedata.Patterns,2);
end
fr_less = cell2mat(fr(proportion_all<=0.2));    
fr_more = cell2mat(fr(proportion_all>0.2 ));    
histogram(proportion_all, 0:0.05:1, 'Facecolor', 'k', 'FaceAlpha', 0.8)
xlabel('Proportion of NS neurons')
histogram(fr_less,0:0.5:20, 'normalization', 'probability')
hold on
histogram(fr_more,0:0.5:20,'normalization', 'probability')
legend({'NS <= 20% of total neurons', 'NS > 20% of total neurons'})
xlabel('Firing rate of NS neurons')
ylabel('porportion')
ranksum(fr_less, fr_more)
%%

