%% set path
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

datapath = '/data/congcong/SqMoPhys_Josh/cNE_analysis';
spkfolder = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh';
figfolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/position';
cd(datapath)
nefiles = dir('*poly*-20dft.mat');

%% plot layout of member neurons
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    cellPosition = nedata.position;
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    Pattern = nedata.Patterns;
    
    for jj = 1:nNE
        member = NEmembers{jj};
        pattern = Pattern(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        if isempty(neg_idx) || length(posi_idx) <=1
            continue
        end
        
        posi_position = cell2mat(cellPosition(posi_idx)');
        neg_position = cell2mat(cellPosition(neg_idx)');
        load(fullfile('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis', 'ElectrodePositions.mat'))
        
        if contains(exp_site_nedata.probetype, 'Tetrode')
            figure
            figuresetup2savepdf(20, 20)
            hold on
            plot(-20*[1 1], [50 600], 'k')
            plot(60*[1 1], [50 600], 'k')
            plot(180*[1 1], [50 600], 'k')
            plot(260*[1 1], [50 600], 'k')
            plot(380*[1 1], [50 600], 'k')
            plot(460*[1 1], [50 600], 'k')
            plot(580*[1 1], [50 600], 'k')
            plot(660*[1 1], [50 600], 'k')
            
            plot([-20, 20], [50 0], 'k')
            plot([20, 60], [0 50], 'k')
            plot([180, 220], [50 0], 'k')
            plot([220, 260], [0 50], 'k')
            plot([380, 420], [50 0], 'k')
            plot([420, 460], [0 50], 'k')
            plot([580, 620], [50 0], 'k')
            plot([620, 660], [0 50], 'k')
            
            x = ElectrodePositions(8).probmat(:,5);
            y = ElectrodePositions(8).probmat(:,1);
            
            scatter(x, y, 100, 'ok')
            xlim([-100 700])
            ylim([-100 700])
        else
            
            figure
            figuresetup2savepdf(10, 30)
            hold on
            plot(-50*[1 1], [50 2500], 'k')
            plot(150*[1 1], [50 2500], 'k')
            
            plot([-50, 50], [50 0], 'k')
            plot([50, 150], [0 50], 'k')
            
            x = ElectrodePositions(10).probmat(:,5);
            y = ElectrodePositions(10).probmat(:,1);
            
            scatter(x, y, 100, 'ok')
            xlim([-100 200])
            ylim([-100 2500])
        end
        
        scatter(posi_position(:,1)+4*rand-2, -posi_position(:,2)+4*rand-2, 'r', 'filled')
        hold on
        scatter(neg_position(:,1)+4*rand-2, -neg_position(:,2)+4*rand-2, 'b', 'filled')
        title(sprintf('%s-cNE#%d', exp_site_nedata.exp, jj), 'interpreter', 'none')
        saveas(gcf, fullfile(figfolder, sprintf('%s-cNE%d', exp_site_nedata.exp, jj)))
        close
    end
end

%% plot summation of ne member distribution
probe = 'poly';
nefiles = dir(sprintf('*%s*-20dft.mat', probe));
if contains(probe, 'Tetrode')
    figure
    figuresetup2savepdf(20, 20)
    hold on
    plot(-20*[1 1], [50 600], 'k')
    plot(60*[1 1], [50 600], 'k')
    plot(180*[1 1], [50 600], 'k')
    plot(260*[1 1], [50 600], 'k')
    plot(380*[1 1], [50 600], 'k')
    plot(460*[1 1], [50 600], 'k')
    plot(580*[1 1], [50 600], 'k')
    plot(660*[1 1], [50 600], 'k')
    
    plot([-20, 20], [50 0], 'k')
    plot([20, 60], [0 50], 'k')
    plot([180, 220], [50 0], 'k')
    plot([220, 260], [0 50], 'k')
    plot([380, 420], [50 0], 'k')
    plot([420, 460], [0 50], 'k')
    plot([580, 620], [50 0], 'k')
    plot([620, 660], [0 50], 'k')
    
    x = ElectrodePositions(8).probmat(:,5);
    y = ElectrodePositions(8).probmat(:,1);
    
    scatter(x, y, 100, 'ok')
    xlim([-100 700])
    ylim([-100 700])
else
    figure
            figuresetup2savepdf(10, 30)
            hold on
            plot(-50*[1 1], [50 2500], 'k')
            plot(150*[1 1], [50 2500], 'k')
            
            plot([-50, 50], [50 0], 'k')
            plot([50, 150], [0 50], 'k')
            
            x = ElectrodePositions(10).probmat(:,5);
            y = ElectrodePositions(10).probmat(:,1);
            
            scatter(x, y, 100, 'ok')
            xlim([-100 200])
            ylim([-100 2500])
end
%
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    cellPosition = nedata.position;
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    Pattern = nedata.Patterns;
    
    for jj = 1:nNE
        member = NEmembers{jj};
        pattern = Pattern(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        if isempty(neg_idx) || length(posi_idx) <=1
            continue
        end
        
        posi_position = cell2mat(cellPosition(posi_idx)');
        neg_position = cell2mat(cellPosition(neg_idx)');
        load(fullfile('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis', 'ElectrodePositions.mat'))
        
        scatter(posi_position(:,1)+15*rand-8, -posi_position(:,2)+15*rand-8, 20,'r', 'filled')
        hold on
        scatter(neg_position(:,1)+15*rand-8, -neg_position(:,2)+15*rand-8, 20, 'b', 'filled')
    end
end
%% distribution of distance
celldist_posi = cell(1, length(nefiles));
celldist_posi_neg = cell(1, length(nefiles));

for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    cellPosition = nedata.position;
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    Pattern = nedata.Patterns;
    
    celldist_posi_tmp = cell(1, nNE);
    celldist_posi_neg_tmp = cell(1, nNE);
    for jj = 1:nNE
        member = NEmembers{jj};
        pattern = Pattern(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        if isempty(neg_idx) || length(posi_idx) <=1
            continue
        end
        cmb_posi = nchoosek(posi_idx, 2);
        celldist_posi_tmp{jj} = zeros(1, length(cmb_posi));
        celldist_posi_neg_tmp{jj} = zeros(length(posi_idx), length(neg_idx));
        for kk = 1:size(cmb_posi,1)
            celldist_posi_tmp{jj}(kk) = sqrt(sum((cellPosition{cmb_posi(kk,1)} - cellPosition{cmb_posi(kk,2)}).^2));
        end
        for kk = 1:length(neg_idx)
            for ll = 1:length(posi_idx)
                celldist_posi_neg_tmp{jj}(ll, kk) = sqrt(sum((cellPosition{neg_idx(kk)} - cellPosition{posi_idx(ll)}).^2));
            end
        end
        celldist_posi_neg_tmp{jj} = celldist_posi_neg_tmp{jj}(:);
        
    end
    celldist_posi{ii} = cell2mat(celldist_posi_tmp);
    celldist_posi_neg{ii} = cell2mat(celldist_posi_neg_tmp');
end

celldist_posi = cell2mat(celldist_posi);
celldist_posi_neg = cell2mat(celldist_posi_neg');
%%
figure
histogram(celldist_posi, 0:100:2000, 'normalization', 'probability', 'FaceColor', 'r')
%histogram(celldist_posi, 0:50:1000, 'normalization', 'probability', 'FaceColor', 'r')
hold on
histogram(celldist_posi_neg, 0:100:2000, 'normalization', 'probability', 'FaceColor', 'b')
%histogram(celldist_posi_neg, 0:50:1000, 'normalization', 'probability', 'FaceColor', 'b')

xlabel('Distance (um)')
ylabel('Proportion')
p = ranksum(celldist_posi, celldist_posi_neg);
title(sprintf('p = %.5f (ranksum)', p))
ylim([0 0.25])
scatter(mean(celldist_posi), 0.24, 'vr')
scatter(median(celldist_posi), 0.24, 'vr', 'filled')
scatter(mean(celldist_posi_neg), 0.23, 'vb')
scatter(median(celldist_posi_neg), 0.23, 'vb', 'filled')
legend({'Distance bwtween positive members', 'Distance bwtween positive and negative members'})


%% distribution of depth(dim = 2)/width (dim = 1)
nefiles = dir('*poly*-20dft.mat');
dim = 2;
celldist_posi = cell(1, length(nefiles));
celldist_posi_neg = cell(1, length(nefiles));

for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    cellPosition = nedata.position;
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    Pattern = nedata.Patterns;
    
    celldist_posi_tmp = cell(1, nNE);
    celldist_posi_neg_tmp = cell(1, nNE);
    for jj = 1:nNE
        member = NEmembers{jj};
        pattern = Pattern(:,jj);
        posi_idx = member(pattern(member) > 0);
        neg_idx = member(pattern(member) < 0);
        if isempty(neg_idx) || length(posi_idx) <=1
            continue
        end
        cmb_posi = nchoosek(posi_idx, 2);
        celldist_posi_tmp{jj} = zeros(1, length(cmb_posi));
        celldist_posi_neg_tmp{jj} = zeros(length(posi_idx), length(neg_idx));
        for kk = 1:size(cmb_posi,1)
            celldist_posi_tmp{jj}(kk) = abs(cellPosition{cmb_posi(kk,1)}(dim) - cellPosition{cmb_posi(kk,2)}(dim));
        end
        for kk = 1:length(neg_idx)
            for ll = 1:length(posi_idx)
                celldist_posi_neg_tmp{jj}(ll, kk) = abs(cellPosition{neg_idx(kk)}(dim) - cellPosition{posi_idx(ll)}(dim));
            end
        end
        celldist_posi_neg_tmp{jj} = celldist_posi_neg_tmp{jj}(:);
        
    end
    celldist_posi{ii} = cell2mat(celldist_posi_tmp);
    celldist_posi_neg{ii} = cell2mat(celldist_posi_neg_tmp');
end

celldist_posi = cell2mat(celldist_posi);
celldist_posi_neg = cell2mat(celldist_posi_neg');
%%
figure
histogram(celldist_posi, 0:100:2000, 'normalization', 'probability', 'FaceColor', 'r')
hold on
histogram(celldist_posi_neg, 0:100:2000, 'normalization', 'probability', 'FaceColor', 'b')
xlabel('Distance (um)')
ylabel('Proportion')
p = ranksum(celldist_posi, celldist_posi_neg);
title(sprintf('p = %.5f (ranksum)', p))
ylim([0 0.25])
scatter(mean(celldist_posi), 0.24, 'vr')
scatter(median(celldist_posi), 0.24, 'vr', 'filled')
scatter(mean(celldist_posi_neg), 0.23, 'vb')
scatter(median(celldist_posi_neg), 0.23, 'vb', 'filled')
legend({'Depth bwtween positive members', 'Depth bwtween positive and negative members'})
%legend({'Width bwtween positive members', 'Width bwtween positive and negative members'})
