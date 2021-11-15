function [Patterns, Activities] = detecting_cell_assemblies_data(Activitymatrix, position)
% detecting_cell_assemblies_data Cell assembly detection using PCA/ICA
% 
%     detecting_cell_assemblies_toy_code(Activitymatrix) estimates the cell
%     assemblies from the spike train matrix Activitymatrix.
%
%     Activitymatrix : spike train matrix. m x n matrix. m = number of neurons, n = number of
%     time bins. Activitymatrix(i,j) = # of spikes by neuron i at time bin
%     j. 
% 
%     Output shows the spike train matrix, the correlation between the
%     different neurons, the cell assemblies, and the activation patterns
%     of each assembly. The activation pattern is the time varying activity
%     of the assembly.
%
%
%    detecting_cell_assemblies_toy_code(Activitymatrix, position) places
%    the electrode recording channel depth on the plots instead of the
%    neuron number. This may help to determine how neurons are synchronized
%    across the cortical column.







% Activitymatrix = toy_simulation(Network_opts, Assembly_opts);

correlationmat = corr(Activitymatrix');
Patterns = assembly_patterns(Activitymatrix);
Activities = assembly_activity(Patterns, Activitymatrix);


figure;

subplot(2,3,1);
imagesc(correlationmat);
if ( nargin == 1 )
    xlabel('Neuron #');
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'xtick', tick, 'xticklabel', position);
    set(gca,'ytick', tick, 'yticklabel', position);
    xlabel('Position (um)');
    ylabel('Position (um)');
%     rotateticklabel(gca);
end
tickpref;


subplot(2,3,[2 3]); 
zSpikeCount = zscore(Activitymatrix');
imagesc(zSpikeCount');
% imagesc(Activitymatrix);
xlim([0 round(size(Activitymatrix,2)/20)]);
xlim([0 500]);
if ( nargin == 1 )
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
tickpref;




subplot(2,3,4);
imagesc(Patterns);
xlabel('Assembly #');
if ( nargin == 1 )
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
tickpref;
set(gca,'xtick', 1:size(Patterns,2), 'xticklabel', 1:size(Patterns,2));


subplot(2,3,[5 6]);
plot(Activities');
xlim([0 round(size(Activitymatrix,2)/20)]);
xlim([0 500]);
xlabel('Time bin');
tickpref;
[nr, nc] = size(Activities');
if ( nc > 0 )
    for i = 1:nc
        leg{i} = num2str(i);
    end
    legend(leg,'Location','Best');
end
set(gcf,'position', [432 334 1121 620]);




% for i = 1:size(Patterns,2)
%     figure;
%     stem(Patterns(:,i));
%     title(sprintf('Ind Comp #%.0f', i));
% end




return;



