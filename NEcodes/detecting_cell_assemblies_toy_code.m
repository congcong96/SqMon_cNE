function detecting_cell_assemblies_toy_code(nbins,nacts,nneurons, assemblies)
% detecting_cell_assemblies_toy_code Example model cell assembly detection code
% 
%     detecting_cell_assemblies_toy_code(nbins,nacts,assemblies)
%     nbins : time bins in spike trains
%     nacts : number of times each assembly fires
%     assemblies : cell array of vectors, where each vector is a cell assembly.
%     Default:     Assembly_opts.assembly_neurons{1} = [1 2 3 4];
%                  Assembly_opts.assembly_neurons{2} = [5 6 7];
% 
%     The simulation has a default of 32 neurons, 10000 bins, and 500 activations.
%     The mean spike rate is 1 spike/bin, and the spike rate for an activation is 3 spikes/bin
%
%     Parameters that affect how many assemblies can be detected:
%     1. Number of time bins
%     2. Number of activations
%     3. Number of total neurons
%     4. Do the assemblies share neurons?
%
%




if ( nargin == 0 )
    Network_opts.nneurons = 32;
    Network_opts.nbins = 10000;
    Assembly_opts.number_of_activations = 500;
    Assembly_opts.assembly_neurons{1} = [1 2 3 4];
    Assembly_opts.assembly_neurons{2} = [5 6 7];
end

if ( nargin == 1 )
    if isempty(nbins)
        Network_opts.nbins = 10000;
    else
        Network_opts.nbins = nbins;
    end
    Network_opts.nneurons = 32;
    Assembly_opts.number_of_activations = 500;
    Assembly_opts.assembly_neurons{1} = [1 2 3 4];
    Assembly_opts.assembly_neurons{2} = [5 6 7];
end

if ( nargin == 2 )
    if isempty(nbins)
        Network_opts.nbins = 10000;
    else
        Network_opts.nbins = nbins;
    end

    if isempty(nacts)
        Assembly_opts.number_of_activations = 500;
    else
        Assembly_opts.number_of_activations = nacts;
    end

    Network_opts.nneurons = 32;
    Assembly_opts.assembly_neurons{1} = [1 2 3 4];
    Assembly_opts.assembly_neurons{2} = [5 6 7];
end

if ( nargin == 3 )

    if isempty(nbins)
        Network_opts.nbins = 10000;
    else
        Network_opts.nbins = nbins;
    end

    if isempty(nacts)
        Assembly_opts.number_of_activations = 500;
    else
        Assembly_opts.number_of_activations = nacts;
    end

    if isempty(nneurons)
        Network_opts.nneurons = 32;
    else
        Network_opts.nneurons = nneurons;
    end

    Assembly_opts.assembly_neurons{1} = [1 2 3 4];
    Assembly_opts.assembly_neurons{2} = [5 6 7];
end

if ( nargin == 4 )

    if isempty(nbins)
        Network_opts.nbins = 10000;
    else
        Network_opts.nbins = nbins;
    end

    if isempty(nacts)
        Assembly_opts.number_of_activations = 500;
    else
        Assembly_opts.number_of_activations = nacts;
    end

    if isempty(nneurons)
        Network_opts.nneurons = 32;
    else
        Network_opts.nneurons = nneurons;
    end

    if isempty(assemblies)
        Assembly_opts.assembly_neurons{1} = [1 2 3 4];
        Assembly_opts.assembly_neurons{2} = [5 6 7];
    else
        Assembly_opts.assembly_neurons = assemblies;
    end

end

fprintf('\n');
fprintf('Nbins = %.0f\n', Network_opts.nbins);
fprintf('Nacts = %.0f\n', Assembly_opts.number_of_activations);
fprintf('Nneurons = %.0f\n', Network_opts.nneurons);
fprintf('\n');




Network_opts.meanspikebin = 1;
Assembly_opts.meanspikerate_activations = 3;



Activitymatrix = toy_simulation(Network_opts, Assembly_opts);
correlationmat = corr(Activitymatrix');
Patterns = assembly_patterns(Activitymatrix);
Activities = assembly_activity(Patterns,Activitymatrix);



figure;

subplot(2,3,1);
imagesc(correlationmat);
xlabel('Neuron #');
ylabel('Neuron #');
tickpref;


subplot(2,3,4);
imagesc(Patterns);
xlabel('Assembly #');
ylabel('Neuron #');
tickpref;
set(gca,'xtick', 1:size(Patterns,2), 'xticklabel', 1:size(Patterns,2));


subplot(2,3,[2 3]); 
imagesc(Activitymatrix);
xlim([0 round(Network_opts.nbins/20)]);
xlim([0 500]);
ylabel('Neuron #');
tickpref;

subplot(2,3,[5 6]);
plot(Activities');
xlim([0 round(Network_opts.nbins/20)]);
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
set(gcf,'position', [496 558 744 420]);


return;



