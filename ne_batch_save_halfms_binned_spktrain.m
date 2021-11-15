function badfiles = ne_batch_save_halfms_binned_spktrain(files, stimfile)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

badfiles = {};

for i = 1:length(files)
    
    fprintf('\nProcessing %s...',files{i})
    
    load(files{i})
    
    if exist('spktrain','var')
        fprintf('\nSpiketrain for %s already calculated...\n', files{i})
        
        clear('spk','trigger','spktrain')
        
        continue
    end
    if ~exist('stimfile', 'var') && ~exist('stim_mat', 'var')
        
        stim = spk.stim;
        stimlength = 10;
        
        stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
        
        stimfilename = sprintf('%s-*%dmin_DFt1_DFf5_matrix.mat', stim, stimlength);
        stimfile = fullfile(stimfolder,stimfilename);
        stimfile = gfn(stimfile,1);
        load(stimfile{i}, 'stim_mat');
    elseif ~exist('stim_mat', 'var')
        load(stimfile, 'stim_mat');
    end
    
    [spktrain, edges] = ne_create_spktrain_from_stim_mat(spk, stim_mat, trigger(1,:));
    
    if isempty(spktrain)
        badfiles = [badfiles files{i}];
        continue
    end
    
    
    fprintf('\nSaving %s...\n', files{i})
    
    save(files{i},'spktrain', 'edges', '-append')
    
    clear('spk','trigger', 'edges')
    
end

fprintf('\n')

