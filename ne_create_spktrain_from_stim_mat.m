function [spktrain_mat, edges, varargout] = ne_create_spktrain_from_stim_mat(spk, stimulus, trigger)

if isstruct(spk) && length(spk) ~= 1
    
    spktrain_mat = zeros(length(spk), size(stimulus,2));
    position = cell(length(spk),1);
    FsAD = spk(1).fs; % A/D system sampling rate
    spktimes = {spk.spiketimes};
    
elseif isstruct(spk) && length(spk) == 1
    
    spktrain_mat = zeros(length(spk.spk), size(stimulus,2));
    position = cell(length(spk.spk), 1);
    FsAD = spk.fs;
    
    spk = spk.spk;
    spktimes = {spk.spiketimes};
    
else
    spktrain_mat = zeros(size(spk,1), size(stimulus,2));
    spktimes = spk;
    FsAD = 20000;
end
        

trig_ms = double(trigger) / double(FsAD) * 1000; % convert trigger to ms
N = size(stimulus,2)/length(trigger);    
    
if mod(N,1) ~= 0
    warning('Wrong number of trigger events')
    spktrain_mat = [];
    edges = [];
    varargout = {[]};
    return
end

max_trig_diff = max(diff(trig_ms));
trig_ms = [trig_ms trig_ms(end) + max_trig_diff];
trig_diff = diff(trig_ms);
trig_step = trig_diff/N;

edges = unique(cell2mat(arrayfun(@(x,y,z) x:z:y, trig_ms(1:end-1),...
    trig_ms(2:end),trig_step,'UniformOutput',0)));

for i = 1:length(spk)
    [spktrain, ~] = histcounts(spktimes{i}, edges);
    spktrain(1:N) = 0;

    spktrain_mat(i,:) = spktrain;
%     if isstruct(spk)
%         position{i} = num2str(spk(i).position);
%     end
end % (for i)

if exist('position','var') && nargout == 3
    varargout{1} = position;
end

end

