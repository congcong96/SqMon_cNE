function NEtrain = ne_upsample_NEact_using_member_neuron_activity(NEact, NEmembers, spktrain, NEthresh)

narginchk(4,4)

% nedata = exp_site_nedata.nedata;
% NEact = nedata.Activities;
% NEthresh = nedata.NEthresh;
% spktrain = nedata.sta_spktrain;
% NEmembers = nedata.NEmembers;
% 
% if isfield(nedata, 'boundaryidx')
%     NEact = NEact(:,1:nedata.boundaryidx-1);
% end

% stadf = 10;
% df = exp_site_nedata.df;
% us_factor = df/stadf;

us_factor = round(size(spktrain,2)/size(NEact,2));

NEtrain = zeros(length(NEmembers), size(NEact,2)* us_factor);

% check for length differences in interpolation
lendiff = size(NEact,2) * us_factor - size(spktrain,2);

if lendiff > 0
    error('Something''s wrong!')
end

for i = 1:size(NEact,1)
        
    eventidx = find(NEact(i,:) > NEthresh(i));
    temp_membertrain = sum(spktrain(NEmembers{i},:), 1);
    
    us_NEtrain = zeros(us_factor, size(NEact, 2));
    membertrain = us_NEtrain;
    
    for j = 1:us_factor
        try
            membertrain(j,:) = temp_membertrain(j:us_factor:end);
        catch % if membertrain lacks one element due to non-divisible length
            membertrain(j,:) = temp_membertrain(j:us_factor:end+lendiff);
        end
    end
    
    for j = 1:length(eventidx)
        
        [maxval, maxidx] = max(membertrain(:,eventidx(j)));
        nummax = sum(maxval == membertrain(:,eventidx(j)));
        
        if nummax > 1
            %maxidx = randi(us_factor,1);
            maxidx = us_factor;
        end  
        us_NEtrain(maxidx, eventidx(j)) = 1;
    end
    
     NEtrain(i,:) = us_NEtrain(:)';
    
end



end

