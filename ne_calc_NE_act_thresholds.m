function [thresh, alpha] = ne_calc_NE_act_thresholds(exp_site_nedata,perm,num_iter,alpha,varargin)

% varargin: Enter shuffled IC weights here if required.

narginchk(1,5)

if nargin == 1
    perm = 'circular';
    num_iter = 20;
    alpha = 99.9;
elseif nargin == 2
    num_iter = 20;
    alpha = 99.9;
elseif nargin == 3
    alpha = 99.9;
end

if (strcmp(perm,'random') | strcmp(perm,'circular')) == 0
    error('Permutation option should either be "random" or "circular".')
end

if isfield(exp_site_nedata, 'nedata')
    nedata = exp_site_nedata.nedata;
else
    nedata = exp_site_nedata;
end
spktrain = nedata.spktrain;


% Use alternative IC weights if input demands
if isempty(varargin)
    try
        ICwt = nedata.Patterns;
    catch
        ICwt = nedata.total_patterns;
    end
else
    ICwt = varargin{1};
end

% thresh = zeros(size(ICwt,2),1);
bigpermmat = zeros([size(ICwt,2),size(spktrain,2),num_iter]);
[M, N] = size(spktrain);


for i = 1:num_iter
    
    fprintf('%d of %d iterations...\n',i,num_iter);
    temp = zeros(M,N);
    
    switch perm
        case 'random'
   
            for j = 1:M
                X = rand(N,1);
                [~,idx] = sort(X);
                temp(j,:) = spktrain(j,idx);
            end

        case 'circular'
        
            shiftidx = random('unid',N,[M 1]);
   
            for j = 1:M
                temp(j,:) = circshift(spktrain(j,:),[0 shiftidx(j)]);               
            end
    end
    
    bigpermmat(:,:,i) = assembly_activity(ICwt,temp);
    
    
end

bigpermmat = reshape(bigpermmat,size(bigpermmat,1),N*num_iter);


thresh = zeros(length(alpha), size(ICwt,2));
% Jermyn 2018 paper
for i = 1:length(alpha)
    thresh(i,:) = prctile(bigpermmat,alpha(i),2);
end

return


