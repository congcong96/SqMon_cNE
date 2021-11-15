cd('/data/congcong/SqMoPhys_Josh/cNE_analysis')
nefiles = dir('*20dft.mat');
Npos = cell(size(nefiles));
Nneg = cell(size(nefiles));
file_number = cell(size(nefiles));
NEnumber = cell(size(nefiles));
for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    Patterns = nedata.Patterns;
    nNE = length(nedata.NEmembers);
    Npos_tmp = zeros(1, nNE);
    Nneg_tmp = zeros(1, nNE);
    
    for jj  = 1:nNE
        pattern = Patterns(:,jj);
        MemberWeights = pattern(nedata.NEmembers{jj});
        Npos_tmp(jj) = sum(MemberWeights > 0);
        Nneg_tmp(jj) = sum(MemberWeights < 0);
    end
    Npos{ii} = Npos_tmp;
    Nneg{ii} = Nneg_tmp;
    file_number{ii} = ones(1, nNE)*ii;
    NEnumber{ii} = 1:nNE;
end

Nneg = cell2mat(Nneg');
Npos = cell2mat(Npos');
file_number = cell2mat(file_number');
NEnumber = cell2mat(NEnumber');
idx = find(Npos== 1);
Nneg(idx) = [];
Npos(idx) = [];
file_number(idx) = [];
NEnumber(idx) = [];
% histogram - distribution of NEs with different number of negative members
[c, edges] = histcounts(Nneg, 'Normalization', 'probability');
centers = (edges(1:end-1) + edges(2:end))/2;
bar(centers, c)
xlabel('Number of negative members')
ylabel('Proportion')
% histogram - proportion of negative members
figure
NegProportion = Nneg./Npos;
NegProportion(NegProportion == 0) = [];
edges = [0:0.1:1 10];
c = histcounts(NegProportion, edges);
bar(0.05:0.1:0.1*length(c), c, 'barwidth', 1)
xticks([0.1:0.2:0.9 1.05])
xticklabels({'0.1' '0.3' '0.5' '0.7' '0.9' '>1'})
xlabel('# negative members / # positive members')
ylabel('Number of cNE')
%%
% histogram - proportion of negative members
NegProportion = Nneg(Nneg > 1)./Npos(Nneg > 1);
edges = [0:0.1:1 10];
c = histcounts(NegProportion, edges);
bar(0.05:0.1:0.1*length(c), c, 'barwidth', 1)
xticks([0.1:0.2:0.9 1.05])
xticklabels({'0.1' '0.3' '0.5' '0.7' '0.9' '>1'})
xlabel('# negative members / # positive members')
ylabel('Number of cNE')