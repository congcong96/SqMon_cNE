cd('/data/congcong/SqMoPhys_Josh/cNE_analysis')
nefiles = dir('*20dft.mat');
ICweights = cell(size(nefiles));
ICweights_proportion = cell(size(nefiles));

for ii = 1:length(nefiles)
    load(nefiles(ii).name, 'exp_site_nedata')
    Patterns = exp_site_nedata.nedata.Patterns;
    Patterns_proportion = zeros(size(Patterns));
    for jj = 1:size(Patterns,2)
        pattern = Patterns(:,jj);
        [~,idx] = max(abs(pattern));
        if pattern(idx) < 0
            pattern = - pattern;
            Patterns(:,jj) = pattern;
        end
        pattern = pattern/max(pattern);
        Patterns_proportion(:,jj) = pattern;
    end
    ICweights{ii} = Patterns;
    ICweights_proportion{ii} = Patterns_proportion;
    exp_site_nedata.nedata.Patterns = Patterns;
    save(nefiles(ii).name, 'exp_site_nedata')
end
ICweights_origin = cell2mat(cellfun(@(x) x(:), ICweights, 'UniformOutput', false));
ICweights_origin = sort(ICweights_origin, 'descend')/max(ICweights_origin);
ICweights_proportion = cell2mat(cellfun(@(x) x(:), ICweights_proportion, 'UniformOutput', false));
ICweights_proportion = sort(ICweights_proportion, 'descend')/max(ICweights_proportion);
%%
subplot(221)
stem(ICweights_origin, 'marker', '.', 'color', [0.2, 0.4 1.0])
xlim([1 length(ICweights_origin)])
ylim([-1 1])
xlabel('IC #')
ylabel('IC Weight')
% curvefit = fit(x,ICweights_origin,'poly9','Normalize','on');
% hold on
% plot(curvefit,x,ICweights_origin, [0.2, 0.4 1.0], 'linewidth', 2)
subplot(222)
x = (1:length(ICweights_origin))';
ICwcdf = x/length(ICweights_origin);
plot(flip(ICweights_origin),ICwcdf, 'color', [0.2, 0.4 1.0], 'linewidth', 2)
title('ICweight CDF')
xlabel('ICweight')
ylabel('cumulative proportion')
% subplot(223)
% [deriv1, deriv2] = differentiate(curvefit, x);
% plot(x, deriv1, 'color', [0.2, 0.4 1.0], 'linewidth', 2)
% xlim([1 length(x)])
% ylabel('1st derivative of ICweight CDF')
% xlabel('ICweight')
% subplot(224)
% plot(x, deriv2, 'color', [0.2, 0.4 1.0], 'linewidth', 2)
% xlim([1 length(x)])
% ylabel('2nd derivative of ICweight CDF')
% xlabel('ICweight')

%%
subplot(223)
%ICweights_proportion2 = ICweights_proportion(ICweights_proportion < 1);
x = (1:length(ICweights_proportion))';
%curvefit = fit(x,ICweights_proportion2,'poly9','Normalize','on');
stem(ICweights_proportion, 'marker', '.', 'color', [1.0, 0.4, 0.2])
hold on
plot(x, 0.37*ones(size(ICweights_proportion)), 'k--')
plot(x, -0.37*ones(size(ICweights_proportion)), 'k--')
%plot(curvefit)
xlim([1 length(ICweights_proportion)])
ylim([-1 1])
xlabel('IC #')
ylabel('relative IC Weight')
subplot(224)
plot(flip(ICweights_proportion), (1:length(ICweights_proportion))/length(ICweights_proportion),  'color', [1.0, 0.4, 0.2], 'linewidth', 2)
title('relative ICweight CDF')
xlabel('relative ICweight')
ylabel('cumulative proportion')
hold on
y = ICwcdf(find(flip(ICweights_proportion) > -0.37, 1));
plot([-1 -0.37], y*[1 1], 'k--')
plot([-0.37 -0.37], y*[0 1], 'k--')
y = ICwcdf(find(flip(ICweights_proportion) > 0.37, 1));
plot([-1 0.37], y*[1 1], 'k--')
plot([0.37 0.37], y*[0 1], 'k--')
% subplot(223)
% [deriv1, deriv2] = differentiate(curvefit, x);
% plot(x, deriv1, 'color', [0.2, 0.4 1.0], 'linewidth', 2)
% xlim([1 length(x)])
% ylabel('1st derivative of ICweight CDF')
% xlabel('ICweight')
% subplot(224)
% plot(x, deriv2, 'color', [0.2, 0.4 1.0], 'linewidth', 2)
% xlim([1 length(x)])
% ylabel('2nd derivative of ICweight CDF')
% xlabel('ICweight')

%% PLOT correlation matrix, IC weights and ICweigths expansion
for ii = 1:length(nefiles)
    nfig = 1;
    load(nefiles(ii).name, 'exp_site_nedata')
    Patterns = exp_site_nedata.nedata.Patterns;
    corrmat = corr(exp_site_nedata.nedata.spktrain');
    corrmat = corrmat - eye(size(corrmat));
    for jj = 1:size(Patterns,2)
        nplot = mod(jj, 4);
        if nplot == 1
            figure
            
        elseif nplot == 0
            nplot = 4;
        end
        pattern = Patterns(:,jj);
        [~, idx] = sort(pattern, 'descend');
        
    end
end