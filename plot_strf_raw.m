function plot_strf_raw(filter, faxis, taxis)

if size(filter, 1) == 1
    filter = reshape(filter, [length(faxis), length(taxis)]);
    filter = fliplr(filter);
end
fticks = [];
fticklabels = [];

for ii = -1:1:5
    fticks = [fticks find(faxis/100 >= 2^ii*10, 1)];
    fticklabels = [fticklabels 2^ii];
end

[fticks, ia] = unique(fticks);
fticklabels = fticklabels(ia);

tticks = [];
tticklabels = [];
for ii = 0:25:100%-1:2:5
    tticks = [tticks find(taxis*1000 >= ii, 1)];
    tticklabels = [tticklabels ii];
end
[tticks, ia] = unique(tticks);
tticklabels = tticklabels(ia);

imagesc(filter)
cmap = cbrewer('div', 'RdBu', 11);
colormap(gca, flipud(cmap))
clim = max(abs(filter(:)));
caxis([-clim clim]);
axis xy

xticks(tticks)
xticklabels(tticklabels)
xlim([tticks(1), tticks(end)])
yticks(fticks)
yticklabels( fticklabels)

xlabel('ms')
ylabel('kHz')