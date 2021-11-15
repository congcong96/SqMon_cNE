%% from the website
% http://stackoverflow.com/questions/10391849/printing-figure-at-same-size-as-figure-window

function [filename] =  printPDFandPSC(fig, fileName, pdfpsc) %printpdf(fig, name)
% printpdf Prints image in PDF format without tons of white space

% The width and height of the figure are found
% The paper is set to be the same width and height as the figure
% The figure's bottom left corner is lined up with
% the paper's bottom left corner

if nargin < 3
    pdfpsc = false; % Natsumi added these lines to only save pdf as a default
end

% Set figure and paper to use the same unit
%set(fig,'Units','inches');
 set(fig, 'Units', 'centimeters');
 set(fig, 'PaperUnits','centimeters');

% Position of figure is of form [left bottom width height]
% We only care about width and height
fig_pos = get(fig,'position');
%pos = get(fig,'Position');

% Set paper size to be same as figure size
%set(fig, 'PaperSize', [pos(3) pos(4)]);

% Set figure to start at bottom left of paper
% This ensures that figure and paper will match up in size
set(fig,'PaperPositionMode','auto'); 
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
% set(fig, 'PaperPositionMode', 'manual');
%set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);

% Set these to get A4 size?
format long
format compact
set(0, 'defaultFigurePaperType', 'A4')
set(0, 'defaultFigurePaperUnits', 'centimeters')
set(0, 'defaultFigurePaperPositionMode', 'auto')


% Print as pdf
print('-dpdf', '-noui','-painters','-r1200', '-bestfit', fileName);
if pdfpsc == true
    print('-dpsc2', '-noui', '-painters', '-bestfit', fileName);% both work well w/ illustrator
end
%print(fig, '-dpdf', name)

% Return full file name
filename = [fileName, '.pdf/ps'];
%filename = [name, '.pdf'];
end