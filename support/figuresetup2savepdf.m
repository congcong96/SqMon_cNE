%% function figuresetup2savepdf(width, height)
% this set up figures for saveing as a pdf with the specified size (width x height in cm).
% this should meet the requrements for publication.
% Natsumi

function figuresetup2savepdf(width, height)
% this set up figure for the best to save a
%set(gcf,'paperpositionmode','manual')
set(gcf, 'Renderer','opengl','Units', 'centimeters' , 'Position', [0 0 width height],'PaperUnits', 'centimeters' , 'PaperSize', [height width])
set(gcf,'color','w');
set(gca,'FontSize',9)
set(gca,'FontName','Arial')
end