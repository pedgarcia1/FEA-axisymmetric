function plotColo(nodes,elements,scalarField)
hold on
% figure
% for iElement=1:size(elements,1)
%     patch('Faces',f(iElement,:),'Vertices',v(f(iElement,:),:),'FaceVertexCData',col(iElement,1),'FaceColor','interp');
% end

for i = 1:size(scalarField, 1)
    patch('Faces', [1 2 3 4], 'Vertices', nodes(elements(i,:)',:), 'FaceVertexCData', scalarField(i,:)', ...
          'FaceColor', 'interp', 'EdgeColor', 'k');
    hold on;
end
colormap(jet(10))
c = colorbar;
c.Label.String = 'Tension [MPa]';
%set(gca,'TickLabelInterpreter','latex','Limits',[min(cg) max(cg)],'box','off','FontSize',25)
%c.TickLabelInterpreter = 'latex';
% c.FontSize = 8;
%c.Limits = [min(cg) max(cg)];
% c.Box = 'off';
% c.Location = 'southoutside';
maxCbar = max(max(scalarField(elements)));
minCbar = min(min(scalarField(elements)));
caxis([minCbar maxCbar])
colormap parula

axis square
set(gca,'XTick',[],'YTick',[])
end