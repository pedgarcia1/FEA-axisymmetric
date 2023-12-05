

% Preasignación de celdas para almacenar los fotogramas del GIF
frames = cell(1, 10);
aux = linspace(1, 200, 10);
% Generar fotogramas para el GIF
for idx = 1:10
    magnificationFactor = aux(idx);
    
    % Deformed plot
    fig = figure('visible', 'on'); % Para evitar que se muestren todas las figuras
    % title(sprintf('Deformed Plot MF: %d', magnificationFactor));
    set(fig, 'Color', 'white');
    
    % Ajustar el tamaño de la figura para mejorar la calidad
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 800, 600]);
    
    % Puedes ajustar el colormap si es posible con la función meshPlot
    meshPlot(elements, nodes + magnificationFactor * reshape(displacementsVector, nDimensions, nNodes)', 'b', 'No');
    
    frames{idx} = getframe(gcf);
    close(fig); % Cerrar la figura después de capturar el cuadro
end

% Guardar los cuadros como un archivo de GIF
filename = 'deformed_plot_high_quality.gif';
for idx = 1:numel(frames)
    [A, map] = rgb2ind(frames{idx}.cdata, 256);
    if idx == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end
