function scalarbandplot(elementos,nodos,variable,lims,lineColor,nColores,flag)

% BANDPLOT  Graficador de variables.
% 
% BANDPLOT(elementos,nodos,variable)
% BANDPLOT(elementos,nodos,variable,lims)
% BANDPLOT(elementos,nodos,variable,lims,lineColor)
% BANDPLOT(elementos,nodos,variable,lims,lineColor,nColores)
% 
% Numeraci�n de los nodos de los elementos:
%  4---7---3
%  |       |
%  8   9   6
%  |       |
%  1---5---2
% 
% elementos: Matriz de conectividades.
% nodos:     Matriz de coordenadas nodales.
% variable:  Matriz m x 1 con la variable a graficar, donde [m n] = size(nodos).
% lims:      Vector [cmin cmax] con los l�mites a utilizar para la variable.
% lineColor: String que especifica el color para los bordes de los
%            lementos. Utilizar 'none' para no mostrar los bordes.
% nColores:  Cantidad de bandas de colores a utilizar.

error(nargchk(3, 7, nargin))


if (nargin < 4) || isempty(lims)
    lims = [min(min(variable)) max(max(variable))];
end

if (nargin < 5) || isempty(lineColor)
    lineColor = 'k';
end


tol = 1E-5;
nMinCol = 3;
if (nargin < 6) || isempty(nColores)
    if diff(lims) < lims(2)*tol
        nColores = nMinCol;
    else
        nColores = 10;
    end
end

if diff(lims) < lims(2)*tol
    lims = lims + [-1 1]*lims(2)*tol;
end

if nColores < nMinCol
    nColores = nMinCol;
end

nNodos = size(elementos,2);
nel = size(elementos,1);

switch nNodos
    case {8,9}
        vNod = [1 5 2 6 3 7 4 8];
    otherwise
        vNod = 1:nNodos;
end

% Creo los patches
switch flag
    case 'interp'
        h1 = patch('Faces',elementos(:,vNod),'Vertices',nodos);
        set(h1,'EdgeColor',lineColor,'FaceVertexCData',variable,'FaceColor','interp');
    case 'flat'
        h1 = patch('Faces',elementos(:,vNod),'Vertices',nodos);
        set(h1,'EdgeColor',lineColor,'FaceVertexCData',variable,'FaceColor','flat');
end

% Acomodo escalas y dem�s
colormap(jet(nColores))
caxis(lims)

if nColores <= 20
    nTicks = nColores;
else
    nTicks = 20;
end

ticks = lims(1):((diff(lims))/nTicks):lims(2);
tickLabels = cell(size(ticks));

for iTick = 1:length(ticks)
    tickLabels{iTick} = sprintf('%6.5E',ticks(iTick));
end

colorbar('YTick',ticks,'YTickLabel',tickLabels);
set(gca,'XTick',[],'YTick',[])
daspect([1 1 1])



