%% run zunchado analitico

E = 200e3; v = 0.3;
a = 300; %mm
sigma_y = 250; %Mpa
p = 100; %Mpa
FS = 1;
tapasFlag = 0;

[espesores, D, interferencia, p_interferencia] = RP_2capa_v2(tapasFlag,p, sigma_y, FS, a, E, v)
b = a + espesores(1);

%% desplazamientos

u = @(r,a,b,pint,pout) (1-v)/E*(a.^2*pint-b.^2*pout)*r/(b.^2-a.^2) + (1+v)/E*(pint-pout)*a.^2*b.^2/(r*(b.^2-a.^2));

% r = 300mm interior
uRa = u(a,a,b,p,p_interferencia)
