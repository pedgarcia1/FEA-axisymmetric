%% run zunchado analitico

E = 200e3; v = 0.3;
a = 300; %mm, radio interior
sigma_y = 250; %Mpa
p = 100; %Mpa
FS = 1;
tapasFlag = 0; % tapas serian los casquetes, si no hay la tension en z es distinta
planeStrainFlag = 0;
% para solo zunchado
% tapasFlag = 0; planeStrainFlag = 0;

if planeStrainFlag
    E = E/(1-nu^2);
    nu = nu/(1-nu);
    fprintf("Caso plane strain \n")
else
    fprintf("Caso plane stress \n")
end

[espesores, D, ~, p_interferencia] = RP_2capa_v2(tapasFlag,planeStrainFlag,p, sigma_y, FS, a, E, v)
b = a + espesores(1);
c = b + espesores(2);

%% desplazamientos

u = @(r,a,b,pint,pout) (1-v)/E*(a.^2*pint-b.^2*pout)*r/(b.^2-a.^2) + (1+v)/E*(pint-pout)*a.^2*b.^2/(r*(b.^2-a.^2));

% r = 300mm interior
uRa = u(a,a,b,p,p_interferencia)

%% ESFERA
% thick walled solution
% chequeado del saad elasticity pag. 388
Ri = 300;
% t = 100;
Ro = @(t) Ri + t;
r = 300;
sigmarr = @(t) p*Ri^3/(Ro(t)^3 - Ri^3)*(1 - Ro(t)^3/r^3);
sigmatt = @(t) p*Ri^3/(Ro(t)^3 - Ri^3)*(1 + 1/2*Ro(t)^3/r^3); % = sigma en la otra cord
svm = @(t) sigma_y - 1/sqrt(2)*sqrt((sigmarr(t)-sigmatt(t))^2 + (sigmatt(t)-sigmarr(t))^2);

options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
tEsfTeorico = fsolve(svm,10,options)

%% Sin Zunchado, con Casquetes (tapas)

syms b 
p_int = 0;
sigma_r1 = -p;
sigma_tita1 = (p./(b.^2-a^2)).*(a.^2 + b.^2) + (p_int/(b.^2-a.^2))*(-2*b.^2);
sigma_z1 = p*a^2./(b.^2-a^2);
sigma1 = sqrt(0.5*((sigma_r1-sigma_tita1).^2 + (sigma_tita1-sigma_z1).^2 + (sigma_z1-sigma_r1).^2 ));
solve(sigma1 == 250,b)

