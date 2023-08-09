%% run zunchado analitico

E = 200e3; v = 0.3;
a = 300; %mm
sigma_y = 250; %Mpa
p = 100; %Mpa
FS = 1;

[espesores, D, interferencia, p_i]= RP_2capa_v2(p, sigma_y, FS, a, E, v)