function [espesores, D, interferencia, p_i]= RP_2capa_v2(p, sigma_y, FS, a, E, v)

e = 0.1:1:250; 
%Fijo estos espesores. El caso de una capa me dio 49,13, obligatoriamente 
%voy a estar abajo de eso se supone, sino no tendría sentido la estructura esta

i = 1;
espesor_1 = zeros(100, length(e));
espesor_2 = zeros(100, length(e));
int = zeros(100, length(e));
%Estas matrices guardan los espesores de de cada tubo, para un caso de
%interferencia (filas) y de b (columnas). Si con esos parámetros, no
%encuentra nada que funcione, usa un -1.

p_int_vec = 0.1:0.1:50; %Mpa
for p_int = p_int_vec
    %Voy a ir cubriendo distintas interferencias
    k = 1;
    for b = a+e %Recorro los distintos valores de b posibles
        c = b+e;
        
%         sigma_z = p*a^2./(c.^2-a^2); 
        sigma_z1 = p*a^2./(b.^2-a^2); 
        sigma_z2 = p*a^2./(c.^2-b^2); 
        
        %Las tensiones se analizan por superposición. El 1 es el cilindro
        %interno y el 2 el externo. Cada uno se ve afectado por la presión
        %de interferencia de su manera correspondiente, y también como un
        %continuo por la presión interna. Se asume un b igual para ambos en
        %la interfase. Se evalúa en el radio interno de cada tubo
        
        %p_int = E*int*(b^2 - a^2)*(c.^2-b^2)./((2-v)*(b^3.*(c.^2-a^2)));
        sigma_r1 = -p;
%         sigma_tita1 = (p./(c.^2-a^2)).*(a^2 + c.^2) + (p_int/(b^2-a^2))*(-2*b^2);
        sigma_tita1 = (p./(b.^2-a^2)).*(a.^2 + b.^2) + (p_int/(b.^2-a.^2))*(-2*b.^2);

        
%         sigma_r2 = (p./(c.^2-a^2)).*(a^2 - (a^2*c.^2)/b^2) - p_int;
        sigma_r2 = -p_int;
%         sigma_tita2 = (p./(c.^2-a^2)).*(a^2 + (a^2*c.^2)/b^2)+(p_int./(c.^2-b^2)).*(b^2+c.^2);
%         sigma_tita2 = (p./(c.^2-b.^2)).*(b.^2 + (a.^2*c.^2)/b^2)+(p_int./(c.^2-b^2)).*(b.^2+c.^2);
        sigma_tita2 = p_int.*(b.^2 + c.^2)./(c.^2 - b.^2);


        %Tresca para ambos tubos
%         sigma1 = max([abs(sigma_z-sigma_r1); abs(sigma_z-sigma_tita1); abs(sigma_tita1-sigma_r1)]);
%         sigma2 = max([abs(sigma_z-sigma_r2); abs(sigma_z-sigma_tita2); abs(sigma_tita2-sigma_r2)]);
                
        % vm para ambos tubos
        sigma1 = sqrt(0.5*((sigma_r1-sigma_tita1).^2 + (sigma_tita1-sigma_z1).^2 + (sigma_z1-sigma_r1).^2 ));
        sigma2 = sqrt(0.5*((sigma_r2-sigma_tita2).^2 + (sigma_tita2-sigma_z2).^2 + (sigma_z2-sigma_r2).^2 ));
        
        %Condiciones para que un valor de c sea válido:
        % -Presión de interferencia < 10 MPa
        % -No pasar la fluencia en ninguno de los dos tubos
%         cond_int = p_int < 10;
        cond_int = p_int > 1;
        cond_y1 = sigma1 < sigma_y/FS;
        cond_y2 = sigma2 < sigma_y/FS;
        
        %Elijo el valor de c más chico que cumpla las condiciones, para un
        %par (interferencia, b) en particular
        c_min = min(c(cond_int & cond_y1 & cond_y2));
        
        %Lleno las matrices
        if ~isempty(c_min)
            espesor_2(i,k) = c_min - b;
            espesor_1(i,k) = b - a;
            int(i,k) = p_int/(E*(b^2 - a^2)*(c_min.^2-b^2)/((2-v)*(b^3*(c_min^2-a^2))));
        else
            espesor_2(i,k) = -1;
            espesor_1(i,k) = -1;
            int(i,k) = -1;
        end
        k = k+1;
    end
    i = i+1;
end

%Busco el elemento cuya suma de espesores sea mínima y positiva
sum_espesores = espesor_2 + espesor_1;
[x, y] = find(sum_espesores == min(sum_espesores(sum_espesores>0)), 1); 

espesores = [espesor_1(x,y), espesor_2(x,y)];
D = (a + min(sum_espesores(sum_espesores>0)))*2;
interferencia = int(x,y);
p_int = p_int_vec;
p_i = p_int(x);
end