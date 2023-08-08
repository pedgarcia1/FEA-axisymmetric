function [kG] = K_tincho_Q8(elementos,nodos,Caxi)
syms ksi eta
X = [1 ksi eta ksi^2 ksi*eta eta^2 ksi^2*eta eta^2*ksi];
fpos = @(ksi,eta) [1 ksi eta ksi^2 ksi*eta eta^2 ksi^2*eta eta^2*ksi];
A = [fpos(-1,-1) 
     fpos(1,-1) % 2
     fpos(1,1) % 3
     fpos(-1,1)  %4 
     fpos(0,-1) %5
     fpos(1,0) %6
     fpos(0,1) %7
     fpos(-1,0)]; %Nodo 8 
 shapefuns = X/A;
%Igual que el ejercicio 0
NL(1,1:2:2*length(shapefuns))=shapefuns;
NL(2,2:2:2*length(shapefuns))=shapefuns; %Tiene la forma de las funciones de forma encontradas en el cook pg 206, ecuacion (6.2-2). Despues veo si me sirven
dNL(1,1:2:2*length(shapefuns))=diff(shapefuns,ksi);
dNL(2,2:2:2*length(shapefuns))=diff(shapefuns,eta);
dNLaux=[diff(shapefuns,ksi);diff(shapefuns,eta)]; %Para calcular jacobiano
N=shapefuns;
dN=[diff(N,ksi);diff(N,eta)];

[Nnod, Ndofpornod] = size(nodos);
[Nelem, Nnodporelem] = size(elementos);
doftot = Nnod*Ndofpornod;
DOF = reshape(1:doftot,Ndofpornod,Nnod)';

W=[5/9 8/9 5/9];
wpg=reshape(W.'*W,1,[]);
npg=9;
a=sqrt(.6);
upg=[-a -a;
     -a 0;
     -a a;
      0 -a;
      0 0;
      0 a;
      a -a;
      a 0;
      a a];

doftot = 2*size(nodos,1);
% Nelem = length(elementos);
Nnodporelem = 8;
kG=zeros(doftot);%Matriz global
jmin=1e-10;
% rigid_time=tic;
for e=1:Nelem
    Ke=zeros(Ndofpornod*Nnodporelem);
    index=elementos(e,:);
    elenod=nodos(index,:);
    A=0;
    for ipg=1:npg
        ksi=upg(ipg,1); eta=upg(ipg,2);
        dNs=eval(subs(dN));
        J=dNs*elenod;
        Djac = det(J);
        Ns = eval(subs(N));
        dNxy=J\dNs; %El double (o eval) para que tarde menos, hace que se saltee la algebra (resultado exacto)
        
        r = Ns*elenod(:,1); 
        
        B = zeros(size(Caxi,2),Ndofpornod*Nnodporelem);
        B(1,1:2:Ndofpornod*Nnodporelem-1) = dNxy(1,:);
        B(2,1:2:Ndofpornod*Nnodporelem-1) = Ns(1,:)/r;
        B(3,2:2:Ndofpornod*Nnodporelem) = dNxy(2,:);
        B(4,1:2:Ndofpornod*Nnodporelem-1) = dNxy(2,:);
        B(4,2:2:Ndofpornod*Nnodporelem) = dNxy(1,:);% extra fila por ser axisim.

        Ke = Ke + B'*Caxi*B*wpg(ipg)*Djac;
        %Térmica
%         KTe = KTe + BT'*Ct*BT*wpg(ipg)*Djac;
        
        A = A + wpg(ipg)*Djac;
        %hago la integral me va a servir despues para calcular las fuerzas
%         int(eleDofs,eleDofs)=int(eleDofs,eleDofs)+Nm'*Nm*det(jac);
        if Djac < jmin
            jmin = Djac;
        end
    end
    u=1;
    meindof=reshape(DOF(index,:)',1,[]);
    kG(meindof,meindof)=kG(meindof,meindof)+Ke;
end
end