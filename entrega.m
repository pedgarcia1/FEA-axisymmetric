
clear;close all;clc

set(0,'DefaultFigureWindowStyle','docked')
%todo en SI
%% item A
fid=fopen('mallas/mshes.txt');
casos=textscan(fid,['%s','%s'],'Delimiter',',');



w = 50e-3;

Mo=1;
qmax=3/2/(25e-3)^2;
f1=figure('Name','BandPlots');
err1=figure('Name','error A');



for caso=1:length(casos{1})
    figure(f1)
    subplot(3,2,caso)
    d = w/1.5;
    y = d/2;
    iz = d^3*2e-3/12;
    
    s(caso) = Mo*y/iz;
    
    
    
    txt=strcat('mallas/',casos{1}{caso},'.dat');
    caseLauncherPlaneProblems
    K(caso)=Sfemx(end)/s(caso);
    figure(err1)
    subplot(3,2,caso)
    errorQ8
    err_a(caso)=etaG;
    
%     Scl(caso,:)=[Svm_m;Svm_b;Svm_F];
    
end

rel=[0.03:0.01:0.28];
Kteo=KUbend(50,1.5,rel);


%% 
figure('Name','Kt graph')
a=str2num(cell2mat(casos{2}));

plot(rel,Kteo,'LineWidth',2)
hold on
plot(a,K,'X','Color','r','MarkerSize',15,'LineWidth',2)
Kinter=interp1(a,K,rel,'spline');
ylim([0,5])
ylabel('Kt')
xlabel('r/d')
xlim([0.025,0.3])
plot(rel,Kinter,'Color','r')
legend('Teorico','FEM','FEM interpolado')


% error cuadratico medio
RMSE = sqrt(mean((Kteo(ismember(rel,a))-K).^2));

%imprimir resultados EJ A

fprintf('para un w/d=1.5 y valores de r/d= \n')
fprintf('  %d  ',a)
fprintf('\n se obtuvioron valores de K \n')
fprintf('  %d  ',K)
fprintf('\n comparados con un Kteorico:\n')
fprintf('  %d  ',Kteo(ismember(rel,a)))
fprintf('\n error cuadratico medio: RMSE=%d\n',RMSE)
fprintf('estimacion de error a posteriori \n ')
fprintf('  %d   ',err_a)
fprintf('\n')

%% item b (SCL)
% err2=figure('Name','error B');
M=  27.164288287444386;  %se obtuvo iterando hasta tener tension de fluencia en el concentrador

    qmax=3*M/2/(25e-3)^2;
    d=w/2;y=d/2;
    iz = d^3*2e-3/12;
    
    s=M*y/iz;

    txt=strcat('mallas/','conc_wd_2_rd_0.15','.dat');
  err2=  figure('Name','BandPlot item b');
  subplot(2,1,1)

    caseLauncherPlaneProblems

    a=msh.cord(SCLnod,2)';
    Smg=Sm(1)*ones(1,length(a));
    Sbg=(1-2*a./a(end))*Sb(1);
figure('Name','SCL graph item b')
    plot(a,Smg,'Color','#D95319');
    hold on
    plot(a,Sbg,'m');
    plot(a,Smg+Sbg,'Color','k')
    plot(a,Sfemx,'Color','b')
    plot(a(end),Sm(1)-Sb(1)+SF(1),'X','Color','r','MarkerSize',15,'LineWidth',2)
    legend({'S_m_x','S_b_x','S_m_x+S_b_x','S_x FEM','S_m_x+S_b_x+S_F_x'},'Location','northwest','FontSize',12.5)
    xlabel('x')
    ylabel('\sigma_x','FontSize',20)
    Kfb=Sfemx(end)/(Sm(1)-Sb(1));
    K=Sfemx(end)/s;
    
    fprintf('\n \n \n Item b:\n')
    fprintf('para w/d=2 y r/d=0.15 \n se obtuvo un K de %d y un Kscl de %d \n',[K,Kfb])
    fprintf(' Smx=%d Sbx=%d SFx=%d \n Smy=%d Sby=%d SFy=%d \n Smxy=%d Sbxy=%d SFxy=%d \n',[Sm(1),Sb(1),SF(1),Sm(2),Sb(2),SF(2),Sm(3),Sb(3),SF(3)])
    fprintf('\n\nSvm_m =%d  Svm_b = %d Svm_F = %d\n\n\n',[Svm_m,Svm_b,Svm_F])
       figure(err2)
         subplot(2,1,2)

    errorQ8 
        err_b=etaG;
fprintf('estimacion de error a posteriori \n ')
fprintf('  %d   ',err_b)
fprintf('\n')

%% Item C

% err3=figure('Name','error C');
    M=23.833893857689834;
    qmax=3*M/2/(25e-3)^2;
    

    d=w/2;
    
    s=6*M/d^2/2e-3; 
  

    txt=strcat('mallas/','un_solo_conc_wd_2_rd_0.15','.dat');
  err3=  figure('Name','BandPlot item C');
  subplot(2,1,1)
    caseLauncherunconcentrador
    a=msh.cord(SCLnod,2)';
    a=a-a(1);
     Smg=Sm(1)*ones(1,length(a));
    Sbg=(1-2*a./a(end))*Sb(1);
figure('Name','SCL graph item C')
  plot(a,Smg,'Color','#D95319');
    hold on
    plot(a,Sbg,'m');
    plot(a,Smg+Sbg,'Color','k')
    plot(a,Sfemx,'Color','b')
    plot(a(end),Sm(1)-Sb(1)+SF(1),'X','Color','r','MarkerSize',15,'LineWidth',2)
    legend({'S_m_x','S_b_x','S_m_x+S_b_x','S_x FEM','S_m_x+S_b_x+S_F_x'},'Location','northwest','FontSize',12.5)
    xlabel('x')
    ylabel('\sigma_x','FontSize',20)
    Kfb=Sfemx(end)/(Sm(1)-Sb(1));
    K=Sfemx(end)/s;
    fprintf('\n \n \n Item C:\n')
    fprintf('para w/d=2 y r/d=0.15 con un solo concentrador \n se obtuvo un K de %d y un Kscl de %d \n',[K,Kfb])
    fprintf(' Smx=%d Sbx=%d SFx=%d \n Smy=%d Sby=%d SFy=%d \n Smxy=%d Sbxy=%d SFxy=%d \n',[Sm(1),Sb(1),SF(1),Sm(2),Sb(2),SF(2),Sm(3),Sb(3),SF(3)])
    fprintf('\n\nSvm_m =%d  Svm_b = %d Svm_F = %d\n\n\n',[Svm_m,Svm_b,Svm_F])
        figure(err3)
          subplot(2,1,2)

    errorQ8
    err_c=etaG;
    fprintf('estimacion de error a posteriori \n ')
fprintf('  %d   ',err_c)
fprintf('\n')
