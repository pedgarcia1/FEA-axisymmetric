
E=200000;
v=0;
a=1000;
b=5000;
p=200;


ur=@(r)  p*(a^2*(1+v)*(b^2+r^2*(1-2*v))/(E*(b^2-a^2)*r));
srr=@(r) p*(a^2/(b^2-a^2))*(1-b^2/r^2);
stt=@(r) p*(a^2/(b^2-a^2))*(1+b^2/r^2);
szz=@(r) p*(2*a^2*v)/(b^2-a^2);