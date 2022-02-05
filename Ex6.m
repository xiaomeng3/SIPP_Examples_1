% Example 6
clear all 
% Solve argmin f(x) 
mpol x 3
f0=(x(1)+x(2)-2)^2+(x(1)-x(2))^2+30*x(3)^2;
K1 = [1e4-x(1)^2>=0, 1e4-x(2)^2>=0, x(3)>=0,x(3)>=x(1)-x(2),...
    x(3)*(x(3)-x(1)+x(2))==0];
P=msdp(min(f0),K1);
[stx,objx]=msol(P); 
xsol = double(x)
obju = -100; 
mpol u 2
iter = 0;
list_obju = [];
while (obju <-1e-8) 
%%%%%%%%%
    g = -xsol(1)*u(2)-xsol(2)*u(1)+1;
    K2= [u(1)>=0, u(1)^2+u(2)^2==1]
    lambda1=-u(2)*xsol(1)/2-u(1)*xsol(2)/2;
    lambda2=u(1)*u(2)*xsol(2)-xsol(1)+u(2)^2*xsol(1)
    newK=[2*u(1)*lambda1+lambda2+xsol(2)==0,...
          2*lambda1*u(2)+xsol(1)==0,...
          lambda2*u(1)==0,...
          lambda2>=0]
    K2=[K2,newK]     
    P2=msdp(min(g),K2,3);
    [stu,obju]=msol(P2); 
    list_obju = [list_obju, obju];
    usol = double(u);
    g1 = -x(1)*usol(2)-x(2)*usol(1)+1;
    K1 = [K1, g1>=0];
    P=msdp(min(f0),K1,2);
    [stx,objx]=msol(P); 
    iter = iter + 1;
    xsol = double(x);
%%%%%%%%%%%
end
iter,
list_obju,
objx
obju = double(g1),
xsol
usol