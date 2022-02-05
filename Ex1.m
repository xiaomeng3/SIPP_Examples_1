clear all
mpol x 2
f0=x(1)^2/3+x(1)/2+x(2)^2-x(2);
K1 = [1e4-x(1)^2>=0, 1e4-x(2)^2>=0];
P=msdp(min(f0),K1,3);
[stx,objx]=msol(P);
xsol = double(x)
obju = -100;
mpol u 1
iter = 0;
list_obju = [];

while (obju < -1e-8) 

    g = -xsol(1)^2-2*xsol(1)*xsol(2)*u^2+u-u^3/6+u^5/120;
    K2 = [u>= 0, 2-u>=0];
    newK = [(1-u/2)*(-4*xsol(1)*xsol(2)+1-u^2/2+u^4/24)>=0,...
       (-u/2)*(-4*xsol(1)*xsol(2)+1-u^2/2+u^4/24)>=0 ];
    K2 = [K2, newK];
    P2=msdp(min(g),K2,3);
    [stu,obju]=msol(P2);
    list_obju = [list_obju, obju];
    usol = double(u);

    g1=-x(1)^2-2*x(1)*x(2)*usol^2+usol-usol^3/6+usol^5/120;
    K1 = [K1, g1>=0]
    P=msdp(min(f0),K1,3);
    [stx,objx]=msol(P);
    iter = iter + 1;
    xsol = double(x);
end
list_obju
xsol
iter
objx
 