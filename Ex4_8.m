clear all
%success
mpol x 1
f0= x^2;
K1 = [x>=1,x<=2];
P=msdp(min(f0),K1);
[stx,objx]=msol(P);
xsol = double(x)
obju = -100;
mpol u 2
iter = 0;
list_obju = [];

while (obju<-1e-8)
    g=xsol*(u(1)-u(2)+1);
    K2=[u(1)^2*(u(1)-u(2)+1)>=0, -7*u(1)*(3*u(1)^2-2*u(1)*u(2))==1,...
     -7*u(1)^3==-1];
    P2=msdp(min(g),K2,2);
    [stu,obju]=msol(P2);
    usol=double(u);
    list_obju=[list_obju,obju];

    g1=x*(usol(1)-usol(2)+1)
    K1=[K1,g1>=0]
    P=msdp(min(f0),K1)
    [stx,objx]=msol(P)
    xsol=double(x);
    iter=iter+1;
end

list_obju
xsol
iter
usol
xsol
