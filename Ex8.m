clear all
mpol x 2
f0=(x(1)+1)^2+(x(2)+1)^2
K1 = [1e4-x(1)^2>=0, 1e4-x(2)^2>=0];
P=msdp(min(f0),K1);
[stx,objx]=msol(P); 
xsol = double(x)
obju = -100; 
list_obju=[];
iter=0;
while (obju < -1e-4)
    mpol u 2
    x=xsol
    p0=-(x(1)^2+u(1)^2*x(2)^2+2*u(1)*u(2)*x(1)*x(2)+x(1)+x(2))
    daoshu1=-2*x(2)^2*u(1)-2*u(2)*x(2)*x(1)
    daoshu2=-2*x(1)*x(2)*u(1)
    K1=[1-u(1)^2>=0, 1-u(2)^2>=0]
    newK=[(-u(1)/2)*daoshu1>=0,  (-u(2)/2)*daoshu2>=0]
    K1=[K1,newK]
    P=msdp(min(p0),K1,2)
    [stu,obju]=msol(P)
    u=double(u)
    u=u(:,:,2)
    list_obju=[list_obju,obju];

    mpol x 2
    f0=(x(1)+1)^2+(x(2)+1)^2
    p0=-(x(1)^2+u(1)^2*x(2)^2+2*u(1)*u(2)*x(1)*x(2)+x(1)+x(2))
    K1 = [1e4-x(1)^2>=0, 1e4-x(2)^2>=0,p0>=0];
    P=msdp(min(f0),K1);
    [stx,objx]=msol(P); 
    xsol = double(x)
    iter=iter+1;

end 

objx
list_obju
iter

%%%%%%%%%
