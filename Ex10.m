clear all 
 
mpol x 2
f0=x(1)+x(2)
%K1 = [1e4-x(1)^2>=0, 1e4-x(2)^2>=0];换成这个就用不了了哦
K1=[x(1)<=100,x(2)<=100,x(1)>=-100,x(2)>=-100]
P=msdp(min(f0),K1);
[stx,objx]=msol(P); 
xsol = double(x)
obju = -100; 
list_obju=[];
iter=0;
mpol u 3

while (obju<-1e-8)
    G=[ 4-xsol(1)^2-xsol(2)^2   xsol(1)   xsol(2);
    xsol(1)           xsol(2)^2-xsol(1)   xsol(1)*xsol(2);
    xsol(2)  xsol(1)*xsol(2)  xsol(1)^2-xsol(2);]
    g=u'*G*u
    K2=[u'*u==1]
    P2=msdp(min(g),K2,4);
    [stu,obju]=msol(P2); 
    usol=double(u)
    usol=usol(:,:,2)
    list_obju=[list_obju,obju]
    iter=iter+1;

    f0=x(1)+x(2)
    G_2=[ 4-x(1)^2-x(2)^2   x(1)   x(2);
        x(1)           x(2)^2-x(1)   x(1)*x(2);
        x(2)  x(1)*x(2)  x(1)^2-x(2);]
    g2=usol'*G_2*usol
    K1 = [1e2-x(1)^2>=0, 1e2-x(2)^2>=0,g2>=0];
    P=msdp(min(f0),K1);
    [stx,objx]=msol(P); 
    xsol = double(x);

end

list_obju
objx
xsol
iter
