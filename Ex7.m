% Example 7
clear all 
% Solve argmin f(x) 
mpol x 2
f0 = x(2)
K1 = [x(1)>=-100,x(1)<=100,x(2)>=-100,x(2)<=100];
P=msdp(min(f0),K1);
[stx,objx]=msol(P); 
xsol = double(x)
obju = -100; 
mpol u 1
iter = 0;
list_obju = [];
U=[];
while (obju <-1e-8) 
%%%%%%%%%
    g = -2*xsol(1)^2*u^2+u^4-xsol(1)^2+xsol(2);
    K2 = [u>= -1, 1 >= u];
    P2=msdp(min(g),K2,3);
    [stu,obju]=msol(P2); 
    list_obju = [list_obju, obju];
    usol = double(u); 
    U = [U,usol];
    for k = 1:length(U)
        u1 =U(k);
        g1 = -2*x(1)^2*u1^2+u1^4-x(1)^2+x(2);
        K1 = [K1, g1>=0];
    end
    P=msdp(min(f0),K1,3);
    [stx,objx]=msol(P); 
    iter = iter + 1;
    xsol = double(x);
%%%%%%%%%%
end
iter,
list_obju,
objx,
obju = double(g1),
xsol
return,