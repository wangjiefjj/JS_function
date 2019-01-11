function [theta_a,theta_e,fd_s]=function_partion_APF(Ns,Ms,K,d,V,lemda,fr)
a = Ns*K;
b = min(Ns,K);
while mod(a,2) ~= 0
    zz = a;
    break;
end
c = a/2;
while mod(c,2) == 0    
    c = c/2;
    if mod(c,Ns) ~= 0 || mod(c,K) ~= 0
        zz = c * 2;
        break;
    end
end
xx = floor(zz*lemda/2/d);
yy = floor(zz*fr/(4*V/lemda));

theta_a = [-1:2/(xx):1];
theta_e = [0:1/(Ms):1]; 
fd_s = [-1:2/(yy):1];