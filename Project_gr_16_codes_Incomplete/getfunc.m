function [kappa,w,D,zeta,B]=getfunc(M,rhoc,r,n,a,R)
format long
if n==0.5
    Y=2.75270
elseif n==1.0
    Y=pi
elseif n==1.5
    Y=3.65375
end
w=sqrt(M/(r.^3))
zeta=M/(pi*rhoc*r.^3)
D=sqrt(r*r-3*M*r+2*a*sqrt(M*r))
B=r*(sqrt(M*r)-a)/D
kappa=4*pi*(R/Y).^2*rhoc.^(1-1/n)
end
