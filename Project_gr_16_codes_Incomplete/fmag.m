function y=fmag(x1,x2,x3,M,r,a,p,qg)
B=r*(sqrt(M*r)-a)/sqrt(r*r-3*M*r+2*a*sqrt(M*r))
w=sqrt(M/(r.^3))
y=2*M*(B/r.^4)*sqrt(1+B.^2/r.^2)*w*(-x1.^3 + x1*(x2.^2 - x3.^2)...
    +(3/4)*p.^2*qg*(x1.^2 - x2.^2));
end