format long
M=1;r=6*M;R=.116478*M;n=0.5;a=M;qs=-0.08;N=50; delta=abs(qs)/40;
fileId=fopen('abcd.dat','w')
syms x1 x2 x3 q1 q2 q3 p C qg rho
for rhoc=.015:.01:.025
    [kappa,omega,D,zeta,B]=getfunc(M,rhoc,r,n,a,R);
    p1=1;
    f(p,q1,q2,q3)=simplify(-2*pi*R*R*rhoc*(1/(p*p))*(1-p*p*(q1.^2+q2.^2+...
        q3.^2)/(3*R*R)));
    for i=1:1:2
        rho=eval((1/(4*pi))*simplify(diff(f(p,q1,q2,q3),q1,2)+...
            diff(f(p,q1,q2,q3),q2,2)+...
            diff(f(p,q1,q2,q3),q3,2)))
        qg2=-diff(f(p,q1,q2,q3),q1)/omega.^2
        qg1=subs(qg2,{p,q1,q2,q3},{p1,0,0,0})
        r1=vpasolve(omega.^2*p.^2*(qg1.^2-(qs-qg1).^2)==kappa*(n+1)*(rhoc.^(1/n)...
            -subs(rho,{q1,q2,q3},{qs,0,0}).^(1/n))+p.^2*(f(p,0,0,0)...
            -f(p,qs,0,0))-ftidal(p*qs,0,0,M,r)-fmag(p*qs,0,0,M,r,a,p,qg1),...
            p)
    end
     
end