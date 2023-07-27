Module rochelimit
IMPLICIT NONE
INTEGER,PARAMETER::DP=KIND(1.0D0)
REAL(DP),PARAMETER::pi=3.14159

CONTAINS

REAL(DP) FUNCTION omega(M,r)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r
omega=Sqrt(M/r*r*r)
END FUNCTION omega

REAL(DP) FUNCTION phitidal(M,r,p,q1,q2,q3)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r,p,q1,q2,q3
phitidal=(M/(2*r*r*r))*(-2*(q1*q1)+q2*q2+q3*q3) &
		-(M/(2*(r**4)))*p*q1*(-2*q1*q1+3*q2*q2+3*q3*q3)&
        -(M/(8*(r**5)))*p*p*(8*(q1**4)+3*(q2**4)+3*(q3**4)-24*(q1*q1*q2*q2+q1*q1*q3*q3)+6*q2*q2*q3*q3)
END FUNCTION phitidal

REAL(DP) FUNCTION Ener(M,r,a)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r,a
Ener=(r*r-2*M*r+a*Sqrt(M*r))/(r*Sqrt(r*r-3*M*r+2*a*Sqrt(M*r)))
END FUNCTION Ener

REAL(DP) FUNCTION L(M,r,a)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r,a
L=Sqrt(M*r)*(r*r-2*a*Sqrt(M*r)+a*a)/(r*Sqrt(r*r-3*M*r+2*a*Sqrt(M*r)))
END FUNCTION L

REAL(DP) FUNCTION B(M,r,a)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r,a
B=L(M,r,a)-a*Ener(M,r,a)
END FUNCTION B

REAL(DP) FUNCTION phimag(M,r,a,qg,p,q1,q2,q3)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r,a,qg,p,q1,q2,q3
REAL(DP)::b1,w
b1=B(M,r,a)
w=omega(M,r)
phimag=2*M*b1*Sqrt(1+b1*b1/(r*r))*w*p*(-q1*q1*q1+q1*(q2*q2-q3*q3)&
		+3*qg*(q1*q1-q2*q2)/4)/(r**4)
END FUNCTION phimag

SUBROUTINE genmat(N,A)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::i=0,j=0,k=0
REAL(DP),INTENT(OUT),DIMENSION(N*N*N,N*N*N)::A

Do i=1,N*N*N
  Do j=1,N*N*N
    A(i,j)=0
  End do
End do

Do i=1,N
  Do j=1,N
    Do k=1,N
      
      A((k-1)*N*N+(j-1)*N+i,(k-1)*N*N+(j-1)*N+i)=-6.D0
      If(i<N) THEN
      	A((k-1)*N*N+(j-1)*N+i,(k-1)*N*N+(j-1)*N+i+1)=1.D0
      End If

      If(i>1) THEN
      	A((k-1)*N*N+(j-1)*N+i,(k-1)*N*N+(j-1)*N+i-1)=1.D0
      End If

      If(j<N) THEN
      	A((k-1)*N*N+(j-1)*N+i,(k-1)*N*N+(j)*N+i)=1.D0
      End If

      If(j>1) THEN
      	A((k-1)*N*N+(j-1)*N+i,(k-1)*N*N+(j-2)*N+i)=1.D0
      End If

      If(k<N) THEN
      	A((k-1)*N*N+(j-1)*N+i,(k)*N*N+(j-1)*N+i)= 1.D0
      End If

      If(k>1) THEN
      	A((k-1)*N*N+(j-1)*N+i,(k-2)*N*N+(j-1)*N+i)=1.D0
      End If

    End do
  End do
End do

END SUBROUTINE genmat

SUBROUTINE ludecomp(N,A)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(DP),INTENT(INOUT),DIMENSION(N,N)::A
INTEGER::i=0,j=0,k=0
Real(DP)::temp=0

Do i=2,N
  Write(*,*) 'row = ',i
  Do j=1,N
    
    If(j==1) Then
      A(i,j)=A(i,j)/A(j,j)
      
    Else If(i>j) Then
      temp=0
      Do k=1,j-1
        temp=temp+A(i,k)*A(k,j)
      End do
      A(i,j)=(A(i,j)-temp)/A(j,j)
      
    Else 
      temp=0
      Do k=1,i-1
        temp=temp+A(i,k)*A(k,j)
      End do
      A(i,j)=A(i,j)-temp
    End If

    
  End do
End do

End SUBROUTINE ludecomp

SUBROUTINE lusolve(N,A,X,b)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(DP),INTENT(IN),DIMENSION(N,N)::A
REAL(DP),INTENT(OUT),DIMENSION(N)::X
REAL(DP),INTENT(IN),DIMENSION(N)::b

INTEGER::i=1,j=1,k=1
REAL(DP),DIMENSION(N)::Y
REAL(DP)::temp=0

Do i=1,N
  temp=0
  Do j=1,i-1
  temp=A(i,j)*Y(j)+temp
  End do
  y(i)=b(i)-temp
End do

Do k=N,1,-1
  temp=0
  Do j=k+1,N
    temp=temp+A(k,j)*X(j)
  End do
  X(k)=(Y(k)-temp)/A(k,k)
End do
End SUBROUTINE lusolve

SUBROUTINE gaussiedel(N,A,X,b,iter)
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iter
REAL(DP),INTENT(IN),DIMENSION(N,N)::A
REAL(DP),INTENT(INOUT),DIMENSION(N)::X
REAL(DP),INTENT(IN),DIMENSION(N)::b
INTEGER::i=1,j=1,k=1,i1=1,j1=1
REAL(DP),DIMENSION(N)::Y
REAL(DP)::temp=0
Do i=1,iter
  Do j=1,N
    Write(*,*) 'Row=',j
    temp=0
    Do k=1,j-1
      temp=temp+A(j,k)*Y(k)
    End do
    Do k=j+1,N
      temp=temp+A(j,k)*X(k)
    End do
    Y(j)=(b(j)-temp)/A(j,j)
  End Do
  Do k=1,N
    X(k)=Y(k)
  End Do
End do
End SUBROUTINE gaussiedel

REAL(DP)FUNCTION rho(x1,x2,x3,qs,rhoc)
REAL(DP),INTENT(IN)::x1,x2,x3,qs,rhoc
If(x1*x1+x2*x2+x3*x3>=25) THEN
  rho=0.0
ELSE IF (x1<qs) THEN
  rho=0
Else
  rho=0.15D0
End if
END FUNCTION rho
End MODULE rochelimit

PROGRAM roche
USE rochelimit
IMPLICIT NONE
INTEGER,PARAMETER::N=10
REAL(DP),PARAMETER::M=1,r=0.5,a=0,rhoc=0.15,qs=-.8,kappa=0.5
INTEGER::i=1,j=1,k=1,i1=1,i2=1,j1=1
REAL(DP),DIMENSION(N*N*N,N*N*N)::A
REAL(DP)::xmin=-1.D0,ymin=-1.D0,zmin=-1.D0,h=0.2D0,mass=0,temp=0,w=0,qg=0
REAL(DP),DIMENSION(N+2)::x,y,z
REAL(DP),DIMENSION(N*N*N)::phi,rhs


Write(*,*) 'Filling xyz'
Do i=1,N+2
  x(i)=xmin+h*(i-1)
  y(i)=ymin+h*(i-1)
  z(i)=zmin+h*(i-1)
End Do
Write(*,*) 'mass calculation'    
Do i=1,N+2
  Do j=1,N+2
    Do k=1,N+2
      mass=h*h*h*rho(x(i),y(i),z(i),qs,rhoc)+mass
    End Do
  End Do
End Do
Write(*,*) 'mass =',mass

Write(*,*) 'RHS calculation'
Do i=1,N
  Do j=1,N
    Do k=1,N
      temp=0
      If(i==1) Then
        temp=temp-mass/Sqrt(x(i)*x(i)+y(j+1)*y(j+1)+z(k+1)*z(k+1))
      Else if(i==N) Then
        temp=temp-mass/Sqrt(x(i+2)*x(i+2)+y(j+1)*y(j+1)+z(k+1)*z(k+1))
      End If
      If(j==1) Then
        temp=temp-mass/Sqrt(x(i+1)*x(i+1)+y(j)*y(j)+z(k+1)*z(k+1))
      Else if(j==N) Then
        temp=temp-mass/Sqrt(x(i+1)*x(i+1)+y(j+2)*y(j+2)+z(k+1)*z(k+1))
      End If
      If(j==1) Then
        temp=temp-mass/Sqrt(x(i+1)*x(i+1)+y(j+1)*y(j+1)+z(k)*z(k))
      Else if(j==N) Then
        temp=temp-mass/Sqrt(x(i+1)*x(i+1)+y(j+2)*y(j+2)+z(k+2)*z(k+2))
      End If
      rhs((k-1)*N*N+(j-1)*N+i)=4*pi*h*h*rho(x(i+1),y(j+1),z(k+1),qs,rhoc)+temp      
    End Do
  End Do
End Do

Do i=1,N
  Do j=1,N
    Do k=1,N
      phi((k-1)*N*N+(j-1)*N+i)=-mass/Sqrt(x(1)*x(1)+y(1)*y(1)+z(1)*z(1))
    End Do
  End Do
End Do

CALL genmat(N,A)
CALL ludecomp(N*N*N,A)
CALL lusolve(N*N*N,A,phi,rhs)
w=omega(M,r)
j1=floor(1-xmin/h)
i1=floor(1-xmin/h)+1
i2=floor(1-xmin/h)-1
temp=(1/(2*h))*(phi((j1-2)*N*N+(j1-2)*N+i1-1)-phi((j1-2)*N*N+(j1-2)*N+i2-1))
qg=-temp/(w*w)

END PROGRAM roche