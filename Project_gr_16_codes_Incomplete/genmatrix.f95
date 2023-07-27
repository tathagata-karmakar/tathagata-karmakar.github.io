MODULE poissolve
IMPLICIT NONE
SAVE

INTEGER,PARAMETER::DP=KIND(1.D0)


CONTAINS

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
INTEGER::i=1,j=1,k=1
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



REAL(DP)FUNCTION rho(x1,x2,x3)
REAL(DP),INTENT(IN)::x1,x2,x3
If(x1*x1+x2*x2+x3*x3>=25) THEN
  rho=0.0
Else
  rho=0.15D0
End if
END FUNCTION rho
 
End Module poissolve

PROGRAM genmatrix
Use poissolve
IMPLICIT NONE
INTEGER,PARAMETER::N=15
INTEGER::i=1,j=1,k=1
REAL(DP),DIMENSION(N*N*N,N*N*N)::A
REAL(DP)::xmin=-10.D0,ymin=-10.D0,zmin=-10.D0,h=2,mass=0,temp=0
REAL(DP),DIMENSION(N+2)::x,y,z
REAL(DP),DIMENSION(N*N*N)::phi,rhs
        
REAL(DP),DIMENSION(3,3)::B=reshape((/1,2,4,3,8,14,2,6,13/),shape(B),order=(/2,1/))
REAL(DP),DIMENSION(3)::X1=(/1,1,1/),b1=(/3,13,4/)

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
      mass=h*h*h*rho(x(i),y(i),z(i))+mass
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
        temp=temp-mass/Sqrt(x(i)*x(i)+y(j+1)*(j+1)+z(k+1)*z(k+1))
      Else if(i==N) Then
        temp=temp-mass/Sqrt(x(i+2)*x(i+2)+y(j+1)*(j+1)+z(k+1)*z(k+1))
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
      rhs((k-1)*N*N+(j-1)*N+i)=h*h*rho(x(i+1),y(j+1),z(k+1))+temp      
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



Write(*,*) 'Genarate A'
Call genmat(N,A)

!Write(*,*) 'LU decomp of A'
CALL ludecomp(N*N*N,A)

Write(*,*) 'Solving for potential'
CALL gaussiedel(N*N*N,A,phi,rhs,1)
      

!Call gaussiedel(3,B,X1,b1,50)
!CALL genmat(N,A)
!CALL ludecomp(3,B)
!CALL lusolve(3,B,X1,b1)
!Do i=1,3
! Do j=1,N*N*N
!    WRITE(*,'(F8.4,1X)',advance='no') X1(i)
!  End do
!  Write(*,*)
!End do
!CALL genmat(N,A)


END PROGRAM genmatrix