MODULE Amat
INTEGER,PARAMETER::DP=KIND(1.0D0)
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
END MODULE Amat

PROGRAM matrix
USE Amat
IMPLICIT NONE
INTEGER,PARAMETER::N=20
INTEGER::i,j,unit,ierror
REAL(DP),DIMENSION(N*N*N,N*N*N)::A
CHARACTER(len=10)::filename
CALL genmat(N,A)
CALL ludecomp(N*N*N,A)
filename='mat_20.dat'
OPEN(UNIT=unit,FILE=filename,STATUS='replace',ACTION='WRITE',IOSTAT=ierror)
100 FORMAT(F15.5,1X)
DO i=1,N*N*N
  Do j=1,N*N*N
    WRITE(unit,100)A(i,j)
  End Do
End DO
End Program matrix