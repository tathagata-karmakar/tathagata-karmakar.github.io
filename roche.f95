Module rochelimit
IMPLICIT NONE
INTEGER,PARAMETER::DP=KIND(1.0D0)
REAL(DP),PARAMETER::pi=3.14159D0

CONTAINS

REAL(DP) FUNCTION omega(M,r)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r
omega=Sqrt(M/(r*r*r))
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

REAL(DP) FUNCTION alpha(kappa,n1,rhoc)
IMPLICIT NONE
REAL(DP),INTENT(IN)::kappa,n1,rhoc
alpha=Sqrt((n1+1)*kappa*(rhoc**(1/n1-1))/(4*pi))
END Function alpha

REAL(DP) FUNCTION lane(n1)
IMPLICIT NONE
REAL(DP)::n1
If(Abs(n1-0.5D0)<.0001) Then
  lane=2.75270D0
Else If(Abs(n1-1.D0)<.0001)Then
  lane=pi
Else If(Abs(n1-0.5D0)<.0001)Then
  lane=lane*3.65375D0
End If
End FUNCTION lane

REAL(DP) FUNCTION phimag(M,r,a,qg,p,q1,q2,q3)
IMPLICIT NONE
REAL(DP),INTENT(IN)::M,r,a,qg,p,q1,q2,q3
REAL(DP)::b1,w
b1=B(M,r,a)
w=omega(M,r)
phimag=2*M*b1*Sqrt(1+b1*b1/(r*r))*w*p*(-q1*q1*q1+q1*(q2*q2-q3*q3)&
		+3*qg*(q1*q1-q2*q2)/4)/(r**4)
END FUNCTION phimag

REAL(DP) FUNCTION AA1(N,l1,m)
IMPLICIT NONE
INTEGER,INTENT(IN)::N,l1,m
INTEGER::i=0,j=0,k=0
!REAL(DP),INTENT(OUT),DIMENSION(N*N*N,N*N*N)::A

!Do i=1,N*N*N
!  Do j=1,N*N*N
!    A(i,j)=0
!  End do
!End do

IF(l1==m) Then
   AA1=-6.D0
   RETURN
END IF 
k=1+l1/(N*N)
j=(l1-(k-1)*N*N)/N+1
i=l1-((k-1)*N*N+(j-1)*N)
      
      If((i<N).and.(m==(k-1)*N*N+(j-1)*N+i+1)) THEN
      	AA1 = 1.D0
        RETURN
      End If

      If((i>1).and.(m==(k-1)*N*N+(j-1)*N+i-1)) THEN
      	AA1=1.D0
        RETURN
      End If

      If((j<N).and.(m==(k-1)*N*N+(j)*N+i)) THEN
      	AA1=1.D0
        RETURN
      End If

      If((j>1).and.(m==(k-1)*N*N+(j-2)*N+i)) THEN
      	AA1=1.D0
        RETURN
      End If

      If((k<N).and.(m==(k)*N*N+(j-1)*N+i)) THEN
      	AA1= 1.D0
        RETURN
      End If

      If((k>1).and.(m==(k-2)*N*N+(j-1)*N+i)) THEN
      	AA1=1.D0
        RETURN
      End If
      
AA1=0.D0

END FUNCTION AA1


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
  !Write(*,*) 'row = ',i
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

SUBROUTINE gaussiedel1(N,X,b,iter,N1)
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iter,N1
!REAL(DP),INTENT(IN),DIMENSION(N,N)::A
REAL(DP),INTENT(INOUT),DIMENSION(N)::X
REAL(DP),INTENT(IN),DIMENSION(N)::b
INTEGER::i=1,j=1,k=1,i1=1,j1=1
REAL(DP),DIMENSION(N)::Y
REAL(DP)::temp=0
Do i=1,iter
  Do j=1,N
    !Write(*,*) 'Row=',j
    temp=0
    Do k=1,j-1
      temp=temp+AA1(N1,j,k)*Y(k)
    End do
    Do k=j+1,N
      temp=temp+AA1(N1,j,k)*X(k)
    End do
    Y(j)=(b(j)-temp)/AA1(N1,j,j)
  End Do
  Do k=1,N
    X(k)=Y(k)
  End Do
End do
End SUBROUTINE gaussiedel1

REAL(DP)FUNCTION rho(x1,x2,x3,qs,rhoc,kappa,n1)
REAL(DP),INTENT(IN)::x1,x2,x3,qs,rhoc,kappa,n1
REAL(DP)::temp,r,temp1
temp=alpha(kappa,n1,rhoc)
r=Sqrt(x1*x1+x2*x2+x3*x3)
If((r/temp)>=lane(n1)) THEN
  rho=0.0
ELSE IF ((Abs(x1)>Abs(qs)).or.(Abs(x2)>Abs(qs)).or.(Abs(x3)>Abs(qs))) THEN
  rho=0.D0
Else
  temp1=Sin(r/temp)*temp/r
  If(temp1>0.D0) Then
  rho=rhoc*(temp1**n1)
  Else 
    rho=0.D0
  End If
End if
END FUNCTION rho
End MODULE rochelimit

MODULE PolynomialRoots


  IMPLICIT NONE
  INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0), DP=KIND(1.0D0)
  REAL(DP),PARAMETER,PRIVATE:: ZERO=0.0D0, FOURTH=0.25D0, HALF=0.5D0
  REAL(DP),PARAMETER,PRIVATE:: ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, FOUR=4.0D0
  COMPLEX(DP),PARAMETER,PRIVATE:: CZERO=(0.D0,0.D0)

  REAL(DP),PARAMETER,PRIVATE:: EPS=EPSILON(ONE)

  CHARACTER(LEN=*),PARAMETER,PUBLIC:: POLYROOTS_VERSION= "1.3 (4 Jan 1999)"
  INTEGER,PRIVATE:: outputCode=100

  PRIVATE:: CubeRoot
  PUBLIC:: LinearRoot
  PRIVATE:: OneLargeTwoSmall
  PUBLIC:: QuadraticRoots
  PUBLIC:: CubicRoots
  PUBLIC:: QuarticRoots
  PUBLIC:: SolvePolynomial
!----------------------------------------------------------------------------

  INTERFACE Swap
    MODULE PROCEDURE SwapDouble, SwapSingle
  END INTERFACE

CONTAINS

!+
FUNCTION CubeRoot(x) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the Cube Root of a REAL(DP) number. If the argument is
!   negative, then the cube root is also negative.

  REAL(DP),INTENT(IN) :: x
  REAL(DP):: f
!----------------------------------------------------------------------------
  IF (x < ZERO) THEN
    f=-EXP(LOG(-x)/THREE)
  ELSE IF (x > ZERO) THEN
    f=EXP(LOG(x)/THREE)
  ELSE
    f=ZERO
  END IF
  RETURN
END Function CubeRoot   ! ---------------------------------------------------

!+
SUBROUTINE LinearRoot(a, z)
IMPLICIT NONE
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z 
!     AND STORES THE RESULTS IN Z. It is assumed that a(2) is non-zero.
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  REAL(DP),INTENT(OUT):: z
!----------------------------------------------------------------------------
  IF (a(2)==0.0) THEN
    z=ZERO
  ELSE
    z=-a(1)/a(2)
  END IF
  RETURN
END Subroutine LinearRoot   ! -----------------------------------------------

!+
SUBROUTINE OneLargeTwoSmall(a1,a2,a4,w, z)
IMPLICIT NONE
! ---------------------------------------------------------------------------
! PURPOSE - Compute the roots of a cubic when one root, w, is known to be
!   much larger in magnitude than the other two

  REAL(DP),INTENT(IN):: a1,a2,a4
  REAL(DP),INTENT(IN):: w
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


  REAL(DP),DIMENSION(3):: aq
!----------------------------------------------------------------------------
  aq(1)=a1
  aq(2)=a2+a1/w
  aq(3)=-a4*w
  CALL QuadraticRoots(aq, z)
  z(3)=CMPLX(w,ZERO,DP)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3)=z(2)
  z(2)=z(1)
  z(1)=CMPLX(w,ZERO,DP)
  RETURN
END Subroutine OneLargeTwoSmall   ! -----------------------------------------

!+
SUBROUTINE QuadraticRoots(a, z)
IMPLICIT NONE
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.

  REAL(DP),INTENT(IN),DIMENSION(:):: a
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


  REAL(DP):: d, r, w, x, y
!----------------------------------------------------------------------------
  IF(a(1)==0.0) THEN     ! EPS is a global module constant (private)
    z(1) = CZERO               ! one root is obviously zero
    z(2) = CMPLX(-a(2)/a(3), ZERO,DP)    ! remainder is a linear eq.
    outputCode=21   ! two identical real roots
    RETURN
  END IF

  d = a(2)*a(2) - FOUR*a(1)*a(3)             ! the discriminant

  IF (ABS(d) <= TWO*eps*a(2)*a(2)) THEN
    z(1) = CMPLX(-HALF*a(2)/a(3), ZERO, DP) ! discriminant is tiny
    z(2) = z(1)
    outputCode=22  ! two distinct real roots
    RETURN
  END IF

  r = SQRT(ABS(d))
  IF (d < ZERO) THEN
    x = -HALF*a(2)/a(3)        ! negative discriminant => roots are complex   
    y = ABS(HALF*r/a(3))
    z(1) = CMPLX(x, y, DP)
    z(2) = CMPLX(x,-y, DP)   ! its conjugate
    outputCode=23                        !  COMPLEX ROOTS
    RETURN
  END IF

  IF (a(2) /= ZERO) THEN              ! see Numerical Recipes, sec. 5.5
    w = -(a(2) + SIGN(r,a(2)))
    z(1) = CMPLX(TWO*a(1)/w,  ZERO, DP)
    z(2) = CMPLX(HALF*w/a(3), ZERO, DP)
    outputCode=22           ! two real roots
    RETURN
  END IF

  x = ABS(HALF*r/a(3))   ! a(2)=0 if you get here
  z(1) = CMPLX( x, ZERO, DP)
  z(2) = CMPLX(-x, ZERO, DP)
  outputCode=22
  RETURN
END Subroutine QuadraticRoots   ! -------------------------------------------

!+
SUBROUTINE CubicRoots(a, z)
IMPLICIT NONE
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z

  REAL(DP),PARAMETER:: RT3=1.7320508075689D0    ! (Sqrt(3)
  REAL (DP) :: aq(3), arg, c, cf, d, p, p1, q, q1
  REAL(DP):: r, ra, rb, rq, rt
  REAL(DP):: r1, s, sf, sq, sum, t, tol, t1, w
  REAL(DP):: w1, w2, x, x1, x2, x3, y, y1, y2, y3

! NOTE -   It is assumed that a(4) is non-zero. No test is made here.
!----------------------------------------------------------------------------
  IF (a(1)==0.0) THEN
    z(1) = CZERO  ! one root is obviously zero
    CALL QuadraticRoots(a(2:4), z(2:3))   ! remaining 2 roots here
    RETURN
  END IF

  p = a(3)/(THREE*a(4))
  q = a(2)/a(4)
  r = a(1)/a(4)
  tol = FOUR*EPS

  c = ZERO
  t = a(2) - p*a(3)
  IF (ABS(t) > tol*ABS(a(2))) c = t/a(4)

  t = TWO*p*p - q
  IF (ABS(t) <= tol*ABS(q)) t = ZERO
  d = r + p*t
  IF (ABS(d) <= tol*ABS(r)) GO TO 110

!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

  s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
  p1 = a(3)/(THREE*s)
  q1 = a(2)/s
  r1 = a(1)/s

  t1 = q - 2.25D0*p*p
  IF (ABS(t1) <= tol*ABS(q)) t1 = ZERO
  w = FOURTH*r1*r1
  w1 = HALF*p1*r1*t
  w2 = q1*q1*t1/27.0D0

  IF (w1 >= ZERO) THEN
    w = w + w1
    sq = w + w2
  ELSE IF (w2 < ZERO) THEN
    sq = w + (w1 + w2)
  ELSE
    w = w + w2
    sq = w + w1
  END IF

  IF (ABS(sq) <= tol*w) sq = ZERO
  rq = ABS(s/a(4))*SQRT(ABS(sq))
  IF (sq >= ZERO) GO TO 40

!                   ALL ROOTS ARE REAL

  arg = ATAN2(rq, -HALF*d)
  cf = COS(arg/THREE)
  sf = SIN(arg/THREE)
  rt = SQRT(-c/THREE)
  y1 = TWO*rt*cf
  y2 = -rt*(cf + rt3*sf)
  y3 = -(d/y1)/y2

  x1 = y1 - p
  x2 = y2 - p
  x3 = y3 - p

  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)

  IF (ABS(x2) > ABS(x3)) CALL Swap(x2,x3)
  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)

  w = x3

  IF (ABS(x2) < 0.1D0*ABS(x3)) GO TO 70
  IF (ABS(x1) < 0.1D0*ABS(x2)) x1 = - (r/x3)/x2
  z(1) = CMPLX(x1, ZERO,DP)
  z(2) = CMPLX(x2, ZERO,DP)
  z(3) = CMPLX(x3, ZERO,DP)
  RETURN

!                  REAL AND COMPLEX ROOTS

40 ra =CubeRoot(-HALF*d - SIGN(rq,d))
  rb = -c/(THREE*ra)
  t = ra + rb
  w = -p
  x = -p
  IF (ABS(t) <= tol*ABS(ra)) GO TO 41
  w = t - p
  x = -HALF*t - p
  IF (ABS(x) <= tol*ABS(p)) x = ZERO
  41 t = ABS(ra - rb)
  y = HALF*rt3*t
  
  IF (t <= tol*ABS(ra)) GO TO 60
  IF (ABS(x) < ABS(y)) GO TO 50
  s = ABS(x)
  t = y/x
  GO TO 51
50 s = ABS(y)
  t = x/y
51 IF (s < 0.1D0*ABS(w)) GO TO 70
  w1 = w/s
  sum = ONE + t*t
  IF (w1*w1 < 0.01D0*sum) w = - ((r/sum)/s)/s
  z(1) = CMPLX(w, ZERO,DP)
  z(2) = CMPLX(x, y,DP)
  z(3) = CMPLX(x,-y,DP)
  RETURN

!               AT LEAST TWO ROOTS ARE EQUAL

60 IF (ABS(x) < ABS(w)) GO TO 61
  IF (ABS(w) < 0.1D0*ABS(x)) w = - (r/x)/x
  z(1) = CMPLX(w, ZERO,DP)
  z(2) = CMPLX(x, ZERO,DP)
  z(3) = z(2)
  RETURN
  61 IF (ABS(x) < 0.1D0*ABS(w)) GO TO 70
  z(1) = CMPLX(x, ZERO,DP)
  z(2) = z(1)
  z(3) = CMPLX(w, ZERO,DP)
  RETURN

!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

70 aq(1) = a(1)
  aq(2) = a(2) + a(1)/w
  aq(3) = -a(4)*w
  CALL QuadraticRoots(aq, z)
  z(3) = CMPLX(w, ZERO,DP)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3) = z(2)
  z(2) = z(1)
  z(1) = CMPLX(w, ZERO,DP)
  RETURN
!-----------------------------------------------------------------------


!                   CASE WHEN D = 0

110 z(1) = CMPLX(-p, ZERO,DP)
  w = SQRT(ABS(c))
  IF (c < ZERO) GO TO 120
  z(2) = CMPLX(-p, w,DP)
  z(3) = CMPLX(-p,-w,DP)
  RETURN

120 IF (p /= ZERO) GO TO 130
  z(2) = CMPLX(w, ZERO,DP)
  z(3) = CMPLX(-w, ZERO,DP)
  RETURN

130 x = -(p + SIGN(w,p))
  z(3) = CMPLX(x, ZERO,DP)
  t = THREE*a(1)/(a(3)*x)
  IF (ABS(p) > ABS(t)) GO TO 131
  z(2) = CMPLX(t, ZERO,DP)
  RETURN
131 z(2) = z(1)
  z(1) = CMPLX(t, ZERO,DP)
  RETURN
END Subroutine CubicRoots   ! -----------------------------------------------


!+
SUBROUTINE QuarticRoots(a,z)
IMPLICIT NONE
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!               A(1) + A(2)*Z + ... + A(5)*Z**4

  REAL(DP), INTENT(IN)     :: a(:)
  COMPLEX(DP), INTENT(OUT) :: z(:)

  COMPLEX(DP) :: w
  REAL(DP):: b,b2, c, d, e, h, p, q, r, t
  REAL(DP),DIMENSION(4):: temp

  REAL(DP):: u, v, v1, v2, x, x1, x2, x3, y


! NOTE - It is assumed that a(5) is non-zero. No test is made here

!----------------------------------------------------------------------------

  IF (a(1)==0.0) THEN
    z(1) = CZERO    !  one root is obviously zero
    CALL CubicRoots(a(2:), z(2:))
    RETURN
  END IF


  b = a(4)/(FOUR*a(5))
  c = a(3)/a(5)
  d = a(2)/a(5)
  e = a(1)/a(5)
  b2 = b*b

  p = HALF*(c - 6.0D0*b2)
  q = d - TWO*b*(c - FOUR*b2)
  r = b2*(c - THREE*b2) - b*d + e

! SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
! NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
! THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE
!     Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
! WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
! THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.

  temp(1) = -q*q/64.0D0
  temp(2) = 0.25D0*(p*p - r)
  temp(3) =  p
  temp(4) = ONE
  CALL CubicRoots(temp,z)
  IF (AIMAG(z(2)) /= ZERO) GO TO 60

!         THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
!         REORDER THE ROOTS IN INCREASING ORDER

  x1 = DBLE(z(1))
  x2 = DBLE(z(2))
  x3 = DBLE(z(3))
  IF (x1 > x2) CALL Swap(x1,x2)
  IF (x2 > x3) CALL Swap(x2,x3)
  IF (x1 > x2) CALL Swap(x1,x2)

  u = ZERO
  IF (x3 > ZERO) u = SQRT(x3)
  IF (x2 <= ZERO) GO TO 41
  IF (x1 >= ZERO) GO TO 30
  IF (ABS(x1) > x2) GO TO 40
  x1 = ZERO

30 x1 = SQRT(x1)
  x2 = SQRT(x2)
  IF (q > ZERO) x1 = -x1
  temp(1) = (( x1 + x2) + u) - b
  temp(2) = ((-x1 - x2) + u) - b
  temp(3) = (( x1 - x2) - u) - b
  temp(4) = ((-x1 + x2) - u) - b
  CALL SelectSort(temp)
  IF (ABS(temp(1)) >= 0.1D0*ABS(temp(4))) GO TO 31
  t = temp(2)*temp(3)*temp(4)
  IF (t /= ZERO) temp(1) = e/t
31 z(1) = CMPLX(temp(1), ZERO,DP)
  z(2) = CMPLX(temp(2), ZERO,DP)
  z(3) = CMPLX(temp(3), ZERO,DP)
  z(4) = CMPLX(temp(4), ZERO,DP)
  RETURN

40 v1 = SQRT(ABS(x1))
v2 = ZERO
GO TO 50
41 v1 = SQRT(ABS(x1))
v2 = SQRT(ABS(x2))
IF (q < ZERO) u = -u

50 x = -u - b
y = v1 - v2
z(1) = CMPLX(x, y,DP)
z(2) = CMPLX(x,-y,DP)
x =  u - b
y = v1 + v2
z(3) = CMPLX(x, y,DP)
z(4) = CMPLX(x,-y,DP)
RETURN

!                THE RESOLVENT CUBIC HAS COMPLEX ROOTS

60 t = DBLE(z(1))
x = ZERO
IF (t < ZERO) THEN
  GO TO 61
ELSE IF (t == ZERO) THEN
  GO TO 70
ELSE
  GO TO 62
END IF
61 h = ABS(DBLE(z(2))) + ABS(AIMAG(z(2)))
IF (ABS(t) <= h) GO TO 70
GO TO 80
62 x = SQRT(t)
IF (q > ZERO) x = -x

70 w = SQRT(z(2))
  u = TWO*DBLE(w)
  v = TWO*ABS(AIMAG(w))
  t =  x - b
  x1 = t + u
  x2 = t - u
  IF (ABS(x1) <= ABS(x2)) GO TO 71
  t = x1
  x1 = x2
  x2 = t
71 u = -x - b

  h = u*u + v*v
  IF (x1*x1 < 0.01D0*MIN(x2*x2,h)) x1 = e/(x2*h)
  z(1) = CMPLX(x1, ZERO,DP)
  z(2) = CMPLX(x2, ZERO,DP)
  z(3) = CMPLX(u, v,DP)
  z(4) = CMPLX(u,-v,DP)
  RETURN

80 v = SQRT(ABS(t))
  z(1) = CMPLX(-b, v,DP)
  z(2) = CMPLX(-b,-v,DP)
  z(3) = z(1)
  z(4) = z(2)
  RETURN

END Subroutine QuarticRoots

!+
SUBROUTINE SelectSort(a)
IMPLICIT NONE
! ---------------------------------------------------------------------------
! PURPOSE - Reorder the elements of in increasing order.
  REAL(DP),INTENT(IN OUT),DIMENSION(:):: a

  INTEGER:: j
  INTEGER,DIMENSION(1):: k
! NOTE - This is a n**2 method. It should only be used for small arrays. <25
!----------------------------------------------------------------------------
  DO j=1,SIZE(a)-1
    k=MINLOC(a(j:))
    IF (j /= k(1)) CALL Swap(a(k(1)),a(j))
  END DO
  RETURN
END Subroutine SelectSort   ! -----------------------------------------------

!+
SUBROUTINE SolvePolynomial(quarticCoeff, cubicCoeff, quadraticCoeff, &
  linearCoeff, constantCoeff, code1, root1,root2,root3,root4)
IMPLICIT NONE
! ---------------------------------------------------------------------------
  REAL(DP),INTENT(IN):: quarticCoeff
  REAL(DP),INTENT(IN):: cubicCoeff, quadraticCoeff
  REAL(DP),INTENT(IN):: linearCoeff, constantCoeff
  INTEGER,INTENT(OUT):: code1
  COMPLEX(DP),INTENT(OUT):: root1,root2,root3,root4
  REAL(DP),DIMENSION(5):: a
  COMPLEX(DP),DIMENSION(5):: z
!----------------------------------------------------------------------------
  a(1)=constantCoeff
  a(2)=linearCoeff
  a(3)=quadraticCoeff
  a(4)=cubicCoeff
  a(5)=quarticCoeff

  IF (quarticCoeff /= ZERO) THEN
    CALL QuarticRoots(a,z)  
  ELSE IF (cubicCoeff /= ZERO) THEN
    CALL CubicRoots(a,z)
  ELSE IF (quadraticCoeff /= ZERO) THEN
    CALL QuadraticRoots(a,z)
  ELSE IF (linearCoeff /= ZERO) THEN
    z(1)=CMPLX(-constantCoeff/linearCoeff, 0, DP)
    outputCode=1
  ELSE
    outputCode=0    !  { no roots }
  END IF

  code1=outputCode
  IF (outputCode > 0) root1=z(1)
  IF (outputCode > 1) root2=z(2)
  IF (outputCode > 23) root3=z(3)
  IF (outputCode > 99) root4=z(4)
  RETURN
END Subroutine SolvePolynomial   ! ------------------------------------------

!+
SUBROUTINE SwapDouble(a,b)
IMPLICIT NONE
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(DP),INTENT(IN OUT):: a,b
  REAL(DP):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine SwapDouble   ! -----------------------------------------------

!+
SUBROUTINE SwapSingle(a,b)
IMPLICIT NONE
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(SP),INTENT(IN OUT):: a,b
  REAL(SP):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine SwapSingle   ! -----------------------------------------------


END Module PolynomialRoots   ! ==============================================



  
PROGRAM roche
USE rochelimit
USE PolynomialRoots
IMPLICIT NONE
INTEGER,PARAMETER::N=20
REAL(DP),PARAMETER::M=1.D0,r=20.D0,a1=0.D0,rhoc=.065D0,qs=-0.80D0,tol=.0001D0,n1=1.0D0,kappa=1+1/n1
INTEGER::i=1,j=1,k=1,code1=0,noi=1,iter=10
INTEGER::i1x=1,i2x=1,i0x=1,i1y=1,i2y=1,i0y=1,i1z=1,i2z=1,i0z=1,iqs=1
REAL(DP),DIMENSION(N*N*N,N*N*N)::A
REAL(DP)::xmin=-1.5D0,ymin=-1.5D0,zmin=-1.5D0,h=0.15D0,mass=0.D0,temp=0.D0,w=0.D0,qg=0.D0,pvalue=1.D0,C
REAL(DP)::rhos,phis,phi0
REAL(DP),DIMENSION(N+2)::x,y,z
REAL(DP),DIMENSION(N*N*N)::phi,rhs,density
COMPLEX(DP),DIMENSION(4)::pvalues
INTEGER::unit=8,ierror
CHARACTER(len=10)::filename
CALL genmat(N,A)
Write(*,*) alpha(kappa,n1,rhoc)
!CALL ludecomp(N*N*N,A)

Write(*,*) 'Filling xyz'
Do i=1,N+2
  x(i)=xmin+h*(i-1)
  y(i)=ymin+h*(i-1)
  z(i)=zmin+h*(i-1)
End Do

Do i=1,N
  Do j=1,N
    Do k=1,N
       density((k-1)*N*N+(j-1)*N+i)=rho(x(i+1),y(j+1),z(k+1),qs,rhoc,kappa,n1)
    End Do
  End Do
End DO
Do noi=1,iter
Write(*,*) 'mass calculation'  
mass=0.D0  
Do i=1,N-1
  Do j=1,N-1
    Do k=1,N-1
      
      !mass=h*h*h*rho(x(i),y(j),z(k),qs,rhoc)+mass
      mass=h*h*h*density((k-1)*N*N+(j-1)*N+i)+mass
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
      rhs((k-1)*N*N+(j-1)*N+i)=4*pi*h*h*density((k-1)*N*N+(j-1)*N+i)-temp      
    End Do
  End Do
End Do

If(noi==1) Then
Do i=1,N
  Do j=1,N
    Do k=1,N
      
      phi((k-1)*N*N+(j-1)*N+i)=-mass/Sqrt(x(1)*x(1)+y(1)*y(1)+z(1)*z(1))
    End Do
  End Do
End Do
End If

CALL gaussiedel(N*N*N,A,phi,rhs,15)
!CALL gaussiedel1(N*N*N,phi,rhs,15,N)
!CALL lusolve(N*N*N,A,phi,rhs) 
w=omega(M,r)
Write(*,*) w
i0x=floor(1-xmin/h)
i1x=i0x+1
i2x=i0x-1
i0y=floor(1-ymin/h)
i1y=i0y+1
i2y=i0y-1
i0z=floor(1-zmin/h)
i1z=i0z+1
i2z=i0z-1
temp=(1/(2*h))*(phi((i0z-2)*N*N+(i0y-2)*N+i1x-1)-phi((i0z-2)*N*N+(i0y-2)*N+i2x-1))
qg=-temp/(w*w)
Write(*,*) qg,qg
iqs=floor(1+(qs-xmin)/h)
Write(*,*) iqs
rhos=density((i0z-2)*N*N+(i0y-2)*N+iqs)
phis=phi((i0z-2)*N*N+(i0y-2)*N+iqs-1)
phi0=phi((i0z-2)*N*N+(i0y-2)*N+i0x-1)
CALL SolvePolynomial(M*(qs**4)/(r**5),-M*(qs**3)/(r**4)-2*M*(3*qg*qs*qs/4-qs**3)*B(M,r,a1)&
	*Sqrt(1+B(M,r,a1)**2/(r**2))/(r**4),M*(qs**2)/(r**3)-M*(qg*qg-(qs-qg)**2)/(2*r**3)+phi0-phis,&
    0.D0,(1+n1)*kappa*(rhoc**(1/n1)-rhos**(1/n1)),code1,&
    pvalues(1),pvalues(2),pvalues(3),pvalues(4))
   Write(*,*) pvalues
!Do i=1,4
!  If(Abs(Dimag(pvalues(i)))<=0.01D0) Then
!    If(Dreal(pvalues(i))>=1.D0) Then
!      pvalue=dreal(pvalues(i))
!    End If
!  End If
!End DO
 
C=w*w*pvalue*pvalue*qg*qg/2-kappa*(n1+1)*(rhoc**(1/n1))-pvalue*pvalue*phi0
Write(*,*)rhos,pvalue,C,i0x
Do i=1,N
  Do j=1,N
    Do k=1,N
      If((Abs(x(i+1))>Abs(qs)).or.(Abs(y(j+1))>Abs(qs)).or.(Abs(z(k+1))>Abs(qs)))Then
       density((k-1)*N*N+(j-1)*N+i)=0.0D0
      Else 
      density((k-1)*N*N+(j-1)*N+i)=((w*w*pvalue*pvalue*((x(i+1)-qg)**2+z(k+1)**2)/2&
      -pvalue*pvalue*phi((k-1)*N*N+(j-1)*N+i)-pvalue*pvalue*(phitidal(M,r,pvalue,x(i+1),y(j+1),z(k+1))&
      +phimag(M,r,a1,qg,pvalue,x(i+1),y(j+1),z(k+1)))-C)/(kappa*(n1+1)))
      End If
      If(density((k-1)*N*N+(j-1)*N+i)<0.D0) Then
        !Write(*,*) i,j,k,density((k-1)*N*N+(j-1)*N+i)
        density((k-1)*N*N+(j-1)*N+i)=0.D0
      End If
      	density((k-1)*N*N+(j-1)*N+i)=density((k-1)*N*N+(j-1)*N+i)**n1
    End Do
  End Do
End Do 
End Do
Write(*,*) w*w/(pi*rhoc)
Write(*,*) w*w*(qs-qg)-(phi((i0z-2)*N*N+(i0y-2)*N+iqs)-phi((i0z-2)*N*N+(i0y-2)*N+iqs-2)&
+phitidal(M,r,pvalue,x(i0x+2),y(i0y+1),z(i0z+1))-phitidal(M,r,pvalue,x(i0x),y(i0y+1),z(i0z+1)))/(2*h)
filename='opoten.dat'
OPEN(UNIT=unit,FILE=filename,STATUS='replace',ACTION='WRITE',IOSTAT=ierror)
100 FORMAT(F15.5,1X,F15.5,1X,F15.5)
IF(ierror==0) Then
  !WRITE(unit,*) N
Do i=1,N
  Do j=1,N
    Write(unit,100) x(i+1),y(j+1),density(i+N*(j-1)+(i0z-1)*N*N)
  End Do
End DO
END IF
CLOSE(UNIT=unit)
Write(*,*) alpha(kappa,n1,rhoc)*lane(n1)
!Write(*,*) L(M,r,a1)
END PROGRAM roche