!WRITTEN BY TATHAGATA KARMAKAR, DEPARTMENT OF PHYSICS,IIT KANPUR
SUBROUTINE poisson(N,del,x,y,z,f,rho)
IMPLICIT NONE
INTEGER, PARAMETER::DP=KIND(1.0D0),SP=KIND(1.0)

INTEGER,INTENT(IN) :: N
REAL(DP),INTENT(IN)::del
REAL(DP),DIMENSION(2*N+1,2*N+1,2*N+1),INTENT(IN) :: rho
REAL(DP),DIMENSION(2*N+1,2*N+1,2*N+1),INTENT(OUT) :: f
REAL(DP),DIMENSION(2*N+1),INTENT(IN)::x,y,z

INTEGER ::i=0,j=0,k=0,l=0
REAL(DP),PARAMETER ::pi=3.141592
REAL(DP)::diff=0,diff1=0,maxdiff=1,tol=0.001,diff2=1,mindiff=2,temp,mass=0

!REAL(DP),DIMENSION(N+1)::y


!OPEN(UNIT=6,FILE="density.dat",STATUS="OLD",ACTION="READ")
Do i=1,2*N+1
	Do j=1,2*N+1
    	Do k=1,2*N+1
        !READ(6,*)rho(i,j,k)
       	!WRITE(*,*) rho(i,j,k)
!        rho(i,j,k)=0.015
         mass=rho(i,j,k)*del*del*del+mass
        end do
    end do
end do

write(*,*) mass,x(1),x(2*N+1),del
Do j=1,2*N+1
  do k=1,2*N+1
    f(1,j,k)=-mass/sqrt(x(1)**2+y(j)**2+z(k)**2)
    f(2*N+1,j,k)=-mass/sqrt(x(2*N+1)**2+y(j)**2+z(k)**2)
    f(j,1,k)=-mass/sqrt(x(j)**2+y(1)**2+z(k)**2)
    f(j,2*N+1,k)=-mass/sqrt(x(j)**2+y(2*N+1)**2+z(k)**2)
    f(j,k,1)=-mass/sqrt(x(j)**2+y(k)**2+z(1)**2)
    f(j,k,2*N+1)=-mass/sqrt(x(j)**2+y(k)**2+z(2*N+1)**2)
  end do
end do

Do i=2,2*N
  Do j=2,2*N
    Do k=2,2*N
      f(i,j,k)=f(1,j,k)+(x(i)-x(1))*(f(2*N+1,j,k)-f(1,j,k))/(x(2*N+1)-x(1))
    end do
  end do
end do

do l=1,10
  do i=2,2*N
	do j=2,2*N
    	do k=2,2*N
     	 temp=(1/6.)*(f(i+1,j,k)+f(i-1,j,k)+f(i,j+1,k)+f(i,j-1,k)+f(i,j,k+1)&
         			+f(i,j,k-1)-(del**2)*4*pi*rho(i,j,k))
         diff=DABS((temp-f(i,j,k))/temp)
         f(i,j,k)=temp
         If (diff>diff1) then
            diff1=diff
         end if
         IF(diff<diff2) then
           diff2=diff
          end if 
    	end do
    end do
  end do
 
  if(diff1<maxdiff) then
    maxdiff=diff1
  end if
  if(diff2<mindiff) then
    mindiff=diff2
  end if
  write(*,*) maxdiff,"  ",mindiff,"  ",f(78,59,68)
  !if(maxdiff<tol)then
   !EXIT
  !end if 
end do    
!CLOSE(UNIT=6)
end SUBROUTINE poisson


PROGRAM roche
IMPLICIT NONE
INTEGER, PARAMETER::DP=KIND(1.0D0),SP=KIND(1.0),N=50
INTEGER ::i,j,k,l
REAL(DP),DIMENSION(2*N+1,2*N+1,2*N+1)::f=0,rho=0
REAL(DP),DIMENSION(2*N+1) x,y,z
REAL(DP):: qs=-.08,M=1.,a=0,r=6,R0=0,del=0,rhoc=0.015
REAL(DP):: xmin=0,ymin=0,zmin=0

del=-qs/40;r=6*M;R0=10*M
xmin=-N*del;ymin=-N*del;zmin=-N*del


Do i=1,2*N+1
	IF(i==N+1) THEN
    x(i)=0
    y(i)=0
    z(i)=0
    ELSE
      x(i)=xmin+del*(i-1)
      y(i)=ymin+del*(i-1)
      z(i)=zmin+del*(i-1)
    END IF
end do

Do i=1,2*N+1
	Do j=1,2*N+1
    	Do k=1,2*N+1
        !READ(6,*)rho(i,j,k)
       	!WRITE(*,*) rho(i,j,k)
         rho(i,j,k)=rhoc
        end do
    end do
end do

CALL poisson(N,del,x,y,z,f,rho)

END PROGRAM roche

