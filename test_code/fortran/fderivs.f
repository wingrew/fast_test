subroutine fderivs(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx,fy,fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1),-1:ex(2),-1:ex(3))   :: fh
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call symmetry_bd(2,ex,f,fh,SoA)

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  fx = ZEO
  fy = ZEO
  fz = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1


! for bam comparison
   if(i+2 <= imax .and. i-2 >= imin .and. &
      j+2 <= jmax .and. j-2 >= jmin .and. &
      k+2 <= kmax .and. k-2 >= kmin) then
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))
      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))
      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))
   elseif(i+1 <= imax .and. i-1 >= imin .and. &
          j+1 <= jmax .and. j-1 >= jmin .and. &
          k+1 <= kmax .and. k-1 >= kmin) then
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))
      fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))
      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))
   endif

  enddo
  enddo
  enddo

  return

  end subroutine fderivs