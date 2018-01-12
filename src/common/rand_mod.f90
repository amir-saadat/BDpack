!Random module
!December 5, 2017

module rand_mod

  use :: prcn_mod

  implicit none

contains

  ! Random numeber seeding (from H. C. Ottinger):
  subroutine ranils(iseed)

    integer,intent(in) :: iseed
    integer,parameter :: in=2147483563,ik=40014,iq=53668,ir=12211,ntab=32
    integer :: iv(ntab),idum,idum2,iy
    integer :: k,j

    common /ranbls/ idum,idum2,iy,iv

    ! Initial seeds for two random number generators
    idum=iseed+123456789
    idum2=idum

    ! Load the shuffle table (after 8 warm-ups)
    do 10 j=ntab+8,1,-1
     k=idum/iq
     idum=ik*(idum-k*iq)-k*ir
     if(idum < 0) idum=idum+in
     if(j <= ntab) iv(j)=idum
     10 continue
     iy=iv(1)
     return
  end subroutine ranils

  ! Uniform random number generator (from H. C. Ottinger):
  real(wp) function ranuls()

   integer,parameter :: in1=2147483563,ik1=40014,iq1=53668,ir1=12211,&
   in2=2147483399,ik2=40692,iq2=52774,ir2=3791 ,&
   ntab=32,inm1=in1-1,ndiv=1+inm1/ntab
   real(wp),parameter :: an=1./in1
   integer :: iv(ntab),idum,idum2,iy
   integer :: k,j

   common /ranbls/ idum,idum2,iy,iv

   ! Linear congruential generator 1
   k=idum/iq1
   idum=ik1*(idum-k*iq1)-k*ir1
   if(idum < 0._wp) idum=idum+in1

   ! Linear congruential generator 2
   k=idum2/iq2
   idum2=ik2*(idum2-k*iq2)-k*ir2
   if(idum2 < 0._wp) idum2=idum2+in2

   !Shuffling and subtracting
   j=1+iy/ndiv
   iy=iv(j)-idum2
   iv(j)=idum
   if(iy < 1) iy=iy+inm1
   ranuls=an*iy
   return

  end function ranuls

  ! Gaussian random number generator (from H. C. Ottinger):
  real(wp) function rangls()

   integer :: iflag
   real(wp) :: gauss2,x1,x2,xsq,aux

   save iflag,gauss2
   data iflag/0/

   if(iflag == 0) then
    10 continue

    ! pair of uniform random numbers in [-1,1]x[-1,1]
    x1=2*ranuls()-1
    x2=2*ranuls()-1

    ! if not in the unit circle, try again
    xsq=x1*x1+x2*x2
    if(xsq >= 1._wp .or. xsq == 0._wp) goto 10
    ! pair of gaussian random numbers; return one and
    ! save the other for next time
    aux=sqrt(-2*log(xsq)/xsq)
    rangls=x1*aux
    gauss2=x2*aux
    iflag=1
  else
    rangls=gauss2
    iflag=0
  endif
  return

  end function rangls

end module rand_mod
