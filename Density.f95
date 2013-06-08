module Density
use General

contains

subroutine NaiveDensity(U,Rho,nelect)
integer, parameter :: dp = selected_real_kind(8)
complex(dp), intent(in), dimension(:,:,:) :: U
real(dp), intent(out), dimension(:) :: Rho
real(dp), dimension(1:nelect) :: magu
integer, intent(in) :: nelect
integer :: NQ = 0
integer :: NK = 0
NQ = size(U,1)
NK = size(U,3)

!THIS IS IMPORTANT!
magu = 0

do i = 1, nelect, 1
 do j = 1, NQ, 1
  do k = 1, NK, 1
   magu(i) = magu(i) + ((U(j,i,k))*conjg(U(j,i,k)))/NK
  end do
 end do
end do
do i = 1, nelect, 1
 magu(i) = sqrt(magu(i))
end do
!magu = 1
!this needs to sum over occupied eigenvectors ONLY
!make density allow TWO electrons per state..
if(mod(nelect,2) .eq. 0) then
 do i = 1, nelect/2, 1
 !this sum over q-waves of occupied eigenvectors
  do j = 1,NQ, 1
 !this sum over k-points...
   do k = 1, nk, 1
    Rho(j) = Rho(j) + 2*((U(i,j,k))*conjg(U(i,j,k)))/(magU(i)*NK)
   end do
  end do
 end do
else if (mod(nelect,2) .ne. 0) then
 do i = 1, (nelect+1)/2, 1
 !this sum over q-waves of occupied eigenvectors
  do j = 1,NQ, 1
 !this sum over k-points...
   do k = 1, nk, 1
    if(i < (nelect+1)/2) then
     Rho(j) = Rho(j) + 2*((U(i,j,k))*conjg(U(i,j,k)))/(magU(i)*NK)
    else
     Rho(j) = Rho(j) + ((U(i,j,k))*conjg(U(i,j,k)))/(magU(i)*NK)
    endif
   end do
  end do
 end do
endif
 

end subroutine NaiveDensity

!subroutine densityfromrecip(uprime, ureg, xvec, qvec, density)
!complex(kind=8), intent(in), dimension(:) :: uprime, ureg
!double precision, intent(in), dimension(:) :: xvec, qvec
!double precision, intent(out), dimension(:) :: density
!integer :: i = 0, j = 0, l=0
!
!do l = 1, NQ, 1
! do i = 1, NQ, 1
!  do j = 1, NQ, 1
!  density(l) = uprime(j)*ureg(i)*(cos((qvec(i)-qvec(j))*xvec(l))  )
!!since this expression should be real, is fine for now. 
!!when I finally have complex wave-functions, I will have to change the
!!variable types here.
!!           I*sin((qvec(i)-qvec(j))*xvec(l))
!  enddo
! enddo
!enddo
!end subroutine
!
!subroutine denscalc(density,H,xvec,qvec,NQ,nelect,nk)
!double precision, intent(in), dimension(:) :: density, xvec, qvec
!integer, intent(in) :: NQ, nelect, nk
!double precision, intent(in), dimension(:,:,:) :: TotH
!integer :: i = 0, k = 0
!
!do k = 1, nk, 1
! do i = 1, nelect
!  call density(TotH
! end do
!do
!
!end subroutine denscalc

end module density
