module General

!this module contains basic things needed by most
!other modules

integer, parameter :: dp = selected_real_kind(8)
!double precision parameters
real(dp), parameter :: pi = 3.141592653589793238462
real(dp), parameter :: hbar = 6.58211928E-16
real(dp), parameter :: c = 3.0E8
real(dp), parameter :: me = 5.11E5
real(dp), parameter :: ang = 1E-10
real(dp), parameter :: Enprefix = (c**2*hbar**2)/(2*me*ang**2)

complex(dp), parameter :: Im = (0,1.00)

contains

!this is separate from qmap because it is better to compare 
!integers than real numbers, as in the delta f'ns in the
!definition of V


integer function imap(i,M)
IMPLICIT NONE
integer, intent(in) :: i, M
!this mapping does NOT include the endpoints at +/- pi/a
!Should it???
if(i >= 1 .AND. i <= 2*M+1) then
!if(i >= 1 .AND. i <= M) then
! imap = -(M-i+1)
!else if (i == M+1) then
! imap = 0
!else if (i > M+1 .AND. i <= 2*M+1) then
 imap = (i-(M+1))
else
 write(1,*) "Bad i in IMAP!"
 stop
endif

end function imap

real(dp) function qmap(i,L,M)
integer, intent(in) :: i, M, L

qmap = (2*imap(i,M)*pi)/L

end function qmap

real(dp) function xmap(i,L,M)
integer, intent(in) :: i, M, L

xmap = (L)*(dble(imap(i,M))/dble(2*M+1))

end function xmap

!this must be different from qmap because the mesh of k-points
!exists on a different scale (inside the 1st BZ) than the
!plane waves
real(dp) function kmap(i,L,n)
integer, intent(in) :: i, n, L
real(dp), parameter :: pi = 3.141592653589793238462
integer :: Nk = 0
Nk = (2*n) + 1
!this mapping includes the endpoints at +/- pi/a
!perhaps there is a better way...
kmap = pi*imap(i,n)/(L*n)
end function kmap


logical function DeltaqqprimeGee(iq,iqprime,igee,M)
integer, intent(in) :: iq, iqprime, igee, M
integer :: NQ
NQ = 2*M+1
DeltaqqprimeGee = .false.

if(imap(iq,M) - imap(iqprime,M) .eq. imap(igee,M)) then
  DeltaqqprimeGee = .true.
else if(imap(iq,M) - imap(iqprime,M) .eq. imap(igee,M) - NQ) then
  DeltaqqprimeGee = .true.
else if(imap(iq,M) - imap(iqprime,M) .eq. NQ + imap(igee,M)) then
  DeltaqqprimeGee = .true.
else
endif

end function Deltaqqprimegee


complex(dp) function FourierXtoQ(func,iq,L,M)
IMPLICIT NONE
complex(dp), intent(in), dimension(:) :: func
complex(dp) :: a = (0,0), b = (0,0)
real(dp) :: h = 0
integer, intent(in) :: M,iq,L
integer :: NQ = 0,j = 0
NQ = 2*M + 1
!NQ = 2*M 

FourierXtoQ = 0

if(NQ .NE. size(func)) then
 print *, "DIFFERENT SIZES IN FUNCTION!!!"
 stop
endif

!plain sum
do j = 1,NQ
!do j = -NQ/2,NQ/2,1
 FourierXtoQ = FourierXtoQ + exp(Im*qmap(iq,L,M)*xmap(j,L,M))* func(j)
enddo

end function FourierXtoQ

function FourierXtoQvec(func,L,M)
IMPLICIT NONE
complex(dp), intent(in), dimension(:) :: func
complex(dp), dimension(size(func)) :: FourierXtoQvec
integer, intent(in) :: M,L
integer :: NQ = 0,j = 0, i=0
NQ = 2*M + 1
!NQ = 2*M 

FourierXtoQvec = 0

if(NQ .NE. size(func)) then
 print *, "DIFFERENT SIZES IN FUNCTION!!!"
 stop
endif

!plain sum
do i = 1, NQ
 do j = 1,NQ
!do j = -NQ/2,NQ/2,1
  FourierXtoQvec(i) = FourierXtoQvec(i) + exp(Im*qmap(i,L,M)*xmap(j,L,M))* func(j)
 enddo
enddo

end function FourierXtoQvec


complex(dp) function FourierQtoX(func,ix,L,M)
IMPLICIT NONE
complex(dp), intent(in), dimension(:) :: func
complex(dp) :: area = (0,0),b=(0,0),a=(0,0)
real(dp) :: h = 0
integer, intent(in) :: M,ix,L
integer :: NQ = 0,j = 0
NQ = 2*M + 1
!NQ = 2*M 

FourierQtoX = 0

if(NQ .NE. size(func)) then
 print *, "DIFFERENT SIZES IN FUNCTION!!!"
 stop
endif

!plain sum
do j = 1,NQ
 FourierQtoX = FourierQtoX + (1.0/(dble(NQ)))*exp(-Im*qmap(j,L,M)*xmap(ix,L,M))*func(j)
enddo

end function FourierQtoX

end module General
