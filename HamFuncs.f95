module Hamfuncs
use General

contains

!this is separate from qmap because it is better to compare 
!integers than real numbers, as in the delta f'ns in the
!definition of V


real(dp) function TotEnergy(nelect,TotEval,nk,NQ)
IMPLICIT NONE
integer, intent(in) :: nelect, nk,nq
real(dp), intent(in), dimension(1:NQ,1:NK) :: TotEval
integer :: i = 0, j = 0

TotEnergy = 0
do i = 1,nelect
 do j = 1, NK
  TotEnergy = TotEnergy + TotEval(i,j)/NK
 end do
end do

end function TotEnergy

!this subroutine defines the kinetic energy matrix!
subroutine Tdef(T,a,M,k)
integer, intent(in) :: M,a
!the kinetic energy matrix in plane-wave space
complex(dp), intent(inout), dimension(2*M+1,2*M+1) :: T
real(dp), intent(in) :: k
integer :: j = 0
 T = 0
 do j = 1, (2*M+1), 1
  T(j,j) = (qmap(j,a,M) + k)**2
 end do
end subroutine Tdef

subroutine Vdef(V, L, M, k, vtype, V_0)
complex(dp), intent(inout), dimension(2*M+1,2*M+1) :: V
complex(dp), dimension(2*M+1) :: Vvec
complex(dp), dimension(2*M+1) :: cVvec
real(dp), intent(in) :: k, V_0
real(dp), parameter :: pi = 3.141592653589793238462
integer, intent(in) :: M, vtype, L
character :: PotFile*5 = "INPOT"
integer :: j = 0, i = 0,tmp = 0, NQ = 0

NQ = 2*M + 1

select case (vtype)
 case(1)
!kronig-penney model V = (V_0)
  V = V_0
 case(2)
!kronig-penney model V = 1/(V_0)
  V = 1/real(V_0)
 case(3)
!this will be a cosine potential,cos(2*pi/a*x)
!so V = V_0*pi*(delta(q-2*pi/a) + delta(q+2*pi/a))
  V = 0
!now, set the elements at V(q-q' = 2*pi/a) = pi*V_0
  do i = 1, NQ, 1
   do j = 1, NQ, 1
    if(abs(i - j) == 2) then
     V(i,j) = pi*V_0
    endif
   enddo
  enddo 
 case(4)
!this will be a cosine potential,cos(2*pi/a*x)
!so V = V_0*pi*(delta(q-2*pi/a) + delta(q+2*pi/a))
  V = 0
!now, set the elements at V(q-q' = 2*pi/a) = pi*V_0
  do i = 1, NQ, 1
   do j = 1, NQ, 1
    if(abs(i - j) == 2) then
     V(i,j) = pi*(1/V_0)
    endif
   enddo
  enddo 
 case(99)
!this needs a lot of work...
!error checking on the file, etc.
  open(unit=10,file=PotFile,action="read")
   read(10,*) tmp
  if(tmp == NQ) then
   do i = 1,NQ,1
    read(10,*) cVvec(i)
   enddo
  else
   print *, "Bad NQ in INPOT"
   close(10)
   stop
  endif
  close(10)
  cVvec = FourierXtoQvec(cVvec,L,M)
  call compVqvectoMat(V,cVvec,NQ)
 case default
  V = 0
end select
end subroutine Vdef

subroutine Vupdate(V,Rho,U_0,mixparam,vtype)
IMPLICIT NONE
integer, intent(in) :: vtype
real(dp), intent(in) :: mixparam,U_0
complex(dp), intent(inout), dimension(:,:) :: V
complex(dp), dimension(:,:),allocatable :: Vtempmat
complex(dp), dimension(:),allocatable :: Vtemp
real(dp), dimension(:),intent(in) :: Rho
integer :: i, j, k, NQ

NQ = size(Rho)

allocate(Vtemp(1:NQ))
allocate(Vtempmat(1:NQ,1:NQ))

Vtemp = 0
Vtempmat = 0

select case(vtype)
!this case is the plain U(rho) = U_0 * Rho**2
case(101:)

do j = 1, NQ, 1
 Vtemp(j) = Vtemp(j) + U_0*Rho(j)
end do

call compVqvectoMat(Vtempmat,Vtemp,NQ)
V = (1-mixparam)*V + mixparam*Vtempmat

case default
 V = V
end select

deallocate(Vtemp)
deallocate(Vtempmat)

end subroutine Vupdate

subroutine VqvectoMat(Vmat,Vvec,NQ)
IMPLICIT NONE
integer, parameter :: dp = selected_real_kind(8)
integer, intent(in) :: NQ
real(dp), dimension(1:NQ,1:NQ),intent(inout) :: Vmat
real(dp), dimension(1:NQ),intent(in) :: Vvec
real(dp), parameter :: tol = 0.0000001
integer :: i = 0, j = 0, k = 0, M = 0
M = (NQ-1)/2

do i = 1, NQ, 1
 do j = 1, NQ, 1
  do k = 1, NQ, 1
!this is the delta((q-q'),g)
!   if(imap(i,M) - imap(j,M) == imap(k,M)) then
   if(DeltaqqprimeGee(i,j,k,M)) then
    Vmat(i,j) = Vvec(k)
   end if 
  end do
 end do
end do

end subroutine VqvectoMat

subroutine compVqvectoMat(Vmat,Vvec,NQ)
IMPLICIT NONE
integer, parameter :: dp = selected_real_kind(8)
integer, intent(in) :: NQ
complex(dp), dimension(1:NQ,1:NQ),intent(inout) :: Vmat
complex(dp), dimension(1:NQ),intent(in) :: Vvec
integer :: i = 0, j = 0, k = 0, M = 0
M = (NQ-1)/2

do i = 1, NQ, 1
 do j = 1, NQ, 1
  do k = 1, NQ, 1
!this is the delta((q-q'),g)
!   if(imap(i,M) - imap(j,M) == imap(k,M)) then
   if(DeltaqqprimeGee(i,j,k,M)) then
    Vmat(i,j) = Vvec(k)
   end if 
  end do
 end do
end do

end subroutine compVqvectoMat

subroutine Vhartree(V,Rho,L)
IMPLICIT NONE
integer, intent(in) :: L
complex(dp), intent(inout), dimension(:,:) :: V
real(dp), dimension(:),intent(in) :: Rho
complex(dp), dimension(:),allocatable :: Vtemp
integer :: i, j, k, M, NQ

NQ = size(Rho)
M = (NQ-1)/2

allocate(Vtemp(1:NQ))

Vtemp = 0
V = 0

do j = 1, NQ, 1

 if(imap(j,M) .NE. 0) then
  Vtemp(j) = Vtemp(j) + (2*pi*Rho(j))/(qmap(j,L,M)**2)
 else
  Vtemp(j) = 0
 endif
 
end do

call compVqvectoMat(V,Vtemp,NQ)

deallocate(Vtemp)

end subroutine Vhartree

subroutine Vprint(vtype)
integer, intent(in) :: vtype

select case(vtype)
 case(1)
  write(1,'(A)') "Using Vtype = 1: constant potential = V_0 "
  write(1,'(A)')
 case(2)
  write(1,'(A)') "Using Vtype = 2: constant potential = 1/(V_0) "
  write(1,'(A)')
 case(3)
  write(1,'(A)') "Using Vtype = 3: cosine potentisl V = pi*V_0*cos(2*pi*x/L)"
  write(1,'(A)')
 case(4)
  write(1,'(A)') "Using Vtype = 4: cosine potentisl V = (pi/V_0)*cos(2*pi*x/L)"
  write(1,'(A)')
 case(199)
  write(1,'(A)') "Using Vtype = 99: Loaded from file INPOT"
  write(1,'(A)') "Also including scf update with V(q) = U_0*U(q)"
  write(1,'(A)')
 case(101)
  write(1,'(A)') "Using Vtype = 101: constant potential = V_0 "
  write(1,'(A)') "Also including scf update with V(q) = U_0*U(q)"
  write(1,'(A)')
 case(102)
  write(1,'(A)') "Using Vtype = 102: constant potential = 1/(V_0) "
  write(1,'(A)') "Also including scf update with V(q) = U_0*U(q)"
  write(1,'(A)')
 case(103)
  write(1,'(A)') "Using Vtype = 103: cosine potentisl V = pi*V_0*cos(2*pi*x/L)"
  write(1,'(A)') "Also including scf update with V(q) = U_0*U(q)"
  write(1,'(A)')
 case default
  write(1,'(A)') "Using Default Potential: V = 0"
  write(1,'(A)')
end select

end subroutine Vprint

end module Hamfuncs
