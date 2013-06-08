module outputfiles
use Density
use General

contains

subroutine WriteHam(TotH,TotEval,kvec,output)
!write this one so that it will just write the Hamiltonian file, either inimat or mat is controlled by Main. 
complex(kind=dp), dimension(:,:,:),intent(in) :: TotH
real(kind=dp), dimension(:,:),intent(in) :: TotEval
real(kind=dp), dimension(:),intent(in) :: kvec
character(len=*),intent(in) :: output
integer :: NQ = 0, NK = 0
NQ = size(TotH,1)
NK = size(kvec)

open(unit = 2,file = output,action="write")

write(2,*) "The final Eigenvectors are: "
do i = 1,nk,1
 write(2,*) "Nk = ", i, kvec(i)
!stored horizontally in H!! with lowest energy first!
 do j = 1, NQ, 1
!  Evec = TotH(:,j,i)
  write(2,*) TotH(j,:,i)  ! this is the order that eigenvectors are stored in H: by row!
 enddo
! write(1,'(7F12.8)') TotH(:,:,i)
 write(2,*) "With Eigenvalues:" 
 write(2,*) TotEval(:,i)
enddo

close(2)

end subroutine WriteHam

subroutine WriteDensity(TotH,nelect,cycles,L,output)
complex(kind=dp), dimension(:,:,:),intent(in) :: TotH
integer, intent(in) :: nelect, cycles, L
integer :: NQ = 0, M=0
real(kind=dp), dimension(:), allocatable :: Rho
real(kind=dp) :: totdens = 0
character(len=*),intent(in) :: output
NQ = size(TotH,1)
M = (NQ-1)/2


if(cycles == 1) then
 open(unit = 2,file = output,action="write")
else if (cycles > 1) then
 open(unit = 2,file = output,action="write",position="append")
endif

write(2,*) "At step ", cycles, ", density is:"

allocate(Rho(1:NQ))

Rho = 0
totdens = 0
!calculate and write out the density!
call NaiveDensity(TotH,Rho,nelect)
do i = 1,NQ,1
 write(2,*) qmap(i,L,M), Rho(i)
 totdens = totdens + Rho(i)
end do
!if there is a better way to write newlines without using a formatted print
!statement, I don't know what it is...
 write(2,*)
 write(2,*) totdens
 write(2,*)

!open .dens file, write to end of file??
! if this, then will have to handle wanting to create a new file... do that in main, or 
! can I check somehow? would still have to tell this function what step the main prog. was on..
! perhaps just send signal if cycles = 1, and then create the file new; otherwise, open from the end only. 
!or change "output" in Main, and write a seperate file for each step?

deallocate(Rho)

close(2)

end subroutine WriteDensity

subroutine WriteBands(TotEval,kvec,output)
real(kind=dp), dimension(:,:),intent(in) :: TotEval
real(kind=dp), dimension(:),intent(in) :: kvec
character(len=*),intent(in) :: output
integer :: NQ = 0, nk = 0

!write eigenvectors! generally only do this at the end, since those are the bands we want. Output could be 
!controlled in Main if we want bands before SCF is finished

NQ = size(TotEval,1)
nk = size(kvec)

open(unit = 2,file = output,action="write")

!should throw an error here if size(kvec) != size(TotEval)

!here, we can write to file
do i = 1,NQ,1
 write(2,*)
 write(2,*)
 write(2,*) "#Band No. ", i
 do j = 1,nk,1
  write(2,*) kvec(j), TotEval(i,j)
 end do
end do

close(2)

end subroutine WriteBands


end module outputfiles
