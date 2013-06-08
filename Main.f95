program main
use POSparser
use HamFuncs
use CompDiagonal
use Density
use Outputfiles
use General
IMPLICIT NONE
!define basis with some function qmap
!define Kinetic Energy Operator with Tdef
!define Potential in VDef
!loop over K-points, and scf cycles
!write to EIGENVAL file
!ADD:
!    normalization of the output vectors? Already normalized
!    updating of the potential! done- naively
!    Some fourier functionality for creating potentials and back for density...
!    density plot - done, in q-space!
!    density of states uggh...


!define INTEGER parameters:
!this is defined in POSPARSER!
type(ParamList) :: Params
!for plane-wave expansion, number of +/-M summed from 0
integer :: M = 5
!Total # of k-points THIS MUST BE ODD FOR NOW!!!
integer :: nk = 51
!type of potential; have to build a 'select case' in Vdef for this
! 1 gives V = V_0, 2 gives V = 1/V_0, 3 gives V(x) = V_0*cos(2*pi/a*x)
! greater than 100 gives self-consistency with Vint = U_0*U(q)!
integer :: vtype = 1
!parameter for amplitude of potential
real(kind=dp) :: V_0 = 0
!parameter for amplitude of interaction potential
real(kind=dp) :: U_0 = 0
!parameter for mixing in self-consistency
real(kind=dp) :: mixparam = 1
!length of unit cell in ANGSTROMS
integer :: L = 4
!maximum # of scf cycles
integer :: maxcycles = 1
!number of "particles" for Density and Fermi Energy, etc. 
integer :: nelect = 4
!negative power of ten to check energy tolerance to
integer :: EnTolPower = 4

!for checking density calculation
real(dp) :: totdens = 0.0


!total # of plane waves; size of H, T, and V
!this cannot be an input parameter! (for safety)
integer :: NQ = 0
!Max EnCut
real(dp) :: Encut
!Total Energy storage
real(dp) :: TotEn = 0, TotEnPrev = 0, EnTol = 0.0001

!variables for counting
!some general ones
integer :: i = 0, j = 0, k = 0, n = 0
!number of scf cycles
integer :: cycles = 0

!VECTORS AND MATRICES
!vector for k-points
real(dp), allocatable, dimension(:) :: kvec 

!vector for q-points
real(dp), allocatable, dimension(:) :: qvec 

!vector for x-points
real(dp), allocatable, dimension(:) :: xvec 

!array for density
real(dp), allocatable, dimension(:) :: Rho

!for Hamiltonian, Kinetic E, Pot. E
complex(dp), allocatable, dimension(:,:) :: H, T, U
complex(dp), allocatable, dimension(:,:) :: Vext,Vint,Vhart
!complex(kind=8), allocatable, dimension(:,:) :: H, T, Vext,Vint
!for work matrix(in DIAG), Eigenvalues
real(dp), allocatable, dimension(:) :: Eval
!complex(dp), allocatable, dimension(:) :: Eval
!complex(kind=8), allocatable, dimension(:) :: Eval, Work

!For storing ALL output across all k-points
!solved Hamiltonians
complex(dp), allocatable, dimension(:,:,:) :: TotH
!complex(kind=8), allocatable, dimension(:,:,:) :: TotH
!eigenvalues (duh)
real(dp), allocatable, dimension(:,:) :: TotEval
!complex(dp), allocatable, dimension(:,:) :: TotEval
!complex(kind=8), allocatable, dimension(:,:) :: TotEval
!eigenvectors (just a check)
!double precision, allocatable, dimension(:) :: Evec

!booleans for loops!
logical :: unfinished = .true.

!output filename
character :: output*32 = "Default"
character :: hart_tag*32 = "on"


!build the list of parameters for input - Has to go after the variable def'ns
call BuildIntList("M",defint=M,list=Params)
call BuildIntList("NK",defint=nk,list=Params)
call BuildIntList("vtype",defint=vtype,list=Params)
call BuildDbleList("V_0",defdble=V_0,list=Params)
call BuildDbleList("U_0",defdble=U_0,list=Params)
call BuildDbleList("mixparam",defdble=mixparam,list=Params)
call BuildIntList("L",defint=L,list=Params)
call BuildIntList("maxsteps",defint=maxcycles,list=Params)
call BuildIntList("nelect",defint=nelect,list=Params)
call BuildIntList("EnTol",defint=EnTolPower,list=Params)
call BuildCharList("output",defchar=output,list=Params)
call BuildCharList("Hart-tag",defchar=hart_tag,list=Params)
if(NeedHelp(Params)) then
 stop
endif

EnTol = 10**(-1*dble(EnTolPower))
NQ = (2*M) + 1
EnCut = EnPrefix*(qmap(NQ,L,M))**2

!allocate space for the matrices!!
allocate(kvec(1:nk))
allocate(H(1:NQ,1:NQ))
allocate(U(1:NQ,1:NQ))
allocate(T(1:NQ,1:NQ))
allocate(Vint(1:NQ,1:NQ))
allocate(Vext(1:NQ,1:NQ))
allocate(Vhart(1:NQ,1:NQ))
allocate(Eval(1:NQ))
allocate(qvec(1:NQ))
allocate(xvec(1:NQ))
allocate(Rho(1:NQ))
allocate(TotH(1:NQ,1:NQ,1:nk))
allocate(TotEval(1:NQ,1:nk))

kvec = 0
H = 0
T = 0
Vint = 0
Vext = 0
Vhart = 0
Eval = 0
qvec = 0
xvec = 0
Rho = 0
TotH = 0
TotEval = 0


open(unit = 1,file = trim(output)//".out",action = "write")
open(unit = 5,file = trim(output)//".inimat",action = "write")
open(unit = 7,file = trim(output)//".vint",action = "write")

write(1,'(A,/,A,I3,/,A,I3/)') "Welcome to the 1-D Kohn-Sham Solver!","You are using ",nk," k-points, with potential type ", vtype
write(1,'(A,I3,/,A,F14.6,A/)') "The length of your unit cell is: ",L,"The cutoff in your plane wave expansion is: ", EnCut, "eV."
write(1,'(A,I3,A,I2,A,I2/)') "There are: ", NQ," plane waves in your expansion from 0 to +/- 2 *",M," * pi/",L
write(1,'(A,I3,A/)') "There will be ",maxcycles," maximum steps in the self-consistency cycle."
write(1,'(A,F14.6/)') "The mizing parameter for the update of V is:",mixparam
write(1,'(A,F14.6,A/)') "Energy mst be converged to ", EnTol, "eV."
write(1,'(A)') "The Hamiltonian at each step is in "//trim(output)//".inimat"
call Vprint(vtype)

do i = 1, nk, 1
 kvec(i) = kmap(i,L,(nk-1)/2)
end do
 write(1,'(A)') "The k-points used are: "
 write(1,'(3F14.6)') kvec

 write(1,'(A)') "The q-vectors used are: "
do i = 1, NQ, 1
 qvec(i) = qmap(i,L,M)
 xvec(i) = dble(i)*dble(L)/(dble(NQ))
 write(1,'(3F14.6)') qmap(i,L,M)
end do

cycles = 1

do while (unfinished)
 unfinished = .true.

!this updates self-consistent interaction potential!
 if(cycles .gt. 1) then
  Rho = 0
  call NaiveDensity(TotH,Rho,nelect)
  call Vupdate(Vint,Rho,U_0,mixparam,vtype)
  Vhart = 0
print *, hart_tag
  if(trim(hart_tag)=="on") call Vhartree(Vhart,Rho,L)

  write(7,'(A,I5,A)') "The Interaction Potential at step ",cycles," is:"
  write(7,'(7F12.8)') Vint(:,:)

 endif
 do i = 1, nk, 1

!this may need to change for non-local functionals...
  call Tdef(T,L,M,kvec(i))

  if(cycles .eq. 1) then
!define Potential Energy matrix
   if(vtype < 100) call Vdef(Vext,L,M,kvec(i),vtype,V_0)
   if(vtype > 100) call Vdef(Vext,L,M,kvec(i),vtype-100,V_0)
  endif

  H = T + Vext + Vint + Vhart

   if(cycles .eq. 1) then
    write(5,'(A)') "The Initial Hamiltonian is: "
    write(5,'(A,I5,F14.6)') "Nk = ", i, kvec(i)
    write(5,'(7F12.8)') H(:,:)
   else if(cycles .gt. 1) then
    write(5,'(A,I5,A)') "The Initial Hamiltonian at step ", cycles," is: "
    write(5,'(A,I5,F18.8)') "Nk = ", i, kvec(i)
    write(5,'(7F14.8)') H(:,:)
   endif

  Eval=0

  call HEigensystem(NQ,H,NQ,Eval,U,NQ,1)

  TotH(:,:,i) = U
  U = 0
  TotEval(:,i) = Eval

 end do

 TotEnPrev = TotEn
 TotEn = EnPrefix*TotEnergy(nelect,TotEval,NK,NQ)

 write(1,'(A,I5,A,F14.6)') "At step ", cycles, " the total energy is: ", TotEn

 call WriteDensity(TotH,nelect,cycles,L,trim(output)//".dens")

!check if done or converged!
 if(cycles == maxcycles) then
  unfinished = .false.
 else if (abs(TotEn - TotEnPrev) < EnTol) then
  write(1,'(A,F14.6)') "Previous Energy ", TotEnPrev
  write(1,'(A,F14.6,A,F14.6)') " matches current energy ",TotEn," within tolerance: ",EnTol
  unfinished = .false.
 endif

 cycles = cycles + 1

end do

!here, write out the Useful Output
call WriteHam(TotH,TotEval,kvec,trim(output)//".mat")
call WriteDensity(TotH,nelect,cycles,L,trim(output)//".dens")
call WriteBands(TotEval,kvec,trim(output)//".eval")
!call WriteDOS(addthislater)
write(1,'(A,I5,A,I5)') "Finished after: ", cycles, "steps out of ",maxcycles
write(1,'(A,F14.6)') "The final Total Energy is: ", TotEn
write(1,'(A,A)') "The final eigenvectors are stored in ", trim(output)//".mat"

!deallocate space for the matrices!!
deallocate(kvec)
deallocate(H)
deallocate(U)
deallocate(T)
deallocate(Vint)
deallocate(Vext)
deallocate(Eval)
deallocate(qvec)
deallocate(xvec)
deallocate(Rho) !
deallocate(TotH)
deallocate(TotEval)

close(1)
close(5)
close(7)
end program main
