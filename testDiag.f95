program maintestImap
use CompDiagonal
use General
use POSparser
integer :: i = 1,j=1,M = 2,NQ=0,L=4,vtype = 99
complex(dp), dimension(:,:),allocatable :: H,U
real, allocatable, dimension(:) :: Eval
real(dp) :: k = 0,V_0=0
type(ParamList) :: Params


call BuildIntList("M",defint=M,list=Params)
call BuildIntList("L",defint=L,list=Params)

NQ = 2*M + 1
allocate(H(1:NQ,1:NQ))
allocate(U(1:NQ,1:NQ))
allocate(Eval(1:NQ))

U = 0
H = 0
Eval = 0

do i = 1,NQ
 do j = 1,NQ
  H(i,j) = 1*(-1)**j
 enddo
enddo


do i = 1, NQ,1
 write(*,'(10F7.3)') H(:,i)
end do

call HEigensystem(NQ,H,NQ,Eval,U,NQ,1)

do i = 1, NQ,1
 write(*,'(10F7.3)') U(:,i)
end do

print *, "With Eigenvalues:"

write(*,'(/10F7.3)') Eval(:)

deallocate(H)
deallocate(U)
deallocate(Eval)

end program maintestImap
