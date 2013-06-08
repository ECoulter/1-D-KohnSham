program maintestImap
use Hamfuncs
use POSparser
integer :: i = 1,j=1,M = 2,NQ=0,L=4,vtype = 99
complex(dp), dimension(:,:),allocatable :: V
real(dp) :: k = 0,V_0=0
type(ParamList) :: Params


call BuildIntList("M",defint=M,list=Params)
call BuildIntList("L",defint=L,list=Params)

NQ = 2*M + 1
allocate(V(1:NQ,1:NQ))

call Vdef(V,L,M,k,vtype,V_0)

do i = 1, NQ,1
 write(*,'(14F7.3)') V(:,i)
end do

deallocate(V)

end program maintestImap
