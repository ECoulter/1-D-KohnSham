program maintestImap
use General
use POSparser
integer :: i = 1,j=1,k=1,M = 3,NQ=0,L=4
type(ParamList) :: Params

call BuildIntList("M",defint=M,list=Params)
call BuildIntList("L",defint=L,list=Params)

NQ = 2*M + 1

do i = 1, (2*M + 1),1
print *, i, " goes to ", imap(i,M)
print *, i, " goes to qvec ", qmap(i,L,M)
print *, i, " goes to kvec ", kmap(i,L,M)
print *, i, " goes to xvec ", xmap(i,L,M)
end do

end program maintestImap
