program maintestImap
use General
use POSparser
integer :: i = 1,j=1,k=1,M = 3,NQ=0,L=4
type(ParamList) :: Params

call BuildIntList("M",defint=M,list=Params)
call BuildIntList("L",defint=L,list=Params)

NQ = 2*M + 1

do i = 1, NQ,1
 do j = 1, NQ,1
  do k = 1, NQ,1
   if(DeltaqqprimeGee(i,j,k,M)) then
    print *, "Delta (",imap(i,M),"-",imap(j,M),"=",imap(i,M)-imap(j,M), ") goes to ", imap(k,M)
   else
!    print *, "Delta (",imap(i,M),"-",imap(j,M),",",imap(k,M),") FAIL"
   endif
  end do
 end do
end do

end program maintestImap
