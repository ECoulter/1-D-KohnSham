program maintestImap
use General
use POSparser
integer :: i = 1,j=1,k=1,M = 3,NQ=0,L=4
complex(dp), dimension(:), allocatable :: Potentialx 
complex(dp), dimension(:), allocatable :: qPotential 
type(ParamList) :: Params

call BuildIntList("M",defint=M,list=Params)
call BuildIntList("L",defint=L,list=Params)

NQ = 2*M + 1
!NQ = 2*M 

allocate(Potentialx(1:NQ))
allocate(qPotential(1:NQ))


!do i = 1,NQ,1
! potentialx(i) = exp(im*pi*2*i/(8))
!enddo
!qpotential = 10
potentialx = (0,0)
potentialx(M+1) = (10,0)
!qpotential(M+1) = 10
!do i = 1,NQ,1
!do i = 1,NQ
! potentialx(i) = sin(xmap(i,L,M))
! qpotential = sin(L/2.0 * (i/M))
!end do
!potentialx(M) = 10
!potentialx(M+1) = 10
!potentialx(M-1) = 10

!do i = 1,NQ,1
! print *, "V(",xmap(i,L,M),") is : ",potentialx(i)
!enddo

!do i = 1,NQ,1
do i = 1,NQ
! print *, "V(",qmap(i,L,M),") is : ",qpotential(i)
 print *, "V(",xmap(i,L,M),") is : ",potentialx(i)
enddo

print *, "Now, go to V(q):"
print *,

do i = 1,NQ,1
qPotential(i) = FourierXtoQ(Potentialx,i,L,M)
 print *, "V(",qmap(i,L,M),") is : ",qPotential(i)
end do

print *, "Now, go back to V(x):"
print *,

do i = 1,NQ,1
 Potentialx(i) = FourierQtoX(qPotential,i,L,M)
 print *, "V(",xmap(i,L,M),") is : ",Potentialx(i)
end do


deallocate(potentialx)
deallocate(qpotential)

end program maintestImap
