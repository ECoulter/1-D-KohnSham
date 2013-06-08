      program testParser
      use POSParser
      IMPLICIT NONE
      type(ParamList) :: Params
      integer :: testint = 0, testint2 = 1
      real(kind=selected_real_kind(8)) :: dubtest = 0.0000
      character :: testchar*32 = "THIS IS THE DEFAULT"

      call BuildIntList("testi1",defint=testint,list=Params)
      call BuildIntList("testi2",defint=testint2,list=Params)
      call BuildDbleList("dubtest",defdble=dubtest,list=Params)
      call BuildCharList("testc",defchar=testchar,list=Params)
      print *, "This is testchar:",testchar
      if(NeedHelp(Params)) GOTO 100

! make a module to recognize three types of input arguments:
! filenames, boolean, and parameter (int or real)
! needs two pieces of info : variable name and type
! just need a subroutine to run through cmdline, match up to the desired name, and set the right variable to the input value. 
! also, output which values are used; note which were left default?
! could also have this as an option?
! build module; when call "getparam", with a param name and default, include in list of parameters!
! PROBLEM WITH LENGTH OF CHARACTER VARIABLES!!!
100   end program testParser
