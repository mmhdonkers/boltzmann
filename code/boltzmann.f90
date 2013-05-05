! Name:      Chiel Donkers
! Course:    International Course on Computational Physics
! Project:   Lattice Botzmann
! 
!Program Summary
!
!  Input:    
!
!  Process:  
!
!  Output:   
!



program boltzmann

  use model
  use plot
 
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! INPUT: Row and column size
  integer,parameter :: Lx = 150, Ly = 50, Nvel = 9
  real(8) :: dens(Lx, Ly, Nvel), newdens(Lx, Ly, Nvel), tau = 1, &
              eqdens(LX, Ly, Nvel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Initiate plot and open files for writing !!
  call plot_init()
  call opentextfiles



!! Close text and plot commands !!
  call closetextfiles
  call plot_close()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine opentextfiles
    integer :: OPEN_STATUS
    OPEN(UNIT=15,FILE="metrop_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_temp file not opened properly------------"
    endif
end subroutine

subroutine closetextfiles
    CLOSE(UNIT=15)
end subroutine

end program boltzmann
