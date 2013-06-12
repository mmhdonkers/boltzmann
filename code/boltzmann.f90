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
  integer,parameter :: Lx = 250, Ly = 100, Nvel = 9, t_final = 10000
  integer :: tt
  real(8),parameter :: tau = 10d0, rho = 1d0
  real(8) :: dens(Lx, Ly, Nvel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Initiate plot and open files for writing !!
  call plot_init()
  call opentextfiles
  
!! Initialize densities for further use !!
  dens = init(Lx, Ly, Nvel, rho) 

!! Main simulation routine !!
  do tt = 1, t_final
      call plot_profile(Lx, Ly, calcavervel(Lx, Ly, Nvel, dens))
    call timestep(Lx, Ly, Nvel, dens, tau)
!    if (modulo(tt, 10) == 0) then
!    end if
  end do

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
