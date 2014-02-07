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
  integer,parameter :: Lx = 40, Ly = 80, Nvel = 9, t_final = 100000
  integer :: tt, j
  real(8),parameter :: tau = 0.9d0, rho = 1.4d0
  real(8) :: dens(Lx, Ly, Nvel), vel(Lx, Ly, 2)

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
    call timestep(Lx, Ly, Nvel, dens, tau)
    if (modulo(tt, 1000) == 0) then
      call plot_profile(Lx, Ly, calcavervel(Lx, Ly, Nvel, dens))
    end if
  end do

  vel = calcavervel(Lx, Ly, Nvel, dens)
  do j = 1, Ly
    WRITE(15,*) vel(1, j, 1)
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
    OPEN(UNIT=15,FILE="velocityprofile.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, file not opened properly------------"
    endif
end subroutine

subroutine closetextfiles
    CLOSE(UNIT=15)
end subroutine

end program boltzmann
