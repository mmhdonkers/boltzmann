module plot
  use plplot

  implicit none

  private numtostr
  public plot_init, plot_close, plot_spin

contains

  subroutine plot_init() 
    call plsdev("xcairo")

    call plscol0(0, 255, 255, 255)  ! white
    call plscol0(1, 255, 0, 0)      ! red
    call plscol0(2, 0, 255, 0)      ! green
    call plscol0(3, 0, 0, 255)      ! blue
    call plscol0(4, 255, 0, 255)    ! magenta
    call plscol0(5, 0, 255, 255)    ! cyan
    call plscol0(6, 255, 255, 0)    ! yellow
    call plscol0(7, 0, 0, 0)        ! black
    call plscol0(8, 255, 77, 0)     ! orange
    call plscol0(9, 128, 128, 128)  ! gray

    call plinit()
  end subroutine

  subroutine plot_close()
    call plspause(.false.)
    call plend()
  end subroutine

  subroutine plot_spin(spin, SIZE, temp)
    integer,intent(in) :: SIZE
    integer,intent(in) :: spin(0:SIZE-1,0:SIZE-1)
    real(8),intent(in) :: temp

    integer :: i, j

    call plcol0(7)
    call plenv(0d0, SIZE*1d0, 0d0, SIZE*1d0, 0, 0)
    call pllab("x", "y", "coupling constant: " // numtostr(1d0/temp))

    do i = 0, SIZE - 1
      do j = 0, SIZE - 1
        if (spin(i,j) .eq. -1) then
          call plcol0(1)
          call plpoin([i + 0.5d0], [j + 0.5d0], 31)
        else if (spin(i,j) .eq. 1) then
          call plcol0(2)
          call plpoin([i + 0.5d0], [j + 0.5d0], 30)
        end if
      end do
    end do

    call plspause(.false.)
  end subroutine

  subroutine plot_profile(Lx, Ly, vel)
    integer,intent(in) :: Lx, Ly
    real(8),intent(in) :: vel(Lx,  Ly, 2)

    integer :: i
    real(8) :: y(Ly)

    y = [(i, i = 1, Ly)]

    call plcol0(7)
    call plenv(0d0, maxval(vel(1, :, 1))*1.1, 1d0, Ly * 1d0, 0, 0)
    call pllab("v", "y", "velocity profile")
    
    call plcol0(1)
    call plline(vel(1, :, 1), y)

    call plspause(.false.)
  end subroutine

  character(len=25) function numtostr(num) result(str)
    real(8),intent(in) :: num

    write(str, '(g12.5)') num
  end function

end module
