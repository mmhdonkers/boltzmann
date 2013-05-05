module model
  implicit none
  
  private calcaverdens, calcavervel
  public equil

contains

  subroutine newequil(Lx, Ly, Nvel, dens)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)
  end subroutine

  real(8) function equil(Lx, Ly, Nvel, dens) result(eqdens(Lx, Ly, Nvel)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: i, j, k, velx, vely
    real(8) :: avervel(Lx, Ly, 2), averdens(Lx, Ly)
    real(8) :: sqvel(Lx, Ly), normvel(Lx, Ly)

    velx = [0,1,1,0,-1,-1,-1,0,1]
    vely = [0,0,1,1,1,0,-1,-1,-1]
    avervel = call calcavervel(Lx, Ly, Nvel, dens) 
    averdens = call calcaverdens(Lx, Ly, Nvel, dens)   

    sqvel = avervel(:, :, 1)**2 + avervel(:, :, 2)**2
    eqdens(:, :, 1) = 4 / 9d0 * averdens * (1 - 1.5d0 * sqvel)

    normvel = avervel(:, :, 1) * velx + avervel(:, :, 2) * vely
    do k = 2, Nvel
      if (modulo(k,2)==0) then
        eqdens(:, :, k) = averdens / 9d0 * &
              (1 + 3 * normvel + 4.5d0 * normvel**2 - 1.5d0 * sqvel)
      else if (modulo(k,2)==1) then
        eqdens(:, :, k) = averdens / 36d0 * &
              (1 + 3 * normvel + 4.5d0 * normvel**2 - 1.5d0 * sqvel)
    end do
  end function

  real(8) function calcaverdens(Lx, Ly, Nvel, dens) result(averdens(Lx, Ly))
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: i, j

    do i = 1, Lx
      do j = 1, Ly
        averdens(i, j) = sum(dens(i, j, :)) / (Nvel * 1d0)
      end do
    end do
  end function

  real(8) function calcavervel(Lx, Ly, Nvel, dens) result(avervel(Lx, Ly, 2))
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: i, j, velx(Nvel), vely(Nvel)
    real(8) :: averdens(Lx, Ly)

    velx = [0,1,1,0,-1,-1,-1,0,1]
    vely = [0,0,1,1,1,0,-1,-1,-1]
    do i = 1, Lx
      do j = 1, Ly
        avervel(i, j, 1) = sum(velx * dens(i, j, :))
        avervel(i, j, 2) = sum(vely * dens(i, j, :))
      end do
    end do
    
    averdens = call calcaverdens(Lx, Ly, Nvel, dens)
    avervel(:, :, 1) = avervel(:, :, 1) / averdens
    avervel(:, :, 2) = avervel(:, :, 2) / averdens
  end function
end module
