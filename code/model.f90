module model
  implicit none
  
  private calcaverdens, equil, applybc
  public timestep, init, calcavervel

contains

  subroutine timestep(Lx, Ly, Nvel, dens, tau)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: tau
    real(8),intent(inout) :: dens(Lx, Ly, Nvel)

    integer :: j
    real(8) :: eqdens(Lx, Ly, Nvel), stepdens(Lx, Ly, Nvel), deltaP

    deltaP = 0.001d0

    stepdens = applybc(Lx, Ly, Nvel, dens)

    do j = 2, Ly - 1
      stepdens(:, j, 2) = stepdens(:, j, 2) + deltaP
      stepdens(:, j, 6) = stepdens(:, j, 6) - deltaP
    end do

    eqdens = equil(Lx, Ly, Nvel, stepdens)
    dens = (1 - 1 / tau) * dens + eqdens / tau
  end subroutine
  
  function init(Lx, Ly, Nvel, rho) result(eqdens)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: rho

    integer :: k
    real(8) :: eqdens(Lx, Ly, Nvel)

    eqdens(:, :, 1) = 4 / 9d0 * rho

    do k = 2, Nvel
      if (modulo(k,2)==0) then
        eqdens(:, :, k) = rho / 9d0
     else if (modulo(k,2)==1) then
        eqdens(:, :, k) = rho / 36d0
      end if
    end do
  end function

  function equil(Lx, Ly, Nvel, dens) result(eqdens)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: k, velx(Nvel), vely(Nvel)
    real(8) :: avervel(Lx, Ly, 2), averdens(Lx, Ly), eqdens(Lx, Ly, Nvel)
    real(8) :: sqvel(Lx, Ly), normvel(Lx, Ly)

    velx = [0,1,1,0,-1,-1,-1,0,1]
    vely = [0,0,1,1,1,0,-1,-1,-1]
    avervel = calcavervel(Lx, Ly, Nvel, dens) 
    averdens = calcaverdens(Lx, Ly, Nvel, dens)   

    sqvel = avervel(:, :, 1)**2 + avervel(:, :, 2)**2
    eqdens(:, :, 1) = 4 / 9d0 * averdens * (1 - 1.5d0 * sqvel)

    do k = 2, Nvel
      normvel(:, :) = avervel(:, :, 1) * velx(k) + avervel(:, :, 2) * vely(k)
      if (modulo(k,2)==0) then
        eqdens(:, :, k) = averdens / 9d0 * &
              (1 + 3 * normvel + 4.5d0 * normvel**2 - 1.5d0 * sqvel)
      else if (modulo(k,2)==1) then
        eqdens(:, :, k) = averdens / 36d0 * &
              (1 + 3 * normvel + 4.5d0 * normvel**2 - 1.5d0 * sqvel)
      end if
    end do
  end function

  function calcaverdens(Lx, Ly, Nvel, dens) result(averdens)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: i, j
    real(8) :: averdens(Lx, Ly)

    do i = 1, Lx
      do j = 1, Ly
        averdens(i, j) = sum(dens(i, j, :))
      end do
    end do
  end function

  function calcavervel(Lx, Ly, Nvel, dens) result(avervel)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: i, j, velx(Nvel), vely(Nvel)
    real(8) :: averdens(Lx, Ly), avervel(Lx, Ly, 2)

    velx = [0,1,1,0,-1,-1,-1,0,1]
    vely = [0,0,1,1,1,0,-1,-1,-1]
    do i = 1, Lx
      do j = 1, Ly
        avervel(i, j, 1) = sum(velx * dens(i, j, :))
        avervel(i, j, 2) = sum(vely * dens(i, j, :))
      end do
    end do
    
    averdens = calcaverdens(Lx, Ly, Nvel, dens)
    avervel(:, :, 1) = avervel(:, :, 1) / (averdens * 1d0)
    avervel(:, :, 2) = avervel(:, :, 2) / (averdens * 1d0)
  end function

  function applybc(Lx, Ly, Nvel, dens) result(newdens)
    integer,intent(in) :: Lx, Ly, Nvel
    real(8),intent(in) :: dens(Lx, Ly, Nvel)

    integer :: i, j, ii, jj
    real(8) :: newdens(Lx, Ly, Nvel)

    newdens = 0d0

!    do i = 1, Lx
!      newdens(i, 1, 4) = dens(i, 2, 8)
!      newdens(i, 1, 5) = dens(modulo(i - 1, Lx), 2, 9)
!      newdens(i, 1, 3) = dens(modulo(i + 1, Lx), 2, 7)
!      newdens(i, Ly, 8) = dens(1, Ly - 1, 4)
!      newdens(i, Ly, 9) = dens(modulo(i + 1, Lx), Ly - 1, 5)
!      newdens(i, Ly, 7) = dens(modulo(i - 1, Lx), Ly - 1, 3)
!    end do
    do i = 1, Lx
      do j = 1, Ly
        ii = modulo(i + 1, Lx)
        jj = modulo(j + 1, Ly)
        newdens(i, j, 1) = dens(i, j, 1)
        newdens(ii, j, 2) = dens(i, j, 2)
        newdens(ii, jj, 3) = dens(i, j, 3)
        newdens(i, jj, 4) = dens(i, j, 4)
        newdens(i, jj, 5) = dens(ii, j, 5)
        newdens(i, j, 6) = dens(ii, j, 6)
        newdens(i, j, 7) = dens(ii, jj, 7)
        newdens(i, j, 8) = dens(i, jj, 8)
        newdens(ii, j, 9) = dens(i, jj, 9)
!        do k = 2, Nvel
!          ii = modulo(i + Lx + velx(k) - 1, Lx) + 1
!          jj = j + vely(k)
!          if (jj < Ly .and. jj > 1) then
!            newdens(ii, jj, k) = dens(i, j, k)
!          end if
!        end do
      end do
    end do

    do i = 1, Lx
      newdens(i, 1, 3) = newdens(i, 1, 7)
      newdens(i, 1, 4) = newdens(i, 1, 8)
      newdens(i, 1, 5) = newdens(i, 1, 9)
      newdens(i, Ly, 7) = newdens(i, Ly, 3)
      newdens(i, Ly, 8) = newdens(i, Ly, 4)
      newdens(i, Ly, 9) = newdens(i, Ly, 5)
    end do
  end function
end module
