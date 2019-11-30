module errcalc
  use const
  use omp_lib
  use BC
  use kernel,          only: get_dim,&
                             get_tasktype
  use neighboursearch, only: getNeibListL1

  implicit none

  public :: err_T0sxsyet, err_infplate, err_sinxet,&
            err_diff_laplace, err_diff_graddiv

  private

contains

  subroutine err_T0sxsyet(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), num(n), t
    real, intent(out)   :: err(n)

    integer             :: i
    real                :: exact

    print*, 'Not ready to NBS will divide to random number'

    !$OMP PARALLEL
    !$OMP DO PRIVATE(exact)
    do i=1,n
      exact = sin(pi*(pos(1,i)+1.)/2.) * sin(pi*(pos(2,i)+1.)/2.) * exp(-2.*(pi/2.)**2 * 0.1 * t)
      err(i) = abs(exact - num(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine err_T0sxsyet

  subroutine err_sinxet(ptype, num, t, err, count)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: num(:)
    real, allocatable, intent(inout) :: err(:)
    integer, intent(inout)           :: count
    real, intent(in)                 :: t

    integer             :: i, n
    real                :: exact

    print*, 'Not ready to NBS will divide to random number'

    n = size(ptype)
    count = 0
    err(1:n) = 0.
    !$omp parallel do default(none) &
    !$omp shared(n,num,err,t,ptype) &
    !$omp private(exact, i)&
    !$omp reduction(+:count)
    do i=1,n
      ! if (ptype(i) /= 0) then
      !   exact = tsin(i) * exp(-pi**2 * t)
      !   err(i) = (exact - num(i))**2
      !   count = count + 1
      ! end if
      err(i) = -100000000.
    end do
    !$omp end parallel do
  end subroutine err_sinxet

  subroutine err_diff_laplace(ptype, x, num, err, count)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: x(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)
    integer, intent(out)             :: count

    integer, allocatable :: nlista(:)
    integer              :: i, j, dim
    real                 :: exact(1:3)

    call get_dim(dim)
    call getNeibListL1(nlista)
    err(:) = 0.
    count = 0
    exact(:) = 0.
    !$omp parallel do default(none) &
    !$omp shared(x, ptype, num, err, dim, nlista) &
    !$omp private(exact, i, j) &
    !$omp reduction(+:count)
    do j = 1,size(nlista)
      i = nlista(j)
      if (ptype(i) == 1) then
        exact(:) = 0.
        if ( dim == 1 ) then
          ! exact(1) = 2*Cos(x(1,i)) - x(1,i)*Sin(x(1,i))
          ! sin
          exact(1) = -sin(x(1,i))
        elseif ( dim == 2 ) then
          ! exact(1) = -x(2,i)*Sin(x(1,i))
          ! exact(2) = -x(1,i)*Sin(x(2,i))
          ! sin
          exact(1) = -sin(x(1,i))
          exact(2) = -sin(x(2,i))
        elseif ( dim == 3 ) then
          ! exact(1) = -(x(2,i)*Sin(x(1,i)))
          ! exact(2) = -(x(3,i)*Sin(x(2,i)))
          ! exact(3) = -(x(1,i)*Sin(x(3,i)))
          ! sin
          exact(1) = -sin(x(1,i))
          exact(2) = -sin(x(2,i))
          exact(3) = -sin(x(3,i))
        end if
        err(i) = dot_product(exact(:) - num(:,i),exact(:) - num(:,i))/dim
        count = count + 1
      endif
    end do
    !$omp end parallel do
  end subroutine err_diff_laplace

  subroutine err_diff_graddiv(ptype, x, num, err, count)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: x(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)
    integer, intent(out)             :: count

    integer             :: n, i, dim, la
    real                :: exact(1:3), xk(3)
    integer, allocatable :: nlista(:)

    call get_dim(dim)
    n = size(ptype)
    count = 0
    err(:) = 0.

    call getNeibListL1(nlista)

    !$omp parallel do default(none) &
    !$omp shared(n,ptype, x,num,err,dim, nlista) &
    !$omp private(exact, i,xk, la) &
    !$omp reduction(+:count)
    do la = 1,size(nlista)
      i = nlista(la)
      if (ptype(i) == 1) then
        ! print*, i, num(:,i)
        exact(:) = 0.
        if (dim == 1) then
          ! exact(1) = 0
          ! exact(1) = 2*Cos(x(1,i)) - (x(1,i))*Sin(x(1,i))
          ! sin
          exact(1) = -sin(x(1,i))
          ! grad only
          ! exact(1) = cos(x(1,i))
        end if
        if (dim == 2) then
          ! exact(1) = 1
          ! exact(2) = 1
          ! exact(1) = Cos(x(2,i)) - x(2,i)*Sin(x(1,i))
          ! exact(2) = Cos(x(1,i)) - x(1,i)*Sin(x(2,i))
          ! sin
          exact(1) = -sin(x(1,i))
          exact(2) = -sin(x(2,i))
          ! grad only
          ! exact(1) = cos(x(1,i))
          ! exact(2) = cos(x(2,i))
        end if
        if (dim == 3) then
          ! exact(1) = x(2,i) + x(3,i)
          ! exact(2) = x(1,i) + x(3,i)
          ! exact(3) = x(1,i) + x(2,i)
          ! exact(1) = Cos(x(3,i)) - (x(2,i)*Sin(x(1,i)))
          ! exact(2) = Cos(x(1,i)) - (x(3,i)*Sin(x(2,i)))
          ! exact(3) = Cos(x(2,i)) - (x(1,i)*Sin(x(3,i)))
          ! sin
          exact(1) = -sin(x(1,i))
          exact(2) = -sin(x(2,i))
          exact(3) = -sin(x(3,i))
          ! grad only
          ! exact(1) = cos(x(1,i))
          ! exact(2) = cos(x(2,i))
          ! exact(3) = cos(x(3,i))
        end if
        ! print*, exact
        ! print*, num(:,i)
        ! print*, '----------'
        ! read*
        err(i) = dot_product(exact(:)-num(:,i),exact(:)-num(:,i))/dim
        count = count + 1
      end if
    end do
    !$omp end parallel do
  end subroutine err_diff_graddiv

  subroutine err_infplate(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), num(n), t
    real, intent(inout) :: err

    integer :: i
    real :: tl, tr, tc, al, ar, ttmp, exact, xm, kl, kr, rhol, rhor, cvl, cvr

    print*, 'Not ready to NBS will divide to random number'

    kl = 1.
    kr = 1.
    rhol = 1.
    rhor = 1.
    cvl = 1.
    cvr = 1.
    tl = 0.
    tr = 1.
    xm = 0.

    al = kl / rhol / cvl
    ar = kr / rhor / cvr

    tc = (tr - tl) * (kr / sqrt(ar)) / (kr / sqrt(ar) + kl / sqrt(al))

    err = 0
    !$omp parallel do default(none) &
    !$omp shared(pos,n,xm,al,ar,kl,kr,tc,tl,num,t) &
    !$omp private(exact, ttmp, i) &
    !$omp reduction(+:err)
    do i=1,n
      if (pos(1,i) < xm) then
        ttmp = erfc((xm-pos(1,i))/(2 * sqrt(al*t)))
      else
        ttmp = 1 + (kl/kr)*sqrt(ar/al)*erf((pos(1,i)-xm)/(2 * sqrt(ar*t)))
      end if
      exact = ttmp * tc + tl
      err = err + (num(i) - exact)**2
    end do
    !$omp end parallel do
    err = sqrt(err/n)
    return
  end subroutine err_infplate
end module errcalc
