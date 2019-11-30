module n2ext
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2Cv, n2R, n2Name, setdimbase

  private

    real :: n2C(3) = (/ 1./120., 7./(478. * pi), 1./(120. * pi) /)
    character (len=10) :: n2Name = ' genesis '
    real :: n2R = 3.0, n2Cv
    integer :: dim
  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    n2Cv = n2C(dim)
  end subroutine

  pure subroutine n2f(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f

    real              :: q
    q = r / h
    if (q >= 3.) then
      f = 0.
    else if (q >= 2.) then
      f = (3. - q)**5
    else if (q >= 1.) then
      f = (3. - q)**5 - 6. * (2. - q)**5
    else if (q >= 0.) then
      f = (3. - q)**5 - 6. * (2. - q)**5 + 15. * (1. - q)**5
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
    end subroutine

  pure subroutine n2df(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df

    real              :: q
    q = r / h
    if (q >= 3.) then
      df = 0.
    else if (q >= 2.) then
      df = -5. * (3. - q)**4 / q
    else if (q >= 1.) then
      df = (-5. * (3. - q)**4 + 30. * (2. - q)**4) / q
    else if (q > 0.) then
      df = (-5. * (3. - q)**4 + 30. * (2. - q)**4 - 75. * (1. - q)**4) / q
    else if (q == 0.) then
      df = 0.
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine

  pure subroutine n2ddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf

    real              :: q
    q = r / h
    if (q >= 3.) then
      ddf = 0.
    else if (q >= 2.) then
      ddf = 20. * (3. - q)**3
    else if (q >= 1.) then
      ddf = 20. * (3. - q)**3 - 120. * (2. - q)**3
    else if (q >= 0.) then
      ddf = 20. * (3. - q)**3 - 120. * (2. - q)**3 + 300. * (1. - q)**3
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine
end module n2ext
