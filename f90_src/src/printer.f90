module printer

  use timing, only: addTime

  implicit none

  public :: Output, AppendLine

  private

  integer, save :: ifile = 0
  integer(8) :: start=0, finish=0


contains
  subroutine Output(time, ptype, x, v, dv, m, den, slen, pres, ien, cf, err)
    real, allocatable, intent(in)    :: x(:,:), v(:,:), dv(:,:), m(:), err(:), &
                                        den(:), slen(:), pres(:), ien(:), cf(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, intent(in)    :: time
    character (len=40)  :: fname
    integer :: iu, j, n

    call system_clock(start)

    n = size(ptype)
    write(fname, "(a,i5.5)") 'output/step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
      ! if (ptype(j) == 1) then
      write(iu, *) x(:,j), v(:,j), dv(:,j), m(j), den(j), slen(j), pres(j), ien(j), cf(j), err(j)
      ! end if
    end do
    close(iu)
    ifile = ifile + 1
    call system_clock(finish)
    call addTime(' printer', finish - start)
  end subroutine Output

  subroutine AppendLine(A, fname, t)
    real, allocatable, intent(in) :: A(:)
    character (len=*), intent(in) :: fname
    integer(8), intent(out) :: t
    integer :: iu
    logical :: exist

    call system_clock(start)

    inquire(file=fname, exist=exist)
    if (exist) then
      open(newunit=iu, file=fname, status='old', form='formatted', access='append')
    else
      open(newunit=iu, file=fname, status='new', form='formatted')
    end if
    write(iu, *) A(:)
    close(iu)

    call system_clock(finish)
    t = finish - start
    call addTime(' printer', t)
  end subroutine AppendLine
end module printer
