module m_matrix_io
  use m_types, only: fp, pp
  implicit none

contains
  
  subroutine readMatrixFromFile(filename, matrix, nrows, ncols, status)
    character(len=*), intent(in) :: filename
    integer(pp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: nrows, ncols
    integer, intent(out) :: status

    integer :: i, j
    integer :: io_status

    integer :: io

    !> Open the file for reading
    open(unit=io, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, 'Error opening file: ', trim(filename)
      status = io_status
      return
    end if

    !> Read the dimensions of the matrix
    read(io, *) nrows, ncols

    !> Allocate the matrix array
    allocate(matrix(nrows, ncols))

    !> Read the matrix elements
    do i = 1, nrows
      read(io, *) (matrix(i, j), j = 1, ncols)
    end do

    !> Close the file
    close(io)

    !> Set the status to indicate success
    status = 0
  end subroutine readMatrixFromFile

end module m_matrix_io
