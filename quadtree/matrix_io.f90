module m_matrix_io
  use m_types, only: fp, pp, dp
  implicit none

contains
  
  subroutine readIntMatrixFromFile(filename, matrix, nrows, ncols, status)
    character(len=*), intent(in) :: filename
    integer(pp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: nrows, ncols
    integer, intent(out) :: status

    integer :: i, j
    integer :: io_status

    integer :: io = 10

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
  end subroutine readIntMatrixFromFile

  subroutine writeRealMatrixToFile(filename, matrix, nrows, ncols, error)
    implicit none
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(in) :: matrix(:,:)
    integer, intent(in) :: nrows, ncols
    integer, intent(out) :: error

    integer :: i, j
    integer :: io_status

    integer :: io = 10

    !> Open the file for reading
    open(unit=io, file=filename, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) then
      print *, 'Error opening file: ', trim(filename)
      error = io_status
      return
    end if

    !> Read the matrix elements
    do i = 1, nrows
      do j = 1, ncols
        write(io, '(F8.2)', advance='no') matrix(i, j)
        if (j < ncols) then
          write(io, '(A)', advance='no') ' '
        else
          write(io, '(A)', advance='yes') ''
        end if
      end do
    end do

    !> Close the file
    close(io)

    !> Set the status to indicate success
    error = 0

  end subroutine writeRealMatrixToFile

  subroutine writeRealMatrixToH5(filename, matrix, nrows, ncols, error)
    use hdf5
    implicit none
    character(len=*), intent(in) :: filename

    real(dp), dimension(:,:), intent(in) :: matrix
    integer, intent(in) :: nrows, ncols
    integer, intent(out) :: error

    character(len=*), parameter :: datasetName = "gridCells"
    integer(hid_t) :: fileID !> file identifier
    integer(hid_t) :: dspaceID !> dataspace identifier
    integer(hid_t) :: dsetID !> dataset identifier

    logical :: datasetExists

    integer(hsize_t), dimension(2) :: dataDims

    dataDims(:) = [nrows, ncols]

    !> initialize Fortran interface
    call h5open_f(error)
    !> open the file and overwrite existing files
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileID, error)


    !> create the dataspace
    call h5screate_simple_f(2, dataDims, dspaceID, error)

    !> create the dataset with default properties
    call h5dcreate_f(fileID, datasetName, H5T_NATIVE_DOUBLE, dspaceID, dsetID, error)

    !> write the data from the matrix to the dataset
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, matrix, dataDims, error)


    !> close all accesses
    call h5dclose_f(dsetID, error)
    call h5sclose_f(dspaceID, error)
    call h5fclose_f(fileID, error)
    call h5close_f(error)
  end subroutine writeRealMatrixToH5

end module m_matrix_io
