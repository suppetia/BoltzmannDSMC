module m_matrix_io  !
  use m_types, only: fp, i4
  implicit none
contains

  subroutine writeRealMatrixToH5Dataset(filename, datasetname, matrix, nrows, ncols, error)
    use hdf5
    implicit none
    character(len=*), intent(in) :: filename, datasetname

    real(fp), dimension(:,:), intent(in) :: matrix
    integer(i4), intent(in) :: nrows, ncols
    integer(i4), intent(out) :: error

    integer(hid_t) :: fileID !> file identifier
    integer(hid_t) :: dspaceID !> dataspace identifier
    integer(hid_t) :: dsetID !> dataset identifier

    logical :: fileExists

    integer(hsize_t), dimension(2) :: dataDims

    dataDims(:) = [ncols, nrows]

    !> initialize Fortran interface
    call h5open_f(error)

    !> check if a file exists
    inquire(file=filename, exist=fileExists)
    if (.not.fileExists) then
      !> create the file if it does exist
      !> the file is now open
      call h5fcreate_f(filename, H5F_ACC_EXCL_F, fileID, error)
    else
      !> otherwise open the existing file
      call h5fopen_f(filename, H5F_ACC_RDWR_F, fileID, error)
    end if

    !> create the dataspace
    call h5screate_simple_f(2, dataDims, dspaceID, error)

    !> create the dataset with default properties
    call h5dcreate_f(fileID, datasetname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, error)

    !> write the data from the matrix to the dataset
    call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, matrix, dataDims, error)


    !> close all accesses
    call h5dclose_f(dsetID, error)
    call h5sclose_f(dspaceID, error)
    call h5fclose_f(fileID, error)
    call h5close_f(error)
  end subroutine writeRealMatrixToH5Dataset
end module m_matrix_io 
