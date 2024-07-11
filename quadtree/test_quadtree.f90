program testQuadTree
  use m_types, only: fp, pp, dp
  use m_quadtree, only: QuadTree, QuadTreeNode, splitNode, deleteSubTree
  use m_quadtree_io, only: createTree, printQuadTree
  use m_matrix_io, only: readMatrixFromFile
  implicit none

  type(QuadTree) :: tree = QuadTree(5)

  integer :: i, j

  integer :: nrows, ncols
  integer :: status

  integer(pp), allocatable :: matrix(:, :)

  character(len=100) :: filename

  filename = "matrix3.txt"

  call readMatrixFromFile(filename, matrix, nrows, ncols, status)

  ! if (status == 0) then
  !   print *, 'Matrix read from file:'
  !   do i = 1, nrows
  !     print *, (matrix(i, j), j = 1, ncols)
  !   end do
  ! else
  !   print *, 'Failed to read the matrix from file.'
  ! end if

  !> basic QuadTree tests

  print *, allocated(tree%root%children)

  call splitNode(tree%root)

  print *, allocated(tree%root%children)

  call splitNode(tree%root%children(1))

  do i = 1, 4
    print *, tree%root%children(i)%value
  end do


  call deleteSubTree(tree%root)

  ! print *, allocated(tree%root%children)


  ! print *, any(matrix(:, :) > 0)

  tree%root%width = size(matrix,2)
  tree%root%height = size(matrix,1)

  ! print *, "root x", tree%root%x
  ! print *, "root y", tree%root%y

  call createTree(tree%root, matrix, 6, real(nrows, dp), real(ncols, dp), 1._dp, 1._dp)
  ! call printQuadTree(tree, matrix)


end program testQuadTree

