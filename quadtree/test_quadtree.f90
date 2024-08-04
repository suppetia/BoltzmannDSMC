program testQuadTree
  use m_types, only: fp, pp, dp
  use m_quadtree, only: QuadTree, QuadTreeNode, splitNode, initializeQuadTree, deleteSubTree, &
    insertElement, removeElementFromNode, removeUnnecessaryNodes
  use m_quadtree_io, only: createTree, printQuadTree, getLeafCells, buildTreeFromMatrix,&
    NodeStack, initializeStack, push, pop, findLeaves
  use m_matrix_io, only: readIntMatrixFromFile, writeRealMatrixToFile, writeRealMatrixToH5
  implicit none

  type(QuadTree) :: tree

  integer :: i, j

  integer :: nrows, ncols
  integer :: status

  integer(pp), allocatable :: matrix(:, :)
  real(dp), allocatable :: cellMatrix(:,:)

  real(dp), dimension(4) :: particle

  character(len=100) :: filename, filenameCells

  type(NodeStack) :: stack
  type(QuadTreeNode), pointer :: n
  

  filename = "matrix8.txt"
  ! filenameCells = "matrix8_cells.txt"
  filenameCells = "matrix8_cells.h5"

  call readIntMatrixFromFile(filename, matrix, nrows, ncols, status)

  ! if (status == 0) then
  !   print *, 'Matrix read from file:'
  !   do i = 1, nrows
  !     print *, (matrix(i, j), j = 1, ncols)
  !   end do
  ! else
  !   print *, 'Failed to read the matrix from file.'
  ! end if

  !> basic QuadTree tests

  ! print *, allocated(tree%root%children)
  !
  ! call splitNode(tree%root)
  !
  ! print *, allocated(tree%root%children)
  !
  ! call splitNode(tree%root%children(1))
  !
  ! do i = 1, 4
  !   print *, tree%root%children(i)%elements(:)
  !   print *, tree%root%children(i)%numberOfElements
  ! end do


  ! call deleteSubTree(tree%root)

  ! print *, allocated(tree%root%children)


  ! print *, any(matrix(:, :) > 0)


  ! tree%root%width = size(matrix,2)
  ! tree%root%height = size(matrix,1)

  ! print *, "root x", tree%root%x
  ! print *, "root y", tree%root%y

  ! call createTree(tree%root, matrix, 15, real(nrows, dp), real(ncols, dp), 1._dp, 1._dp)
  ! call printQuadTree(tree, matrix)


  call initializeQuadTree(tree, 15, 20, real(size(matrix, 2), dp), real(size(matrix, 1), dp))
  print *, tree%root%numberOfElements
  call buildTreeFromMatrix(tree, matrix)
  print *, tree%root%numberOfElements
  
  ! allocate(cellMatrix(100000, 5 + 4*tree%maxElementsPerCell))
  call getLeafCells(tree, cellMatrix, j)
  
  print *, j
  deallocate(cellMatrix)
  
  do i=1,1000000
    call random_number(particle)
    particle(1) = particle(1) * real(size(matrix, 2), dp)
    particle(2) = particle(2) * real(size(matrix, 1), dp)
    ! print *, particle
    call insertElement(tree%root, particle)
  end do
  
  call getLeafCells(tree, cellMatrix, j)
  
  print *, j
  deallocate(cellMatrix)
  ! do i = 1, j
  !   print *, cellMatrix(i, :)
  ! end do 
  
  !> remove 100 points
  call initializeStack(stack, 100000)
  call findLeaves(tree%root, stack)
  
  j = stack%topIndex
  i = 0
  do while (i < 300)
    call pop(stack, n)
    ! print *, loc(n)
    if (n%numberOfElements > 0) then
      do j = 1,n%numberOfElements
        call removeElementFromNode(n, 1)
        i = i + 1
      end do
    end if
  end do
  
  n => tree%root
  call removeUnnecessaryNodes(n)
  
  call getLeafCells(tree, cellMatrix, j)
  print *, j
  
  ! call writeRealMatrixToFile(filenameCells, cellMatrix, j, 5+ 4*tree%maxElementsPerCell, status)
  call writeRealMatrixToH5(filenameCells, cellMatrix, j, 5+ 4*tree%maxElementsPerCell, status)
  
  call deleteSubTree(tree%root)
  
  deallocate(cellMatrix)
  deallocate(matrix)
  

end program testQuadTree

