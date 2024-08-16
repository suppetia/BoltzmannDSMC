program test_simulation
  use m_quadtree, only: QuadTree, QuadTreeNode, initializeQuadtree, insertElement, deleteTree, removeUnnecessaryNodes, &
    NodeStack, initializeStack, push, pop
  use m_quadtree_io, only: createTree, getLeafCells, buildTreeFromMatrix
  use m_simulation, only: moveParticles
  use m_matrix_io, only: readIntMatrixFromFile, writeRealMatrixToH5
  use m_types, only: dp, i4, pp
  implicit none

  type(QuadTree), target :: tree
  type(QuadTree), pointer :: pTree

  integer :: i, j, it, writeStateFrequency, numTimeSteps

  integer :: nrows, ncols
  integer :: status

  integer(pp), allocatable :: matrix(:, :)
  real(dp), allocatable :: cellMatrix(:,:)

  real(dp), dimension(4) :: particle

  character(len=100) :: filename
  character(len=5) :: itCounter

  type(NodeStack) :: stack
  type(QuadTreeNode), pointer :: n
  
  real(dp) :: dt
  

  filename = "data/matrix8"
  writeStateFrequency = 10
  numTimeSteps = 250
  dt=5._dp


  ! call initializeQuadTree(tree, 2, 1, 1._dp, 1._dp)
  !
  ! call insertElement(tree%root, [0._dp,0._dp,1._dp,0._dp])
  ! call insertElement(tree%root, [0._dp,0.8_dp,1._dp,0._dp])
  !
  ! print *, tree%root%children(1)%numberOfElements
  ! 
  ! do j = 1, 5
  !   do i = 1, 4
  !     print *, tree%root%children(i)%elements
  !   end do
  !  
  !   call moveParticles(tree, .3_dp)
  ! end do
  !
  ! call deleteTree(tree)


  call readIntMatrixFromFile(trim(filename)//".txt", matrix, nrows, ncols, status)
  call initializeQuadTree(tree, 10, 15, real(size(matrix, 2), dp), real(size(matrix, 1), dp))
  print *, tree%width

  pTree => tree
  print *, pTree%width
  call buildTreeFromMatrix(pTree, matrix)

  print *, pTree%width
  
  do i=1,1000
    call random_number(particle)
    particle(1) = particle(1) * real(size(matrix, 2), dp)
    particle(2) = particle(2) * real(size(matrix, 1), dp)
    ! print *, particle
    call insertElement(tree%root, particle, pTree)
  end do

  print *, pTree%width
  
  call getLeafCells(pTree, cellMatrix, j)
  write(itCounter, "(i5.5)") 0
  call writeRealMatrixToH5(trim(filename)//"_"//itCounter//".h5", cellMatrix, j, 5+4*tree%maxElementsPerCell, status)
  deallocate(cellMatrix)
  
  do it=1,numTimeSteps
    call moveParticles(dt, pTree)
    call removeUnnecessaryNodes(tree%root)
  
    if (mod(it, writeStateFrequency) == 0) then
      call getLeafCells(pTree, cellMatrix, j)
      write(itCounter, "(i5.5)") it
      call writeRealMatrixToH5(trim(filename)//"_"//itCounter//".h5", cellMatrix, j, 5+4*tree%maxElementsPerCell, status)
      deallocate(cellMatrix)
    end if
  end do 


  call deleteTree(tree)

  deallocate(matrix)


  

end program test_simulation
