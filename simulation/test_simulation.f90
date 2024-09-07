program test_simulation
  use m_quadtree, only: QuadTree, QuadTreeNode, initializeQuadtree, insertElement, deleteTree, removeUnnecessaryNodes, &
    NodeStack, initializeStack, push, pop
  use m_quadtree_io, only: createTree, getLeafCells, buildTreeFromMatrix, findLeaves
  use m_data_structures, only: addParticleCount
  use m_simulation, only: moveParticles, selectCollisionPairs, collide
  use m_matrix_io, only: readIntMatrixFromFile, writeRealMatrixToH5
  use m_types, only: dp, i4, pp, i2
  use m_util, only: SimulationParameters
  implicit none

  type(QuadTree), target :: tree
  type(QuadTree), pointer :: pTree

  integer :: i, j, it, writeStateFrequency, numTimeSteps

  integer :: nrows, ncols
  integer :: status

  integer(pp), allocatable :: matrix(:, :)
  real(dp), allocatable :: cellMatrix(:,:)

  real(dp), dimension(5) :: particle
  integer(i2), dimension(:), pointer :: collisionPairs

  character(len=100) :: filename
  character(len=5) :: itCounter

  type(NodeStack) :: stack
  type(QuadTreeNode), pointer :: n
  
  real(dp) :: dt, dimX, dimY

  type(SimulationParameters) :: params
  

  filename = "data/matrix8"
  writeStateFrequency = 10
  numTimeSteps = 25000
  dt=.25e-4_dp
  dimX = 1._dp
  dimY = 1._dp

  params%V_c = dimX * dimY * 1.
  allocate(params%d_ref(1))
  allocate(params%m(1))
  params%d_ref(1) = 5e-10
  params%m(1) = 1e-25
  params%dt = dt
  params%F_N = 1e17
  params%collisionModel = 1
  params%cellHistoryLength = 100

  
  


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


  ! call readIntMatrixFromFile(trim(filename)//".txt", matrix, nrows, ncols, status)
  ! call initializeQuadTree(tree, 10, 30, real(size(matrix, 2), dp), real(size(matrix, 1), dp))

  call initializeQuadTree(tree, 10, 30, dimX, dimY, params)
  pTree => tree
  ! call buildTreeFromMatrix(pTree, matrix)

  print *, pTree%width
  
  do i=1,1000
    call random_number(particle)
    particle(1) = particle(1) * dimX
    particle(2) = particle(2) * dimY
    particle(3:) = (particle(3:) - .5)*10
    ! print *, particle
    call insertElement(tree%root, particle, pTree)
  end do

  print *, pTree%width
  
  call getLeafCells(pTree, cellMatrix, j)
  write(itCounter, "(i5.5)") 0
  call writeRealMatrixToH5(trim(filename)//"_"//itCounter//".h5", cellMatrix, j, 5+5*tree%maxElementsPerCell, status)
  deallocate(cellMatrix)
  
  do it=1,numTimeSteps
    call moveParticles(dt, pTree)
    call removeUnnecessaryNodes(tree%root)

    call initializeStack(stack)
    call findLeaves(tree%root, stack)
    do i = 1,stack%topIndex
      call pop(stack, n)
      call addParticleCount(n%stats%particleCounter, n%numberOfElements)!> move this somewhere else
      ! print *, n%numberOfElements
      call selectCollisionPairs(n, params, collisionPairs, j)
      ! print *, collisionPairs(:)
      call collide(n, collisionPairs, j, params)
    end do
    
  
    if (mod(it, writeStateFrequency) == 0) then
      call getLeafCells(pTree, cellMatrix, j)
      write(itCounter, "(i5.5)") it
      call writeRealMatrixToH5(trim(filename)//"_"//itCounter//".h5", cellMatrix, j, 5+5*tree%maxElementsPerCell, status)
      deallocate(cellMatrix)
    end if
  end do 


  call deleteTree(tree)

  ! deallocate(matrix)


  

end program test_simulation
