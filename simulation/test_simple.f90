program test_simple
  use m_quadtree, only: QuadTree, QuadTreeNode, initializeQuadtree, insertElement, deleteTree, removeUnnecessaryNodes, &
    NodeStack, initializeStack, push, pop, cellX, cellY, cellWidth, cellHeight, &
    rightNeighborID, cellLevel
  use m_quadtree_io, only: createTree, getLeafCells, buildTreeFromMatrix, findLeaves, createTreeFromFile
  use m_simulation, only: moveParticles, selectCollisionPairs, collide, collideWithStructures
  use m_matrix_io, only: readIntMatrixFromFile, writeRealMatrixToH5, writeRealMatrixToH5Dataset
  use m_util, only: SimulationParameters
  use m_data_structures, only: addParticleCount
  use m_types, only: dp, i4, pp, i8, i2
  implicit none


  type(QuadTree), target :: tree
  type(QuadTree), pointer :: pTree
  type(QuadTreeNode), pointer :: n
  type(NodeStack) :: stack

  character(len=100) :: filename
  character(len=5) :: itCounter
  integer :: status

  integer :: i,j, it, writeStateFrequency, numTimeSteps

  real(dp) :: dt, dimX, dimY
  real(dp), allocatable :: cellMatrix(:,:)


  type(SimulationParameters), pointer :: params

  real(dp), dimension(5) :: particle
  integer(i2), dimension(:), pointer :: collisionPairs

  character(8) :: fmt

  ! fmt = "(B64.64)"
  !
  ! write(*,fmt) ibset(0_i8, 63)
  ! print *, ibset(0_i8, 63)
  ! write(*,fmt) ibset(0_i8, 62)
  ! print *, ibset(0_i8, 62)


  real(dp), dimension(4,5) :: particles

  ! dt=.25e-4_dp
  dt = .1e-1_dp * 10
  dimX = 1._dp
  dimY = 1._dp
  filename = "data/empty1x1"
  numTimeSteps = 20
  writeStateFrequency = 1

  allocate(params)

  params%V_c = dimX * dimY * 1.
  allocate(params%d_ref(1))
  allocate(params%m(1))
  params%d_ref(1) = 5e-10
  params%m(1) = 1e-25
  params%dt = dt
  params%F_N = 1e17
  params%collisionModel = 1
  params%cellHistoryLength = 3
  params%maxElementsPerCell = 1
  params%width = dimX
  params%height = dimY


  ! call initializeQuadTree(tree, 2, params)
  pTree => tree
  call createTreeFromFile(pTree, params, trim(filename)//".tree")

  allocate(tree%root%structures(1,6))
  ! tree%root%structures(:,:) = reshape([0._dp,0._dp,1._dp,1._dp, 0._dp, 0._dp], [1,6])
  tree%root%structures(:,:) = reshape([.5_dp,0._dp,.5_dp,1._dp, -1._dp, 0._dp], [1,6])
  call insertElement(tree%root, [0._dp,.2_dp, 1._dp, 0._dp, 0._dp], pTree)
  print *, tree%root%numberOfElements
  call insertElement(tree%root, [0._dp,.7_dp, 1._dp, 0._dp, 0._dp], pTree)

  print *, tree%root%numberOfElements

  do i = 1,4
    print *, tree%root%children(i)%numberOfElements
    if (associated(tree%root%children(i)%structures)) then
      print *, tree%root%children(i)%structures
    end if
  end do

  !
  ! pTree => tree
  ! call createTreeFromFile(pTree, params, trim(filename)//".tree")
  ! ! call initializeQuadTree(tree, 2, 1, dimX, dimY, params)
  !
  ! ! print *, tree%leafs%first%data%cellID
  ! ! print *, tree%leafs%last%data%cellID
  !
  ! particles = transpose(reshape(&
  !   [0._dp,0.6_dp,0.1_dp,0._dp, 0._dp&
  !   ,0._dp,0.7_dp,0.1_dp,0._dp, 0._dp&
  !   ,0.3_dp,0.1_dp,0.0_dp,0.1_dp, 0._dp&
  !   ,0.3_dp,0.9_dp,0.0_dp,-0.1_dp, 0.0_dp&
  !   ], [5,size(particles,1)]))
  ! ! call insertElement(tree%root, [0._dp,0.1_dp,0.1_dp,0._dp, 0._dp], pTree)
  ! do i = 1, size(particles,1)
  !   particle = particles(i, :)
  !   particle(1:4:2) = particle(1:4:2) * params%width
  !   particle(2:4:2) = particle(2:4:2) * params%height
  !   call insertElement(tree%root, particle, pTree)
  ! end do
  ! ! call insertElement(tree%root, [0._dp,0.6_dp,0.1_dp,0._dp, 0._dp], pTree)
  ! ! call insertElement(tree%root, [0._dp,0.73_dp,0.1_dp,0._dp, 0._dp], pTree)
  ! ! call insertElement(tree%root, [0._dp,0.3_dp,0.1_dp,0._dp, 0._dp], pTree)
  ! ! call insertElement(tree%root, [0.3_dp,0.1_dp,0.0_dp,0.1_dp, 0.0_dp], pTree)
  ! ! call insertElement(tree%root, [0.3_dp,0.9_dp,0.0_dp,-0.1_dp, 0.0_dp], pTree)
  ! ! call insertElement(tree%root, [0.3_dp,0.3_dp,0.1_dp,0._dp, 0._dp], pTree)
  !
  ! ! print *, cellLevel(tree%root)
  ! ! print *, cellLevel(tree%root%children(1))
  ! !
  ! ! print *, cellX(tree%root) * tree%width
  ! ! print *, cellY(tree%root) * tree%height
  ! !
  ! ! print *, cellX(tree%root%children(3))
  ! ! print *, cellX(tree%root%children(4))
  ! ! print *, cellY(tree%root%children(4))
  ! !
  ! ! print *, cellX(tree%root%children(1)%children(2)) * tree%width
  ! !
  ! !
  ! ! do i = 1,4
  ! !   print *, cellX(tree%root%children(i)) * tree%width
  ! !   print *, cellY(tree%root%children(i)) * tree%height
  ! ! end do
  ! !
  ! ! print *, tree%root%children(1)%numberOfElements
  !
  !> delete the last file
  open(unit=1234, iostat=status, file=trim(filename)//".h5", status='old')
  if (status == 0) close(1234, status='delete')
  
  
  call getLeafCells(pTree, cellMatrix, j)
  write(itCounter, "(i5.5)") 0
  ! call writeRealMatrixToH5Dataset(trim(filename)//"_"//itCounter//".h5", cellMatrix, j, 5+5*tree%maxElementsPerCell, status)
  call writeRealMatrixToH5Dataset(trim(filename)//".h5", itCounter, cellMatrix, j, 5+5*tree%maxElementsPerCell, status)
  deallocate(cellMatrix)
  
  do it = 1, numTimeSteps
    ! do i = 1, 4
    !   print *, tree%root%children(i)%elements
    ! end do
  
    ! print *, it
    call moveParticles(params%dt, pTree)
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
      ! call collideWithStructures(n, params)
    end do
   
  
    if (mod(it, writeStateFrequency) == 0) then
      call getLeafCells(pTree, cellMatrix, j)
      write(itCounter, "(i5.5)") it
      ! call writeRealMatrixToH5(trim(filename)//"_"//itCounter//".h5", cellMatrix, j, 5+5*tree%maxElementsPerCell, status)
      call writeRealMatrixToH5Dataset(trim(filename)//".h5", itCounter, cellMatrix, j, 5+5*tree%maxElementsPerCell, status)
      deallocate(cellMatrix)
      
      ! print *, (tree%root%children(i)%numberOfElements, i=1, 4)
      ! print *, (tree%root%children(i)%elements(:4), i=1, 4)

    end if
  end do
  ! !
  ! ! write(*,fmt) tree%root%children(1)%cellID, rightNeighborID(tree%root%children(1))
  ! ! write(*,fmt) tree%root%children(2)%cellID, rightNeighborID(tree%root%children(2))
  ! ! write(*,fmt) tree%root%children(3)%cellID, rightNeighborID(tree%root%children(3))
  ! ! write(*,fmt) tree%root%children(1)%children(1)%cellID, rightNeighborID(tree%root%children(1)%children(1))
  ! ! write(*,fmt) tree%root%children(1)%children(2)%cellID, rightNeighborID(tree%root%children(1)%children(2))
  ! !
  ! !
  ! ! print *, tree%root%children(2)%stats%particleCounter%average
  ! ! print *, tree%root%children(2)%stats%particleCounter%history
  ! ! print *, tree%root%children(2)%stats%particleCounter%currentPosition
  ! ! call addParticleCount(tree%root%children(2)%stats%particleCounter, 10)
  ! ! print *, tree%root%children(2)%stats%particleCounter%average
  ! ! print *, tree%root%children(2)%stats%particleCounter%history
  ! ! print *, tree%root%children(2)%stats%particleCounter%currentPosition
  ! ! call addParticleCount(tree%root%children(2)%stats%particleCounter, 20)
  ! ! print *, tree%root%children(2)%stats%particleCounter%average
  ! ! print *, tree%root%children(2)%stats%particleCounter%history
  ! ! print *, tree%root%children(2)%stats%particleCounter%currentPosition
  ! ! call addParticleCount(tree%root%children(2)%stats%particleCounter, 30)
  ! ! print *, tree%root%children(2)%stats%particleCounter%average
  ! ! print *, tree%root%children(2)%stats%particleCounter%history
  ! ! print *, tree%root%children(2)%stats%particleCounter%currentPosition
  ! ! call addParticleCount(tree%root%children(2)%stats%particleCounter, 40)
  ! ! print *, tree%root%children(2)%stats%particleCounter%average
  ! ! print *, tree%root%children(2)%stats%particleCounter%history
  ! ! print *, tree%root%children(2)%stats%particleCounter%currentPosition
  !
  ! call deleteTree(tree)
  ! deallocate(params)
  !
  ! print *, "hi"
  !


end program test_simple
