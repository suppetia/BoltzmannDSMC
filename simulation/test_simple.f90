program test_simple
  use m_quadtree, only: QuadTree, QuadTreeNode, initializeQuadtree, insertElement, deleteTree, removeUnnecessaryNodes, &
    NodeStack, initializeStack, push, pop, cellX, cellY, cellWidth, cellHeight, &
    rightNeighborID, cellLevel
  use m_quadtree_io, only: createTree, getLeafCells, buildTreeFromMatrix
  use m_simulation, only: moveParticles
  use m_matrix_io, only: readIntMatrixFromFile, writeRealMatrixToH5
  use m_types, only: dp, i4, pp, i8
  implicit none


  type(QuadTree), target :: tree
  type(QuadTree), pointer :: pTree
  type(QuadTreeNode), pointer :: n
  type(NodeStack) :: stack

  integer :: i,j

  character(8) :: fmt

  fmt = "(B64.64)"

  write(*,fmt) ibset(0_i8, 63)
  print *, ibset(0_i8, 63)
  write(*,fmt) ibset(0_i8, 62)
  print *, ibset(0_i8, 62)


  pTree => tree
  call initializeQuadTree(tree, 2, 1, 1._dp, 1._dp)

  print *, tree%leafs%first%data%cellID
  print *, tree%leafs%last%data%cellID

  
  call insertElement(tree%root, [0._dp,0._dp,1._dp,0._dp], pTree)
  call insertElement(tree%root, [0._dp,0.8_dp,1._dp,0._dp], pTree)
  call insertElement(tree%root, [0._dp,0.3_dp,1._dp,0._dp], pTree)

  print *, cellLevel(tree%root)
  print *, cellLevel(tree%root%children(1))
  
  print *, cellX(tree%root) * tree%width
  print *, cellY(tree%root) * tree%height

  print *, cellX(tree%root%children(1)%children(2)) * tree%width
  

  do i = 1,4
    print *, cellX(tree%root%children(i)) * tree%width
    print *, cellY(tree%root%children(i)) * tree%height
  end do
  
  print *, tree%root%children(1)%numberOfElements
  
  ! do j = 1, 5
  !   do i = 1, 4
  !     print *, tree%root%children(i)%elements
  !   end do
  ! 
  !   call moveParticles(.3_dp, tree)
  ! end do

  write(*,fmt) tree%root%children(1)%cellID, rightNeighborID(tree%root%children(1))
  write(*,fmt) tree%root%children(2)%cellID, rightNeighborID(tree%root%children(2))
  write(*,fmt) tree%root%children(3)%cellID, rightNeighborID(tree%root%children(3))
  write(*,fmt) tree%root%children(1)%children(1)%cellID, rightNeighborID(tree%root%children(1)%children(1))
  write(*,fmt) tree%root%children(1)%children(2)%cellID, rightNeighborID(tree%root%children(1)%children(2))

  
  
  call deleteTree(tree)



end program test_simple
