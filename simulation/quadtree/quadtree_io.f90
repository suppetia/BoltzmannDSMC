module m_quadtree_io
  use m_quadtree, only: QuadTreeNode, QuadTree, splitNode, NodeStack, push, pop, deleteStack, initializeStack, &
    cellWidth, cellHeight, cellX, cellY
  use m_types, only: pp, dp, i4
  implicit none

contains

  !> create a QuadTree from a matrix
  !> each node corresponds to a cell
  !> a cell is split if a value different from 0 is contained and maxDepth was not reached
  recursive subroutine createTree(tree, root, matrix, maxDepth, maxElementsPerCell, height, width, y, x)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(QuadTreeNode), pointer, intent(inout) :: root
    integer, intent(in) :: maxDepth, maxElementsPerCell
    integer(pp), allocatable, intent(in) :: matrix(:,:)
    real(dp), intent(in) :: width, height, x, y
    type(QuadTreeNode), pointer :: childNode

    integer :: rowOffset, colOffset
    ! integer :: i, j
    real(dp) :: newWidth, newHeight

    if (.not.ALLOCATED(matrix)) then
      print *, "invalid matrix"
      return
    end if

    if (maxDepth == 1) then
      root%isCollapsable = .false.
      return
    end if


    !> +1 since arrays start with 1
    rowOffset = ceiling(y)+1
    colOffset = ceiling(x)+1
    newWidth = floor(width)
    newHeight = floor(height)

    ! print *, height, floor(height)
    ! print *, y, ceiling(y)

    ! print *, all(matrix(rowOffset:rowOffset + floor(height) - 1, colOffset: colOffset + floor(width) - 1) > 0)
    ! print *, any(matrix(rowOffset:rowOffset + floor(height) - 1, colOffset: colOffset + floor(width) - 1) > 0)

    !> if all matrix elements in this cell are > 0, the cell won't be split
    if (all(matrix(rowOffset:rowOffset + floor(height) - 1, colOffset: colOffset + floor(width) - 1) > 0)) then
      return

    !> if a matrix element in this cell is /= 0 split the cell
    else if (any(matrix(rowOffset:rowOffset + floor(height) - 1, colOffset: colOffset + floor(width) - 1) > 0)) then
      call splitNode(root, tree)

      newWidth = width / 2
      newHeight = height / 2
      childNode => root%children(1)
      call createTree(tree, childNode, matrix, maxDepth-1, maxElementsPerCell, newHeight, newWidth, y, x)
      childNode => root%children(2)
      call createTree(tree, childNode, matrix, maxDepth-1, maxElementsPerCell, newHeight, newWidth, y, x+newWidth)
      childNode => root%children(3)
      call createTree(tree, childNode, matrix, maxDepth-1, maxElementsPerCell, newHeight, newWidth, y+newHeight, x)
      childNode => root%children(4)
      call createTree(tree, childNode, matrix, maxDepth-1, maxElementsPerCell, newHeight, newWidth, y+newHeight, x+newWidth)
    else
      root%isCollapsable = .false.
    end if

  end subroutine createTree

  subroutine buildTreeFromMatrix(tree, matrix)
    implicit none

    type(QuadTree), pointer, intent(inout) :: tree
    integer(pp), allocatable, dimension(:,:), intent(in) :: matrix

    call createTree(tree, tree%root, matrix, tree%maxDepth, tree%maxElementsPerCell,&
      real(size(matrix, 1), dp), real(size(matrix, 2), dp), 0._dp, 0._dp)

  end subroutine buildTreeFromMatrix


  recursive subroutine findLeaves(root, stack)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: root
    type(NodeStack), intent(inout) :: stack

    integer :: i

    if (root%numberOfElements < 0) then
      do i=1,4
        call findLeaves(root%children(i), stack)
      end do
    else
      ! if (ASSOCIATED(stack%top)) then
      ! print *, stack%top%data%numberOfElements
      ! end if
      call push(stack, root)
    end if
  end subroutine findLeaves

  subroutine getLeafCells(tree, arrLeaves, numLeaves)
    implicit none
    type(QuadTree), pointer, intent(in) :: tree
    real(dp), dimension(:,:), allocatable, intent(out) :: arrLeaves
    integer(i4), intent(out) :: numLeaves

    type(QuadTreeNode), pointer :: n
    type(NodeStack) :: stack

    integer(i4) :: i

    call initializeStack(stack)

    !> find all leaf nodes and store them in the stack
    call findLeaves(tree%root, stack)

    numLeaves = stack%topIndex

    !> retrieve the information how many elements per cell are stored and allocate the cellMatrix
    allocate(arrLeaves(numLeaves, 5+5*tree%maxElementsPerCell))

    do i = 1, numLeaves
      call pop(stack, n)
      ! arrLeaves(i, 1:5) = [n%x, n%y, n%width, n%height, real(n%numberOfElements, dp)]
      arrLeaves(i, 1:5) = [cellX(n) * tree%width, cellY(n)*tree%height, cellWidth(n)*tree%width, &
        cellHeight(n)*tree%height, real(n%numberOfElements, dp)]
      arrLeaves(i, 6:) = n%elements(:)
    end do

    call deleteStack(stack)

  end subroutine getLeafCells


  !
  ! !> plotting does not work accordingly with fractional height/width
  ! subroutine printQuadTree(tree, matrix)
  !   type(QuadTree), intent(in) :: tree
  !   integer(pp), allocatable, intent(in) :: matrix(:,:)
  !
  !   ! integer :: x,y, width, height, layer
  !   logical, allocatable :: gridMatrix(:, :, :)
  !   integer :: i,j
  !
  !   type(QuadTreeNode), pointer :: n
  !   type(NodeStack) :: childrenStack
  !
  !   call initializeStack(childrenStack)
  !
  !   call findLeaves(tree%root, childrenStack)
  !
  !   allocate(gridMatrix(size(matrix, 1), size(matrix, 2), 2))
  !
  !   gridMatrix = .false. !> initialize all values false
  !
  !
  !   do while (childrenStack%topIndex > 0)
  !     call pop(childrenStack, n)
  !
  !     ! print *, n%x, n%y, n%width, n%height
  !
  !     !> set upper border
  !     gridMatrix(ceiling(n%y), ceiling(n%x) : ceiling(n%x)+floor(n%width)-1, 1) = .true.
  !     !> set left border
  !     gridMatrix(ceiling(n%y) : ceiling(n%y)+floor(n%height)-1, ceiling(n%x), 2) = .true.
  !   end do
  !
  !   !> display the quadtree structure around the data matrix
  !   do i=1, size(matrix, 1)
  !     write(*, '(A1)', advance="no") '|'
  !     do j=1, size(matrix, 2)
  !       if (gridMatrix(i, j, 1)) then
  !         write(*, '(A3)', advance="no") '---'
  !       else
  !         write(*, '(A3)', advance="no") '   '
  !       end if 
  !       if (j == size(matrix, 2)) then
  !         exit
  !       end if
  !       write(*, '(A1)', advance="no") '-'
  !     end do
  !     write(*, '(A1)') '|'
  !     do j=1, size(matrix, 2)
  !       if (gridMatrix(i, j, 2)) then
  !         write(*, '(A1)', advance="no") '|'
  !       else
  !         write(*, '(A1)', advance="no") ' '
  !       end if
  !       write(*, '(1x,I1,1x)', advance="no") matrix(i, j)
  !     end do
  !     write(*, '(A1)') '|'
  !   end do
  !   write(*, '(A1)', advance="no") '|'
  !   do j=1, size(matrix, 2)-1
  !     write(*, '(A4)', advance="no") '----'
  !   end do
  !   write(*, '(A4)') '---|'
  !
  !   deallocate(gridMatrix)
  !   call deleteStack(childrenStack)
  !
  ! end subroutine printQuadTree
  !


end module m_quadtree_io

