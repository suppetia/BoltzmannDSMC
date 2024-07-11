module m_quadtree_io
  use m_quadtree, only: QuadTreeNode, QuadTree, splitNode
  use m_types, only: pp
  implicit none

  type :: NodeStack
    integer :: stackSize

    integer :: topIndex = 0
    type(QuadTreeNode), allocatable :: elements(:)

  end type NodeStack

contains
  
  subroutine initializeStack(stack, stackSize)
    type(NodeStack), intent(inout) :: stack
    integer, intent(in) :: stackSize

    allocate(stack%elements(stackSize))
    stack%stackSize = stackSize
    stack%topIndex = 0

  end subroutine initializeStack

  subroutine deleteStack(stack)
    type(NodeStack), intent(inout) :: stack

    deallocate(stack%elements)
    stack%stackSize = 0
    stack%topIndex = 0
  end subroutine deleteStack

  subroutine push(stack, newValue)
    type(NodeStack), intent(inout) :: stack
    type(QuadTreeNode), intent(in) :: newValue

    if (stack%topIndex >= stack%stackSize) then
      print *, "Stack full"
      return
    end if

    stack%topIndex = stack%topIndex + 1
    stack%elements(stack%topIndex) = newValue
  end subroutine push

  subroutine pop(stack, value)
    type(NodeStack), intent(inout) :: stack
    type(QuadTreeNode), intent(out) :: value

    if (stack%topIndex <= 0) then
      print *, "Stack empty"
      return
    end if

    value = stack%elements(stack%topIndex)
    stack%topIndex = stack%topIndex - 1

  end subroutine pop


  !> create a QuadTree from a matrix
  !> each node corresponds to a cell
  !> a cell is split if a value different from 0 is contained and maxDepth was not reached
  recursive subroutine createTree(root, matrix, maxDepth, height, width, rowOffset, colOffset)
    implicit none
    type(QuadTreeNode), intent(inout) :: root
    integer, intent(in) :: maxDepth
    integer(pp), allocatable, intent(in) :: matrix(:,:)
    integer, intent(in) :: width, height, rowOffset, colOffset

    ! integer :: i, j
    integer :: newWidth, newHeight

    if (.not.ALLOCATED(matrix)) then
      print *, "invalid matrix"
      return
    end if

    if (maxDepth == 1) then
      return
    end if

    !> if all matrix elements in this cell are > 0, the cell won't be split
    if (all(matrix(rowOffset:rowOffset + height - 1, colOffset: colOffset + width - 1) > 0)) then
      return

    !> if a matrix element in this cell is /= 0 split the cell
    else if (any(matrix(rowOffset:rowOffset + height - 1, colOffset: colOffset + width - 1) > 0)) then
      call splitNode(root)
      newWidth = width / 2
      newHeight = height / 2
      call createTree(root%children(1), matrix, maxDepth-1, newHeight, newWidth, rowOffset, colOffset)
      call createTree(root%children(2), matrix, maxDepth-1, newHeight, newWidth, rowOffset, colOffset+newWidth)
      call createTree(root%children(3), matrix, maxDepth-1, newHeight, newWidth, rowOffset+newHeight, colOffset)
      call createTree(root%children(4), matrix, maxDepth-1, newHeight, newWidth, rowOffset+newHeight, colOffset+newWidth)
    end if

  end subroutine createTree


  recursive subroutine findChildren(root, stack)
    implicit none
    type(QuadTreeNode), intent(in) :: root
    type(NodeStack), intent(inout) :: stack

    integer :: i

    if (root%value < 0) then
      do i=1,4
        call findChildren(root%children(i), stack)
      end do
    else
      call push(stack, root)
    end if
  end subroutine findChildren


  subroutine printQuadTree(tree, matrix)
    type(QuadTree), intent(in) :: tree
    integer(pp), allocatable, intent(in) :: matrix(:,:)

    ! integer :: x,y, width, height, layer
    logical, allocatable :: gridMatrix(:, :, :)
    integer :: i,j

    type(QuadTreeNode) :: n
    type(NodeStack) :: childrenStack

    call initializeStack(childrenStack, 1000) !> use a reasonable heuristic here

    call findChildren(tree%root, childrenStack)

    allocate(gridMatrix(size(matrix, 1), size(matrix, 2), 2))

    gridMatrix = .false. !> initialize all values false


    do while (childrenStack%topIndex > 0)
      call pop(childrenStack, n)

      print *, n%x, n%y, n%width, n%height

      !> set upper border
      gridMatrix(n%y, n%x : n%x+n%width-1, 1) = .true.
      !> set left border
      gridMatrix(n%y : n%y+n%height-1, n%x, 2) = .true.
    end do

    !> display the quadtree structure around the data matrix
    do i=1, size(matrix, 1)
      write(*, '(A1)', advance="no") '|'
      do j=1, size(matrix, 2)
        if (gridMatrix(i, j, 1)) then
          write(*, '(A3)', advance="no") '---'
        else
          write(*, '(A3)', advance="no") '   '
        end if 
        if (j == size(matrix, 2)) then
          exit
        end if
        write(*, '(A1)', advance="no") '-'
      end do
      write(*, '(A1)') '|'
      do j=1, size(matrix, 2)
        if (gridMatrix(i, j, 2)) then
          write(*, '(A1)', advance="no") '|'
        else
          write(*, '(A1)', advance="no") ' '
        end if
        write(*, '(1x,I1,1x)', advance="no") matrix(i, j)
      end do
      write(*, '(A1)') '|'
    end do
    write(*, '(A1)', advance="no") '|'
    do j=1, size(matrix, 2)-1
      write(*, '(A4)', advance="no") '----'
    end do
    write(*, '(A4)') '---|'

    deallocate(gridMatrix)
    call deleteStack(childrenStack)

  end subroutine printQuadTree



end module m_quadtree_io

