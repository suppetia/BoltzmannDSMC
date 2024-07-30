module m_quadtree
  use m_types, only: fp, pp, i4, dp
  implicit none

  type :: QuadTreeNode
    integer(pp) :: numberOfElements = 0
    real(dp), allocatable, dimension(:) :: elements

    !> it should be possible to calculate this on the fly
    !> but I'm not capable of wrapping my head around this
    real(dp) :: x, y
    real(dp) :: width, height

    type(QuadTreeNode), allocatable, dimension(:) :: children

  end type QuadTreeNode

  type :: QuadTree
    integer :: maxDepth
    integer :: maxElementsPerCell = 1

    type(QuadTreeNode) :: root = QuadTreeNode(x=0, y=0, width=0, height=0)

  end type QuadTree

contains

  !> initialize a new QuadTreeNode
  subroutine initializeQuadTreeNode(node, maxElements, x, y, width, height)
    implicit none
    type(QuadTreeNode), intent(out) :: node
    integer, intent(in) :: maxElements
    real(dp), intent(in) :: x, y, width, height

    node%x = x
    node%y = y
    node%width = width
    node%height = height

    node%numberOfElements = 0
    allocate(node%elements(4*maxElements))
    node%elements(:) = -1

  end subroutine initializeQuadTreeNode
    

  !> initialize a new QuadTree
  subroutine initializeQuadTree(tree, maxDepth, maxElementsPerCell, width, height)
    implicit none
    integer, intent(in) :: maxDepth
    integer, intent(in) :: maxElementsPerCell
    real(dp), intent(in) :: width, height
  
    type(QuadTree), intent(out) :: tree

    tree%maxDepth = maxDepth
    tree%maxElementsPerCell = maxElementsPerCell
    call initializeQuadTreeNode(tree%root, maxElementsPerCell, 0._dp, 0._dp, width, height)
  end subroutine initializeQuadTree

  !> add children to a node
  subroutine splitNode(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node

    real(dp) :: newWidth, newHeight
    integer :: i, j
    integer, dimension(4) :: idx

    allocate(node%children(4))
    idx(:) = 0

    newWidth = node%width / 2
    newHeight = node%height / 2

    call initializeQuadTreeNode(node%children(1), size(node%elements,1)/4, node%x, node%y, newWidth, newHeight)
    call initializeQuadTreeNode(node%children(2), size(node%elements,1)/4, node%x+newWidth, node%y, newWidth, newHeight)
    call initializeQuadTreeNode(node%children(3), size(node%elements,1)/4, node%x, node%y+newHeight, newWidth, newHeight)
    call initializeQuadTreeNode(node%children(4), size(node%elements,1)/4, node%x+newWidth, node%y+newHeight, newWidth, newHeight)

    do i = 0,node%numberOfElements-1
      !> find the child where to sort the element in
      if (node%elements(4*i+1) < node%x + newWidth) then
        !> left side
        j = 1
      else
        !> right side
        j = 2
      end if
      if (node%elements(4*i+2) > node%y + newHeight) then
        !> lower row
        j = j+2
      end if
      !> copy the element entries to the corresponding child
      node%children(j)%elements(4*idx(j)+1:4*(idx(j)+1)) = node%elements(4*i+1:4*(i+1))
      idx(j) = idx(j)+1
    end do

    do i=1,4
      node%children(i)%numberOfElements = idx(i)
    end do
    
    !> mark that the node is no leave
    node%numberOfElements= -1
    deallocate(node%elements)

  end subroutine splitNode

  !> merge the children making node a leaf
  subroutine mergeChildren(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node

    integer :: i,j, idx

    !> get the number of elements from a child node
    allocate(node%elements(size(node%children(1)%elements,1)))
    node%elements(:) = -1

    idx = 0

    do i=1,4
      do j=0,node%children(i)%numberOfElements-1
        node%elements(4*idx+1:4*(idx+1)) = node%children(i)%elements(4*j+1:4*(j+1))
        deallocate(node%elements)
        idx = idx+1
      end do
    end do
    node%numberOfElements = idx
    deallocate(node%children)
  end subroutine mergeChildren


  recursive subroutine removeUnnecessaryNodes(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node

    integer :: i,numElements

    if (node%numberOfElements >= 0) then
      return
    end if

    numElements = 0
    do i=1,4
      call removeUnnecessaryNodes(node%children(i))
      if (node%children(i)%numberOfElements < 0) then
        return
      end if
      numElements = numElements + node%children(i)%numberOfElements
    end do
    if (numElements <= size(node%children(1)%elements,1)/4) then
      call mergeChildren(node)
    end if
  end subroutine removeUnnecessaryNodes

  !> delete a QuadTree by freeing all memory
  recursive subroutine deleteSubTree(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node
    integer :: i
    
    if (allocated(node%children)) then
      do i = 1,4
        call deleteSubTree(node%children(i))
      end do
      deallocate(node%children)
    end if
  end subroutine deleteSubTree

  recursive subroutine insertElement(node, element)
    implicit none
    type(QuadTreeNode), intent(inout) :: node
    real(dp), dimension(4), intent(in) :: element
    integer :: i

    if (node%numberOfElements >= 0) then
      if (node%numberOfElements == size(node%elements, 1)/4) then
        !> if the cell is full, split it
        call splitNode(node)
        call insertElement(node, element)
      else
        !> otherwise insert the new element
        i = node%numberOfElements
        node%elements(4*i+1:4*(i+1)) = element(:)
        node%numberOfElements = i+1
      end if
    else
      if (element(1) < node%x + node%width/2) then
        i = 1
      else
        i = 2
      end if
      if (element(2) >= node%y + node%height/2) then
        i = i + 2
      end if
      call insertElement(node%children(i), element)
    end if

  end subroutine insertElement

!   subroutine insertElement(tree, element)
!     implicit none
!     type(QuadTree), intent(inout) :: tree
!     real(dp), dimension(4), intent(in) :: element
!     type(QuadTreeNode) :: node
!     integer :: i
!
!     print *, loc(tree%root)
!     node = tree%root
!     print *, loc(node)
!     do while (.true.)
!       print *, node%numberOfElements
!       if (node%numberOfElements >= 0) then
!         if (node%numberOfElements == size(node%elements, 1)/4) then
!           !> if the cell is full, split it
!           call splitNode(node)
!         else
!           !> otherwise insert the new element
!           i = node%numberOfElements
!           node%elements(4*i+1:4*(i+1)) = element(:)
!           node%numberOfElements = i+1
!           exit
!         end if
!       end if
!       !> find the last child node
!       ! print *, node%numberOfElements
!       if (element(1) < node%x + node%width/2) then
!         i = 1
!       else
!         i = 2
!       end if
!       if (element(2) >= node%y + node%height/2) then
!         i = i + 2
!       end if
!       !> this is not working since node is not actually overwritten
!       node = node%children(i)
!     end do
!
!   end subroutine insertElement

  subroutine removeElementFromNode(node, i)
    implicit none
    type(QuadTreeNode), intent(inout) :: node
    integer, intent(in) :: i

    integer :: j

    j = node%numberOfElements

    if (i > j) then
      print *, "requested element ID not in node"
      return
    end if
    !> copy the last element to index i
    node%elements(4*(i-1)+1:4*i) = node%elements(4*(j-1)+1:4*j)
    !> remove the last element
    node%elements(4*(j-1)+1:4*j) = -1
    node%numberOfElements = j - 1

  end subroutine removeElementFromNode

end module m_quadtree

