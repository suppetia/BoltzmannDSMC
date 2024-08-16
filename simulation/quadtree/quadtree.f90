module m_quadtree
  use m_types, only: fp, pp, i4, dp, i8, i2
  implicit none

  type :: QuadTreeNode
    integer(pp) :: numberOfElements

    ! !> it should be possible to calculate this on the fly
    ! !> but I'm not capable of wrapping my head around this
    ! real(dp) :: x, y
    ! real(dp) :: width, height

    integer(i8) :: cellID = 0
    integer(i2) :: cellLevel = 0

    !> whether a grid cell is allowed to be merged
    !> false if it is given from the structure
    logical :: isCollapsable = .true.

    real(dp), pointer, dimension(:) :: elements => null()
    type(QuadTreeNode), pointer, dimension(:) :: children => null()

  end type QuadTreeNode

  type :: StackNode
    type(QuadTreeNode), pointer :: data
    type(StackNode), pointer :: next => null()
  end type StackNode

  type :: NodeStack
    integer :: topIndex = 0
    type(StackNode), pointer :: first => null()
    type(StackNode), pointer :: top => null()
  end type NodeStack

  type :: QuadTree
    integer :: maxDepth
    integer :: maxElementsPerCell = 1
    real(dp) :: width, height
    type(NodeStack), pointer :: leafs

    type(QuadTreeNode), pointer :: root
  end type QuadTree

contains

  function cellWidth(node) result (width)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: width
    
    ! print *, node%cellLevel
    width = 2._dp**(-node%cellLevel)
  end function cellWidth

  function cellHeight(node) result (height)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: height
    
    height = 2._dp**(-node%cellLevel)
  end function cellHeight

  function cellX(node) result (x)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: width
    integer(i2) :: level
    real(dp) :: x

    width = 1._dp
    x = 0._dp

    do level = 0,node%cellLevel-1_i2
      width = width * .5
      x = x + width * ibits(node%cellID, level*2, 1)
    end do
  end function cellX
    
  function cellY(node) result (y)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: height
    integer(i2) :: level
    real(dp) :: y

    height = 1._dp
    y = 0._dp

    do level = 0,node%cellLevel-1_i2
      height = height * .5
      y = y + height * ibits(node%cellID, level*2+1, 1)
    end do
  end function cellY
    

  subroutine initializeStack(stack)
    type(NodeStack), intent(inout) :: stack

    nullify(stack%first)
    nullify(stack%top)
    stack%topIndex = 0
  end subroutine initializeStack

  subroutine deleteStack(stack)
    type(NodeStack), intent(inout) :: stack
    integer :: i
    type(StackNode), pointer :: n1, n2

    if (associated(stack%first)) then
      n1 => stack%first
      do i=1, stack%topIndex
        n2 => n1%next
        deallocate(n1)
        n1 => n2
      end do
    end if

    nullify(stack%first)
    nullify(stack%top)
  end subroutine deleteStack

  subroutine push(stack, newValue)
    type(NodeStack), intent(inout) :: stack
    type(QuadTreeNode), pointer, intent(in) :: newValue
    type(StackNode), pointer :: newNode

    allocate(newNode)
    newNode%data => newValue
    if (associated(stack%top)) then
      newNode%next => stack%top
    end if 
    if (.not.associated(stack%first)) then
      stack%first => newNode
    end if
    stack%top => newNode
    stack%topIndex = stack%topIndex + 1
  end subroutine push

  subroutine pop(stack, value)
    type(NodeStack), intent(inout) :: stack
    type(QuadTreeNode), pointer, intent(out) :: value
    type(StackNode), pointer :: tempNode

    if (associated(stack%top)) then
      value => stack%top%data
      tempNode => stack%top
      stack%top => stack%top%next
      deallocate(tempNode)
    else
      nullify(value)
    end if
    stack%topIndex = stack%topIndex - 1
  end subroutine pop

  !> initialize a new QuadTree
  subroutine initializeQuadTree(tree, maxDepth, maxElementsPerCell, width, height)
    implicit none
    integer, intent(in) :: maxDepth
    integer, intent(in) :: maxElementsPerCell
    real(dp), intent(in) :: width, height
  
    type(QuadTree), intent(out) :: tree

    tree%width = width
    tree%height = height
  
    allocate(tree%root)
    ! tree%root%x = 0._dp
    ! tree%root%y = 0._dp
    ! tree%root%width = width
    ! tree%root%height = height
    tree%root%numberOfElements = 0
    allocate(tree%root%elements(4*maxElementsPerCell))
    tree%root%elements(:) = -1
    nullify(tree%root%children)
  
  
    tree%maxDepth = maxDepth
    tree%maxElementsPerCell = maxElementsPerCell
    ! call initializeQuadTreeNode(tree%root, maxElementsPerCell, 0._dp, 0._dp, width, height)
  end subroutine initializeQuadTree

  !> add children to a node
  subroutine splitNode(node, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree

    real(dp) :: newWidth, newHeight, x,y
    integer :: i, j
    integer, dimension(4) :: idx

    allocate(node%children(4))
    idx(:) = 0

    newWidth = cellWidth(node) * tree%width / 2
    newHeight = cellHeight(node) * tree%height / 2
    x = cellX(node) * tree%width
    y = cellY(node) * tree%height

    ! x = node%x
    ! y = node%y
    ! print *, newWidth
    ! newWidth = node%width / 2
    ! newHeight = node%height / 2
    ! print *, newWidth

    do i = 1,4
      node%children(i)%numberOfElements = 0
      allocate(node%children(i)%elements(size(node%elements, 1)))
      node%children(i)%elements(:) = -1
      nullify(node%children(i)%children)
      ! node%children(i)%x = node%x
      ! node%children(i)%y = node%y
      ! node%children(i)%width = newWidth
      ! node%children(i)%height = newHeight
      
      node%children(i)%cellID = node%cellID + shiftl(i-1, node%cellLevel*2)
      node%children(i)%cellLevel = node%cellLevel + 1_i2
    end do
    ! node%children(2)%x = node%x + newWidth
    ! node%children(3)%y = node%y + newHeight
    ! node%children(4)%x = node%x + newWidth
    ! node%children(4)%y = node%y + newHeight



    ! call initializeQuadTreeNode(node%children(i), size(node%elements,1)/4, node%x, node%y, newWidth, newHeight)
    ! call initializeQuadTreeNode(childNode, size(node%elements,1)/4, node%x+newWidth, node%y, newWidth, newHeight)
    ! call initializeQuadTreeNode(childNode, size(node%elements,1)/4, node%x, node%y+newHeight, newWidth, newHeight)
    ! call initializeQuadTreeNode(childNode, size(node%elements,1)/4, node%x+newWidth, node%y+newHeight, newWidth, newHeight)

    do i = 0,node%numberOfElements-1
      !> find the child where to sort the element in
      if (node%elements(4*i+1) < x + newWidth) then
        !> left side
        j = 1
      else
        !> right side
        j = 2
      end if
      if (node%elements(4*i+2) > y + newHeight) then
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
    node%elements => null()

  end subroutine splitNode

  !> merge the children making node a leaf
  subroutine mergeChildren(node)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTreeNode), pointer :: childNode

    integer :: i,j, idx

    !> get the number of elements from a child node
    allocate(node%elements(size(node%children(1)%elements,1)))
    node%elements(:) = -1

    idx = 0

    do i=1,4
      childNode => node%children(i)
      do j=0,childNode%numberOfElements-1
        node%elements(4*idx+1:4*(idx+1)) = childNode%elements(4*j+1:4*(j+1))
        idx = idx+1
      end do
      deallocate(childNode%elements)
      ! deallocate(childNode)
    end do
    node%numberOfElements = idx
    deallocate(node%children)
  end subroutine mergeChildren


  recursive subroutine removeUnnecessaryNodes(node)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node

    type(QuadTreeNode), pointer :: childNode

    integer :: i,numElements

    logical :: isMergeable

    isMergeable = .true.

    if (node%numberOfElements >= 0) then
      return
    end if

    numElements = 0
    do i=1,4
      childNode => node%children(i)
      call removeUnnecessaryNodes(childNode)
      if (.not.childNode%isCollapsable .or. childNode%numberOfElements < 0) then
        isMergeable = .false.
      else 
        numElements = numElements + childNode%numberOfElements
      end if
    end do
    if (isMergeable .and. numElements <= size(node%children(1)%elements,1)/4) then
      call mergeChildren(node)
    end if
  end subroutine removeUnnecessaryNodes

  subroutine deleteTree(tree)
    implicit none
    type(QuadTree), intent(inout) :: tree

    call deleteSubTree(tree%root)

    deallocate(tree%root)

  end subroutine deleteTree

  !> delete a QuadTree by freeing all memory
  recursive subroutine deleteSubTree(node)
    implicit none
    ! type(QuadTreeNode), intent(inout) :: node
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTreeNode), pointer :: childNode
    integer :: i
   
    if (associated(node%children)) then
      do i = 1,4
        childNode => node%children(i)
        call deleteSubTree(childNode)
        ! call deleteSubTree(node%children(i))
      end do
      deallocate(node%children)
    end if
    ! print *, associated(node%elements)
    if (associated(node%elements)) then
      deallocate(node%elements)
    end if
    ! deallocate(node)
  end subroutine deleteSubTree

  recursive subroutine insertElement(node, element, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    real(dp), dimension(4), intent(in) :: element
    type(QuadTree), pointer, intent(inout) :: tree
    integer :: i
    type(QuadTreeNode), pointer :: nextNode => null()

    if (node%numberOfElements >= 0) then
      if (node%numberOfElements == size(node%elements, 1)/4) then
        !> if the cell is full, split it
        call splitNode(node, tree)
        call insertElement(node, element, tree)
      else
        !> otherwise insert the new element
        i = node%numberOfElements
        node%elements(4*i+1:4*(i+1)) = element(:)
        node%numberOfElements = i+1
      end if
    else
      if (element(1) < (cellX(node) + cellWidth(node)/2)*tree%width) then
        i = 1
      else
        i = 2
      end if
      if (element(2) >= (cellY(node) + cellHeight(node)/2)*tree%height) then
        i = i + 2
      end if
      
      nextNode => node%children(i)

      call insertElement(nextNode, element, tree)
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
    type(QuadTreeNode), pointer :: node
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
    ! print *, node%numberOfElements

  end subroutine removeElementFromNode

end module m_quadtree

