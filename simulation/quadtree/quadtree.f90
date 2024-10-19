module m_quadtree
  use m_types, only: fp, pp, i4, dp, i8, i1, i2
  use m_data_structures, only: ParticleAverageCounter, initializeParticleAverageCounter, deleteParticleAverageCounter
  use m_util, only: SimulationParameters, lineIntersection
  implicit none

  type :: CellStats
    type(ParticleAverageCounter), pointer :: particleCounter
    integer(i4), pointer :: N_avg !> average number of particles in the cell

    real(dp) :: maxSigmaC !> maxValue of sigma_T * c_r

  end type CellStats

  type :: QuadTreeNode
    integer(pp) :: numberOfElements

    !> encoding for the position of the cell
    !> 10 | 11
    !> -------
    !> 00 | 01
    !> binary representation of each level
    !> the first level is stored on the most right bits progressing to the left
    !> in the bits 58-63 the level of the cell is stored (bit 64 is used for the sign to mark invalid values)
    !> therefore the maximum level can be (8*8-6)/2=58/2=29
    integer(i8) :: cellID = 0

    !> store stats for each cell
    type(CellStats), pointer :: stats

    !> whether a grid cell is allowed to be merged
    !> false if it is given from the structure
    logical :: isCollapsable = .true.

    !> an structure is a line reaching from one end of the cell to another
    !> format: [left_x, left_y, right_x, right_y, sin(angle), cos(angle)]
    !> multiple lines are possible
    real(dp), pointer, dimension(:,:) :: structures => null()

    real(dp), pointer, dimension(:) :: elements => null()
    type(QuadTreeNode), pointer, dimension(:) :: children => null()
  end type QuadTreeNode

  type :: ListNode
    type(QuadTreeNode), pointer :: data
    type(ListNode), pointer :: next => null()
    type(ListNode), pointer :: previous => null()
  end type ListNode

  type :: NodeLinkedList
    integer :: length = 0
    type(ListNode), pointer :: first => null()
    type(ListNode), pointer :: last => null()
  end type NodeLinkedList

  type :: NodeLinkedListIterator
    type(NodeLinkedList), pointer :: list
    type(ListNode), pointer :: previous, current, next
  end type NodeLinkedListIterator


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
    type(NodeLinkedList), pointer :: leafs
    type(SimulationParameters) :: params

    type(QuadTreeNode), pointer :: root
  end type QuadTree

contains

  function cellLevel(node) result (level)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    integer(i1) :: level

    level = int(shiftr(node%cellID, 58), i1)
  end function cellLevel

  function cellWidth(node) result (width)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: width
    integer(i1) :: level

    level = cellLevel(node)
    
    width = 2._dp**(-level)
  end function cellWidth

  function cellHeight(node) result (height)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: height
    integer(i1) :: level

    level = cellLevel(node)
    
    height = 2._dp**(-level)
  end function cellHeight

  function cellX(node) result (x)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: width
    integer(i1) :: level, currentLevel
    real(dp) :: x

    width = 1._dp
    x = 0._dp

    level = cellLevel(node)

    do currentLevel = 0,level-1_i1
      width = width * .5
      x = x + width * ibits(node%cellID, currentLevel*2, 1)
    end do
  end function cellX
    
  function cellY(node) result (y)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(dp) :: height
    integer(i1) :: level, currentLevel
    real(dp) :: y

    height = 1._dp
    y = 0._dp

    level = cellLevel(node)

    do currentLevel = 0,level-1_i1
      height = height * .5
      y = y + height * ibits(node%cellID, currentLevel*2+1, 1)
    end do
  end function cellY

  !> return the cellID of a theoretical right neighbor cell of the same level
  !> if this is already a right most cell, return -1 as an error indication
  function rightNeighborID(node) result (cellID)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    integer(i8) :: cellID
    integer(i1) :: level, levelOfCell

    cellID = node%cellID
    levelOfCell = cellLevel(node)

    !> iterate starting from the lowest level
    do level = levelOfCell-1_i1,0,-1
      !> if the x-bit at the current level indicates a left cell, set it to the right cell
      !> and return the value
      if (ibits(node%cellID, 2*level, 1) == 0) then
        cellID = ibset(cellID, 2*level)
        return
      end if
      !> if the x-bit indicates a right cell, set it to a left cell and repeat with one level up
      cellID = ibclr(cellID, 2*level)
    end do
    !> if all levels have been right cells, there is no right neighbor cell
    !> set the return value to an error state
    cellID = -1_i1
  end function rightNeighborID

  !> return the cellID of a theoretical left neighbor cell of the same level
  !> if this is already a left most cell, return -1 as an error indication
  function leftNeighborID(node) result (cellID)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    integer(i8) :: cellID
    integer(i1) :: level, levelOfCell

    cellID = node%cellID
    levelOfCell = cellLevel(node)

    !> iterate starting from the lowest level
    do level = levelOfCell-1_i1,0,-1
      !> if the x-bit at the current level indicates a right cell, set it to the left cell
      !> and return the value
      if (ibits(node%cellID, 2*level, 1) == 1) then
        cellID = ibclr(cellID, 2*level)
        return
      end if
      !> if the x-bit indicates a left cell, set it to a right cell and repeat with one level up
      cellID = ibset(cellID, 2*level)
    end do
    !> if all levels have been left cells, there is no left neighbor cell
    !> set the return value to an error state
    cellID = -1_i1
  end function leftNeighborID

  !> return the cellID of a theoretical lower neighbor cell of the same level
  !> if this is already a lower most cell, return -1 as an error indication
  function lowerNeighborID(node) result (cellID)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    integer(i8) :: cellID
    integer(i1) :: level, levelOfCell

    cellID = node%cellID
    levelOfCell = cellLevel(node)

    !> iterate starting from the lowest level
    do level = levelOfCell-1_i1,0,-1
      !> if the x-bit at the current level indicates an upper cell, set it to the lower cell
      !> and return the value
      if (ibits(node%cellID, 2*level+1, 1) == 0) then
        cellID = ibset(cellID, 2*level+1)
        return
      end if
      !> if the x-bit indicates a lower cell, set it to a upper cell and repeat with one level up
      cellID = ibclr(cellID, 2*level+1)
    end do
    !> if all levels have been lower cells, there is no lower neighbor cell
    !> set the return value to an error state
    cellID = -1_i1
  end function lowerNeighborID

  !> return the cellID of a theoretical upper neighbor cell of the same level
  !> if this is already an upper most cell, return -1 as an error indication
  function upperNeighborID(node) result (cellID)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    integer(i8) :: cellID
    integer(i1) :: level, levelOfCell

    cellID = node%cellID
    levelOfCell = cellLevel(node)

    !> iterate starting from the lowest level
    do level = levelOfCell-1_i1,0,-1
      !> if the x-bit at the current level indicates a lower cell, set it to the upper cell
      !> and return the value
      if (ibits(node%cellID, 2*level+1, 1) == 1) then
        cellID = ibclr(cellID, 2*level+1)
        return
      end if
      !> if the x-bit indicates a left cell, set it to a lower cell and repeat with one level up
      cellID = ibset(cellID, 2*level+1)
    end do
    !> if all levels have been upper cells, there is no upper neighbor cell
    !> set the return value to an error state
    cellID = -1_i1
  end function upperNeighborID

  subroutine initializeLinkedList(list)
    implicit none
    type(NodeLinkedList), intent(inout) :: list

    nullify(list%first)
    nullify(list%last)
    list%length = 0
  end subroutine initializeLinkedList

  subroutine deleteLinkedList(list)
    implicit none
    type(NodeLinkedList), intent(inout) :: list

    type(ListNode), pointer :: n1,n2

    if (associated(list%first)) then
      n1 => list%first

      do while (associated(n1%next))
        n2 => n1%next
        deallocate(n1)
        n1 => n2
      end do
      deallocate(n1)
    end if

    nullify(list%first)
    nullify(list%last)
    list%length = 0
  end subroutine deleteLinkedList

  subroutine insertAfter(list, newValue, currentItem)
    implicit none
    type(NodeLinkedList), intent(inout) :: list
    type(QuadTreeNode), pointer, intent(in) :: newValue
    type(ListNode), pointer, intent(in), optional :: currentItem

    type(ListNode), pointer :: newItem

    allocate(newItem)
    newItem%data => newValue

    if (present(currentItem)) then
      if (associated(currentItem%next)) then
        newItem%next => currentItem%next
      end if
      currentItem%next => newItem
      newItem%previous => currentItem

      !> check if currentItem is the last item in the list and if yes update the list%last reference
      if (.not.associated(newItem%next)) then
        list%last => newItem
      end if
    else
      !> if no currentItem is given, append to the end of the list
      if (associated(list%last)) then
        list%last%next => newItem
      else
        !> if no item was in the list
        list%first => newItem
      end if
      list%last => newItem
    end if

  end subroutine insertAfter

  subroutine deleteListNode(list, item)
    implicit none
    type(NodeLinkedList), intent(inout) :: list
    type(ListNode), pointer, intent(inout) :: item

    if (associated(item%previous)) then
      item%previous%next => item%next
    else
      list%first => item%next
    end if
    if (associated(item%next)) then
      item%next%previous => item%previous
    else
      list%last => item%previous
    end if
    deallocate(item)

  end subroutine deleteListNode

  subroutine initializeLinkedListIterator(it, list)
    implicit none
    type(NodeLinkedListIterator), intent(inout) :: it
    type(NodeLinkedList), pointer, intent(in) :: list

    nullify(it%previous)
    nullify(it%current)
    if (associated(list%first)) then
      it%next => list%first
    else
      nullify(it%next)
    end if
  end subroutine initializeLinkedListIterator
      
  subroutine next(it)
    implicit none
    type(NodeLinkedListIterator), intent(inout) :: it

    if (.not.associated(it%current)) then
      print *, "reached end of list"
      return
    end if
    it%previous => it%current
    it%current => it%next
    if (associated(it%next)) then
      it%next => it%next%next
    end if 
  end subroutine next

  subroutine previous(it)
    implicit none
    type(NodeLinkedListIterator), intent(inout) :: it

    if (.not.associated(it%current)) then
      print *, "reached start of list"
      return
    end if
    it%next => it%current
    it%current => it%previous
    if (associated(it%previous)) then
      it%previous => it%previous%previous
    end if
  end subroutine previous


  subroutine getCurrentNode(it, node)
    implicit none

    type(NodeLinkedListIterator), intent(in) :: it
    type(ListNode), pointer, intent(out) :: node

    node => it%current

  end subroutine getCurrentNode



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


  subroutine initializeCellStats(stats, historyLength)
    implicit none
    type(CellStats), pointer, intent(inout) :: stats
    integer(i2), intent(in) :: historyLength

    allocate(stats)
    call initializeParticleAverageCounter(stats%particleCounter, historyLength)
    ! print *, allocated(stats%particleCounter%history), size(stats%particleCounter%history, 1)
    stats%maxSigmaC = 1e-16

  end subroutine initializeCellStats

  subroutine deleteCellStats(stats)
    implicit none
    type(CellStats), pointer, intent(inout):: stats

    call deleteParticleAverageCounter(stats%particleCounter)

    deallocate(stats)
  end subroutine deleteCellStats


  !> initialize a new QuadTree
  subroutine initializeQuadTree(tree, maxDepth, params)
    implicit none
    integer, intent(in) :: maxDepth
    type(SimulationParameters), intent(in) :: params
  
    type(QuadTree), pointer, intent(out) :: tree
    
    tree%params = params 

    tree%width = params%width
    tree%height = params%height
  
    allocate(tree%root)
    tree%root%cellID = 0_i8
    tree%root%numberOfElements = 0
    allocate(tree%root%elements(5*params%maxElementsPerCell))
    tree%root%elements(:) = -1
    nullify(tree%root%children)
    call initializeCellStats(tree%root%stats, params%cellHistoryLength)

    allocate(tree%leafs)
    call initializeLinkedList(tree%leafs)
    call insertAfter(tree%leafs, tree%root)

  
    tree%maxDepth = maxDepth
    tree%maxElementsPerCell = params%maxElementsPerCell
    ! call initializeQuadTreeNode(tree%root, maxElementsPerCell, 0._dp, 0._dp, width, height)
  end subroutine initializeQuadTree

  !> add children to a node
  subroutine splitNode(node, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree

    real(dp) :: newWidth, newHeight, x,y
    ! real(dp) :: obsMidX, obsMidY
    integer :: structureDirectionX, structureDirectionY
    real(dp), dimension(2) :: verticalIntersection, horizontalIntersection
    real(dp), dimension(:, :, :), allocatable :: tmpStructures
    integer, dimension(4) :: numStructures
    integer :: i, j
    integer, dimension(4) :: idx
    integer(i8) :: level, childLevel

    allocate(node%children(4))

    newWidth = cellWidth(node) * tree%width / 2
    newHeight = cellHeight(node) * tree%height / 2
    x = cellX(node) * tree%width
    y = cellY(node) * tree%height

    ! print *, newHeight, newWidth

    level = cellLevel(node)
    childLevel = ibits(node%cellID, 58, 6) + 1_i8

    do i = 1,4
      node%children(i)%numberOfElements = 0
      allocate(node%children(i)%elements(size(node%elements, 1)))
      node%children(i)%elements(:) = -1
      nullify(node%children(i)%children)

      node%children(i)%cellID = node%cellID + shiftl(i-1, level*2)
      call mvbits(childLevel, 0, 6, node%children(i)%cellID, 58)

      call initializeCellStats(node%children(i)%stats, node%stats%particleCounter%historyLength)
      node%children(i)%stats%particleCounter%history(:) = node%stats%particleCounter%history(:) / 4
      node%children(i)%stats%particleCounter%average = node%stats%particleCounter%average / 4
      node%children(i)%stats%particleCounter%currentPosition = node%stats%particleCounter%currentPosition
    end do

    ! !> split the structures
    ! if (associated(node%structures)) then
    !   do i = 0, size(node%structures/5)-1
    !     !> find the y coordinate of the crossing
    !     obsMidY = (node%structures(5*i+4)-node%structures(5*i+2))/(node%structures(5*i+3)-node%structures(5*i+1)) &
    !       * (x+newWidth - node%structures(5*i+1)) + node%structures(5*i+2)
    !     obsMidX = (node%structures(5*i+3)-node%structures(5*i+1))/(node%structures(5*i+4)-node%structures(5*i+2)) &
    !       * (y+newHeight - node%structures(5*i+2)) + node%structures(5*i+4)
    !  
    !     if (.not.all(obsMidX == -1))
    !   
    !     end if
    !   end do
    ! end if

    idx(:) = 1
    if (associated(node%structures)) then
      allocate(tmpStructures(4, size(node%structures,1), 6))
      numStructures = 0
      do i = 1, size(node%structures,1)
    
        !> store the orientation of the structure line to find the corresponding subcells later
        if (node%structures(i, 1) - node%structures(i, 3) < 0._dp) then
          structureDirectionX = 1
        else
          structureDirectionX = -1
        end if
        if (node%structures(i, 2) - node%structures(i, 4) < 0._dp) then
          structureDirectionY = 1
        else
          structureDirectionY = -1
        end if
        
        ! !> in the case that a structure is exactly aligning with the middle
        ! !> store the structure for both halfs
        ! if (all(abs(node%structure(:4) - [x+newWidth, y, x+newWidth, y+newHeight*2]) < 1e-10)) then
        !   !> if it aligns with the vertical line
        !   do j = 1,2
        !   tmpStructures(j, idx(j), :) = 
    
        verticalIntersection = lineIntersection([x+newWidth, y, x+newWidth, y+newHeight*2], node%structures(i, :4))
        horizontalIntersection = lineIntersection([x, y+newHeight, x+newWidth*2, y+newHeight], node%structures(i, :4))
        ! print *, node%structures(i, :4)
        ! print *, [x+newWidth, y, x+newWidth, y+newHeight*2]
        ! print *, [x, y+newHeight, x+newWidth*2, y+newHeight]
        ! print *, verticalIntersection(:)
        ! print *, horizontalIntersection(:)
    
        !> determine the subcells which the line crosses
        if (node%structures(i, 1) >= x + newWidth) then
          j = 2
        else
          j = 1
        end if
        if (node%structures(i,2) >= y + newHeight) then
          j = j + 2
        end if
        tmpStructures(j, idx(j), 1:2) = node%structures(i,1:2)
    
        if (any(verticalIntersection /= -1._dp)) then
          if (any(horizontalIntersection /= -1._dp)) then
            if ((verticalIntersection(1)-horizontalIntersection(1))*structureDirectionX < 0) then
              !> if the vertical intersection happens first
              tmpStructures(j,idx(j),3:4) = verticalIntersection
              tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
              idx(j) = idx(j) + 1
              j = j + structureDirectionX
              tmpStructures(j,idx(j),1:2) = verticalIntersection
    
              tmpStructures(j,idx(j),3:4) = horizontalIntersection
              tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
              idx(j) = idx(j) + 1
              j = j + 2*structureDirectionY
              tmpStructures(j,idx(j),1:2) = horizontalIntersection
    
            else
              !> if the horizontalIntersection happens first
              tmpStructures(j,idx(j),3:4) = horizontalIntersection
              tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
              idx(j) = idx(j) + 1
              j = j + 2*structureDirectionY
              tmpStructures(j,idx(j),1:2) = horizontalIntersection
             
              tmpStructures(j,idx(j),3:4) = verticalIntersection
              tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
              idx(j) = idx(j) + 1
              j = j + structureDirectionX
              tmpStructures(j,idx(j),1:2) = verticalIntersection
            end if
          else
            !> only the verticalIntersection happens
            tmpStructures(j,idx(j),3:4) = verticalIntersection
            tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
            idx(j) = idx(j) + 1
            j = j + structureDirectionX
            tmpStructures(j,idx(j),1:2) = verticalIntersection
          end if
        elseif (any(horizontalIntersection /= -1._dp)) then
          !> only the horizontalIntersection happens
          tmpStructures(j,idx(j),3:4) = horizontalIntersection
          tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
          idx(j) = idx(j) + 1
          j = j + 2*structureDirectionY
          tmpStructures(j,idx(j),1:2) = horizontalIntersection
        end if
        tmpStructures(j,idx(j),3:4) = node%structures(i, 3:4)
        tmpStructures(j,idx(j),5:6) = node%structures(i,5:6)
        idx(j) = idx(j) + 1
    
        ! if (any(verticalIntersection /= -1._dp)) then
        !   if (any(horizontalIntersection /= -1._dp)) then
        !     if (node%structures(i))
        !
        !
        !     if (verticalIntersection(2) > x + newHeight) then
        !
        !
        !     end if
        !
        !   else
        !
        !   end if
        ! elseif (any(horizontalIntersection /= -1)) then
        !
        ! else
        !
        ! end if
      end do
    
      do j = 1,4
        if (idx(j)-1 > 0) then
          allocate(node%children(j)%structures(idx(j)-1,6))
          node%children(j)%structures(:idx(j)-1, :) = tmpStructures(j, :idx(j)-1, :)
        end if
      end do
      ! print *, tmpStructures
    end if


    idx(:) = 0
    do i = 0,node%numberOfElements-1
      !> find the child where to sort the element in
      if (node%elements(5*i+1) < x + newWidth) then
        !> left side
        j = 1
      else
        !> right side
        j = 2
      end if
      if (node%elements(5*i+2) > y + newHeight) then
        !> lower row
        j = j+2
      end if
      !> copy the element entries to the corresponding child
      node%children(j)%elements(5*idx(j)+1:5*(idx(j)+1)) = node%elements(5*i+1:5*(i+1))
      idx(j) = idx(j)+1
    end do

    do i=1,4
      node%children(i)%numberOfElements = idx(i)
    end do
    
    !> mark that the node is no leave
    node%numberOfElements= -1
    deallocate(node%elements)
    node%elements => null()
    deallocate(node%stats)
    node%stats => null()

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

    call initializeCellStats(node%stats, node%children(1)%stats%particleCounter%historyLength)
    node%stats%particleCounter%currentPosition = node%children(1)%stats%particleCounter%currentPosition

    idx = 0

    do i=1,4
      childNode => node%children(i)
      node%stats%particleCounter%history(:) = node%stats%particleCounter%history(:) + childNode%stats%particleCounter%history(:)
      node%stats%particleCounter%average = node%stats%particleCounter%average + childNode%stats%particleCounter%average
      do j=0,childNode%numberOfElements-1
        node%elements(5*idx+1:5*(idx+1)) = childNode%elements(5*j+1:5*(j+1))
        idx = idx+1
      end do
      deallocate(childNode%elements)
      call deleteCellStats(childNode%stats)
      ! deallocate(childNode%stats)
      if (associated(childNode%structures)) then
        deallocate(childNode%structures)
      end if
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
    if (isMergeable .and. numElements <= size(node%children(1)%elements,1)/5) then
      call mergeChildren(node)
    end if
  end subroutine removeUnnecessaryNodes

  subroutine deleteTree(tree)
    implicit none
    type(QuadTree), intent(inout) :: tree

    call deleteLinkedList(tree%leafs)
    deallocate(tree%leafs)
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

  subroutine insertElementWithCellIDHint(cellID, element, tree)
    implicit none
    integer(i8), intent(in) :: cellID
    real(dp), dimension(5), intent(in) :: element
    type(QuadTree), pointer, intent(inout) :: tree

    type(QuadTreeNode), pointer :: node
    integer(i1) :: level, levelOfCell
    real(dp) :: x,y,width,height

    node => tree%root

    levelOfCell = int(shiftr(cellID, 58), i1)

    do level = 0,levelOfCell
      if (associated(node%children)) then
        node => node%children(ibits(cellID, 2*level, 2))
      else
        exit
      end if
    end do
    x = cellX(node) * tree%width
    y = celly(node) * tree%height
    width = cellWidth(node) * tree%width
    height = cellHeight(node) * tree%height
    if (x <= element(1) .and. element(1) < x+width &
        .and. y <= element(2) .and. element(2) < y+height) then
      call insertElement(node, element, tree)
    else
      !> if the hint was not correct, use a usual insert starting from the root node
      call insertElement(tree%root, element, tree)
    end if

  end subroutine insertElementWithCellIDHint

  recursive subroutine insertElement(node, element, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    real(dp), dimension(5), intent(in) :: element
    type(QuadTree), pointer, intent(inout) :: tree
    integer :: i
    type(QuadTreeNode), pointer :: nextNode => null()

    ! print *, node%numberOfElements
    ! write (*, "(B64.64)") node%cellID
    ! print *, size(node%elements, 1)

    if (node%numberOfElements >= 0) then
      if (node%numberOfElements == size(node%elements, 1)/5) then
        !> if the cell is full, split it
        call splitNode(node, tree)
        call insertElement(node, element, tree)
      else
        !> otherwise insert the new element
        i = node%numberOfElements
        node%elements(5*i+1:5*(i+1)) = element(:)
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
    node%elements(5*(i-1)+1:5*i) = node%elements(5*(j-1)+1:5*j)
    !> remove the last element
    node%elements(5*(j-1)+1:5*j) = -1
    node%numberOfElements = j - 1
    ! print *, node%numberOfElements

  end subroutine removeElementFromNode

end module m_quadtree

