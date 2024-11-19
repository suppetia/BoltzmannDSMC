module m_quadtree
  use m_types, only: fp, i1, i4, i8
  use m_datastructures, only: QuadTreeParameters, SimulationParameters, StructuresPointer, CellStats, &
    addParticleCount, initializeCellStats, deleteCellStats
  use m_util, only: lineIntersection, sort!,conditionedMergeSort

  use omp_lib

  implicit none

  type QuadTreeNode

    !> encoding for the position of the cell
    !> 10 | 11
    !> -------
    !> 00 | 01
    !> binary representation of each level
    !> the first level is stored on the most right bits progressing to the left
    !> in the bits 58-63 the level of the cell is stored (bit 64 is used for the sign to mark invalid values)
    !> therefore the maximum level can be (8*8-6)/2=58/2=29
    integer(i8) :: cellID = 0

    !> index of the node in the data arrays
    !> -1 if not a leaf
    integer(i4) :: nodeIdx = -1

    !> stats for each node
    type(CellStats), pointer :: stats => null()

    !> the structures are stored in a 2D-array for each node
    !> the structure array just stores a pointer to that array
    !> ---
    !> a structure is a line reaching from one end of the cell to another
    !> format: [left_x, left_y, right_x, right_y, n_x, n_y], where n is the (unit) normal vector
    !> multiple lines are possible
    real(fp), pointer, dimension(:,:) :: structures => null()

    !> whether a grid cell is allowed to be merged
    !> false if it is given from the structure
    logical :: isCollapsable = .true.
    
    !> child nodes (if present)
    type(QuadTreeNode), pointer, dimension(:) :: children => null()
    !> parent node (if present)
    type(QuadTreeNode), pointer :: parent => null()
 
    !> TODO: check if there is another way to do this
    !> when splitting a cell, store the number of elements for the subcells
    !> negative values can be used as flags: -1: no recent updates, -2: freshly split 
    integer(i4) :: tmpParticleCount = -1
    !> index of the particles associated with this node starting in the tree particle array
    integer(i4) :: tmpParticleIndex = -1
    !> array of particles belonging to the cell
    real(fp), pointer, dimension(:,:) :: tmpParticles => null()
    integer(i1), pointer, dimension(:) :: tmpParticleTypes => null()

  end type QuadTreeNode

  !> wrapper type to store the pointer to a QuadTreeNodePointer
  type QuadTreeNodePointer
    type(QuadTreeNode), pointer :: node => null()
  end type QuadTreeNodePointer

  type QuadTree

    !> metadata of the QuadTree
    type(QuadTreeParameters), pointer :: treeParams => null()

    !> root node
    type(QuadTreeNode), pointer :: root => null()

    !> array where the pointer to the corresponding leaf node is stored
    !> TODO: maybe change to the cellID ??
    type(QuadTreeNodePointer), pointer, dimension(:) :: leafs => null()

    !> number of the leafs
    !> this is modified during the child merge process and hence temporarily does not equal the length of the leafs array
    integer(i4) :: leafNumber = 0

    !> array of all the particles
    !> dimensions: (N + overhead) x 5
    real(fp), pointer, dimension(:,:) :: particles => null()
    !> array of the types of particles
    integer(i1), pointer, dimension(:) :: particleTypes => null()
    !> array where the number of particles of each leaf node is stored
    integer(i4), pointer, dimension(:) :: particleNumbers => null()
    !> array where the starting index of the memory of the particles for each leaf is stored
    integer(i4), pointer, dimension(:) :: particleStartIndices => null()

    !> array of pointers associated with the leaf node's structures
    type(StructuresPointer), pointer, dimension(:) :: structures => null()

  end type QuadTree

  !> data structure to temporarily store references to nodes
  !> elements of such stack
  type :: StackNode
    type(QuadTreeNode), pointer :: data
    type(StackNode), pointer :: next => null()
  end type StackNode

  !> the stack itself
  type :: NodeStack
    integer :: topIndex = 0
    type(StackNode), pointer :: first => null()
    type(StackNode), pointer :: top => null()
  end type NodeStack


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
    real(fp) :: width
    integer(i1) :: level

    level = cellLevel(node)
    
    width = 2._fp**(-level)
  end function cellWidth

  function cellHeight(node) result (height)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(fp) :: height
    integer(i1) :: level

    level = cellLevel(node)
    
    height = 2._fp**(-level)
  end function cellHeight

  function cellX(node) result (x)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(fp) :: width
    integer(i1) :: level, currentLevel
    real(fp) :: x

    width = 1._fp
    x = 0._fp

    level = cellLevel(node)

    do currentLevel = 0,level-1_i1
      width = width * .5
      x = x + width * ibits(node%cellID, currentLevel*2, 1)
    end do
  end function cellX
    
  function cellY(node) result (y)
    implicit none
    type(QuadTreeNode), pointer, intent(in) :: node
    real(fp) :: height
    integer(i1) :: level, currentLevel
    real(fp) :: y

    height = 1._fp
    y = 0._fp

    level = cellLevel(node)

    do currentLevel = 0,level-1_i1
      height = height * .5
      y = y + height * ibits(node%cellID, currentLevel*2+1, 1)
    end do
  end function cellY

  subroutine initializeTree(tree, params)
    implicit none
    type(QuadTree), pointer, intent(out) :: tree
    type(QuadTreeParameters), pointer, intent(inout) :: params
    
    allocate(tree)
    tree%treeParams => params

    allocate(tree%root)
    call initializeCellStats(tree%root%stats, params)
    tree%root%nodeIdx = 1
    tree%leafNumber = 1
    allocate(tree%leafs(tree%leafNumber))
    tree%leafs(1)%node => tree%root
    allocate(tree%particles(tree%treeParams%elementChunkSize, 5))
    allocate(tree%particleTypes(tree%treeParams%elementChunkSize))
    allocate(tree%particleNumbers(tree%leafNumber))
    allocate(tree%particleStartIndices(tree%leafNumber))
    allocate(tree%structures(1))
    tree%particleNumbers(1) = 0
    tree%particleStartIndices(1) = 1
  end subroutine initializeTree 


  recursive subroutine splitNode(node, tree, stack)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(NodeStack), intent(inout) :: stack

    real(fp), pointer, dimension(:,:) :: tmpParticles
    integer(i1), pointer, dimension(:) :: tmpParticleTypes
    integer(i4) :: numParticles, idxStartParticles
    integer(i1), dimension(:), pointer :: particleLocations
    integer(i4), dimension(4) :: numParticlesPerChildCell

    real(fp), dimension(:,:,:), allocatable :: tmpStructures
    real(fp), dimension(2) :: verticalIntersection, horizontalIntersection
    integer(i1) :: structureDirectionX, structureDirectionY

    real(fp) :: newWidth, newHeight, x,y
    integer(i8) :: level

    integer(i1) :: j
    integer(i4) :: i
    integer(i4), dimension(4) :: idx

    type(QuadTreeNode), pointer :: childNode

    allocate(node%children(4))

    !> update the number of leafs in the tree
    tree%leafNumber = tree%leafNumber + 3

    newWidth = cellWidth(node) * tree%treeParams%width / 2
    newHeight = cellHeight(node) * tree%treeParams%height / 2
    x = cellX(node) * tree%treeParams%width
    y = cellY(node) * tree%treeParams%height

    level = ibits(node%cellID, 58, 6)

    if (node%tmpParticleCount < 0) then
      !> if the node has not been recently created by splitting another node
      numParticles = tree%particleNumbers(node%nodeIdx)
      idxStartParticles = tree%particleStartIndices(node%nodeIdx)
    else
      !> if the node is "freshly" created, the tree was not updated yet
      !> in that case take the number of particles from the tmp variable
      numParticles = node%tmpParticleCount
      idxStartParticles = node%tmpParticleIndex
    end if

    !> the node is not a leaf anymore
    node%nodeIdx = -1


    !> initialize the child nodes
    do i = 1,4
      node%children(i)%cellID = node%cellID + shiftl(i-1, level*2)
      !> update the cellID by the cellLevel
      call mvbits(level + 1_i8, 0, 6, node%children(i)%cellID, 58)
      !> store the information of the parent cell
      node%children(i)%parent => node

      !> split the cellstats evenly for all child nodes
      call initializeCellStats(node%children(i)%stats, tree%treeParams)
      node%children(i)%stats%particleCounter%history = node%stats%particleCounter%history / 4
      node%children(i)%stats%particleCounter%average = node%stats%particleCounter%average / 4
      node%children(i)%stats%particleCounter%currentPosition = node%stats%particleCounter%currentPosition
    end do
    call deleteCellStats(node%stats)

    !> split structures if present
    if (associated(node%structures)) then
      idx(:) = 1
      allocate(tmpStructures(4, 6, size(node%structures,2)))
      do i = 1, size(node%structures,2)
    
        !> store the orientation of the structure line to find the corresponding subcells later
        if (node%structures(1, i) - node%structures(3, i) < 0._fp) then
          structureDirectionX = 1_i1
        else
          structureDirectionX = -1_i1
        end if
        if (node%structures(2, i) - node%structures(4, i) < 0._fp) then
          structureDirectionY = 1_i1
        else
          structureDirectionY = -1_i1
        end if
    
        verticalIntersection = lineIntersection([x+newWidth, y, x+newWidth, y+newHeight*2], node%structures(:4,i))
        horizontalIntersection = lineIntersection([x, y+newHeight, x+newWidth*2, y+newHeight], node%structures(:4,i))

        !> determine the subcells which the line crosses
        if (node%structures(1,i) >= x + newWidth) then
          j = 2_i1
        else
          j = 1_i1
        end if
        if (node%structures(2,i) >= y + newHeight) then
          j = j + 2_i1
        end if
        tmpStructures(j, 1:2, idx(j)) = node%structures(1:2,i)
    
        if (any(verticalIntersection /= -1._fp)) then
          if (any(horizontalIntersection /= -1._fp)) then
            !> TODO: handle the case where the structure is exactly aligning with the middle
            if ((verticalIntersection(1)-horizontalIntersection(1))*structureDirectionX < 0) then
              !> if the vertical intersection happens first
              tmpStructures(j,3:4,idx(j)) = verticalIntersection
              tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
              idx(j) = idx(j) + 1
              j = j + structureDirectionX
              tmpStructures(j,1:2,idx(j)) = verticalIntersection
    
              tmpStructures(j,3:4,idx(j)) = horizontalIntersection
              tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
              idx(j) = idx(j) + 1
              j = j + 2_i1 * structureDirectionY
              tmpStructures(j,1:2,idx(j)) = horizontalIntersection
    
            else
              !> if the horizontalIntersection happens first
              tmpStructures(j,3:4,idx(j)) = horizontalIntersection
              tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
              idx(j) = idx(j) + 1
              j = j + 2_i1 * structureDirectionY
              tmpStructures(j,1:2,idx(j)) = horizontalIntersection
             
              tmpStructures(j,3:4,idx(j)) = verticalIntersection
              tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
              idx(j) = idx(j) + 1
              j = j + structureDirectionX
              tmpStructures(j,1:2,idx(j)) = verticalIntersection
            end if
          else
            !> only the verticalIntersection happens
            tmpStructures(j,3:4,idx(j)) = verticalIntersection
            tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
            idx(j) = idx(j) + 1
            j = j + structureDirectionX
            tmpStructures(j,1:2,idx(j)) = verticalIntersection
          end if
        elseif (any(horizontalIntersection /= -1._fp)) then
          !> only the horizontalIntersection happens
          tmpStructures(j,3:4,idx(j)) = horizontalIntersection
          tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
          idx(j) = idx(j) + 1
          j = j + 2_i1 * structureDirectionY
          tmpStructures(j,1:2,idx(j)) = horizontalIntersection
        end if
        tmpStructures(j,3:4,idx(j)) = node%structures(3:4,i)
        tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
        idx(j) = idx(j) + 1
    
      end do
    
      do j = 1_i1,4_i1
        if (idx(j)-1 > 0) then
          allocate(node%children(j)%structures(6, idx(j)-1))
          node%children(j)%structures(:, :idx(j)-1) = tmpStructures(j, :, :idx(j)-1)
        end if
      end do
    end if

    !> handle the particles
    if (numParticles == 0) then
      do j = 1_i1, 4_i1
        node%children(j)%tmpParticleCount = 0
      end do 
      return
    end if
    
    !> sort the particles into the childnodes
    numParticlesPerChildCell = 0

    allocate(particleLocations(numParticles))
    allocate(tmpParticles(numParticles, 5))
    ! allocate(tmpParticleTypes(numParticles))

    ! !> copy the particles to a temporary array for fast access and easy sorting
    ! tmpParticles = tree%particles(idxStartParticles:idxStartParticles+numParticles-1, :)
    ! tmpParticleTypes = tree%particleTypes(idxStartParticles:idxStartParticles+numParticles-1)
    !> copy the particles to a temporary array for fast access and easy sorting
    !> TODO: check whether this is faster on the original arrays
    ! print *, node%tmpParticleCount, associated(node%tmpParticles)
    if (node%tmpParticleCount < 0) then
      tmpParticles = tree%particles(idxStartParticles:idxStartParticles+numParticles-1, :)
      ! tmpParticleTypes = tree%particleTypes(idxStartParticles:idxStartParticles+numParticles-1)
    else
      tmpParticles = node%tmpParticles
      ! tmpParticleTypes = node%tmpParticleTypes
      deallocate(node%tmpParticles)
      ! deallocate(node%tmpParticleTypes)
      node%tmpParticleCount = -1
    end if
    !> mark node as freshly split
    ! node%tmpParticleCount = -1

    !> find out which particles go into which child node
    do i = 1, numParticles
      !> find the child where to sort the element in
      if (tmpParticles(i, 1) < x + newWidth) then
        !> left side
        j = 1_i1
      else
        !> right side
        j = 2_i1
      end if
      if (tmpParticles(i, 2) > y + newHeight) then
        !> lower row
        j = j+2_i1
      end if
      particleLocations(i) = j
      numParticlesPerChildCell(j) = numParticlesPerChildCell(j) + 1
    end do

    !> find the starting indices for each sub node
    ! idx = idxStartParticles
    idx = 1
    do j=1,4
      ! idx(:i) = idx(:i) + numParticlesPerChildCell(i)
      allocate(node%children(j)%tmpParticles(numParticlesPerChildCell(j), 5))
      node%children(j)%tmpParticleCount = numParticlesPerChildCell(j)
    end do
      

    !> actually sort the particles into the tree array
    !> TODO: can be parallelized
    do i = 1, numParticles
      ! tree%particles(:, idx(particleLocations(i))) = tmpParticles(:, i)
      ! tree%particleTypes(idx(particleLocations(i))) = tmpParticleTypes(i)

      node%children(particleLocations(i))%tmpParticles(idx(particleLocations(i)),:) = tmpParticles(i, :)

      idx(particleLocations(i)) = idx(particleLocations(i)) + 1
    end do
    deallocate(tmpParticles)
    ! deallocate(tmpParticleTypes)
    deallocate(particleLocations)

    !> check if all the child cells fulfill the threshold condition
    if (level < 63) then
      do j = 1,4
        childNode => node%children(j)
        if (childNode%tmpParticleCount > tree%treeParams%elementSplitThreshold) then
          call splitNode(childNode, tree, stack)
        ! else
        !   call push(stack, childNode)
        end if
      end do
    end if 
    
  end subroutine splitNode

  !> split all nodes that exceed the threshold defined in the tree parameters
  subroutine splitNodes(tree, stack)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(NodeStack), intent(inout) :: stack

    ! type(QuadTreeNode), pointer :: node, childNode
    integer(i4) :: i

    !> TODO: can be parallelized
    do i = 1, size(tree%leafs)
      !> use node%tmpParticleCount as this metric is updated when inserting particles
      if (max(tree%particleNumbers(i),tree%leafs(i)%node%tmpParticleCount) &
        >= tree%treeParams%elementSplitThreshold) then
        call splitNode(tree%leafs(i)%node, tree, stack)
        nullify(tree%leafs(i)%node)
      end if
    end do
  end subroutine splitNodes

  !> merge the child nodes making node a leaf
  subroutine mergeChildNodes(node, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(QuadTreeNode), pointer :: childNode

    integer(i1) :: i
    integer(i4) :: j, idx, numParticles

    !> update the number of leafs in the tree
    tree%leafNumber = tree%leafNumber - 3

    !> get the number of elements from the child nodes
    numParticles = 0
    !> TODO: can be parallelized
    do i = 1,4
      if (node%children(i)%tmpParticleCount < 0) then
        numParticles = numParticles + tree%particleNumbers(node%children(i)%nodeIdx)
      else
        numParticles = numParticles + node%children(i)%tmpParticleCount
      end if
    end do
    node%tmpParticleCount = numParticles

    allocate(node%tmpParticles(numParticles, 5))

    !> merge the cell stats from the child nodes
    call initializeCellStats(node%stats, tree%treeParams)
    node%stats%particleCounter%currentPosition = node%children(1)%stats%particleCounter%currentPosition

    idx = 1

    do i=1,4
      childNode => node%children(i)
      node%stats%particleCounter%history = node%stats%particleCounter%history + childNode%stats%particleCounter%history
      node%stats%particleCounter%average = node%stats%particleCounter%average + childNode%stats%particleCounter%average

      if (childNode%tmpParticleCount < 0) then !> if the node was a leaf in the last iteration
        numParticles = tree%particleNumbers(childNode%nodeIdx)
        j = tree%particleStartIndices(childNode%nodeIdx)
        node%tmpParticles(idx:idx+numParticles-1, :) = tree%particles(j:j+numParticles-1, :)

        !> remove the pointer to this childNode from the tree
        tree%leafs(childNode%nodeIdx)%node => null()
      else
        numParticles = childNode%tmpParticleCount
        node%tmpParticles(idx:idx+numParticles-1, :) = childNode%tmpParticles
        deallocate(childNode%tmpParticles)
        ! if (associated(childNode%tmpParticles)) then
        !   if (numParticles > 0) then
        !     node%tmpParticles(:, idx:idx+numParticles-1) = childNode%tmpParticles
        !   end if 
        !   deallocate(childNode%tmpParticles)
        ! end if
      end if
      idx = idx + numParticles
      nullify(childNode%parent)
      ! call deleteCellStats(childNode%stats)
      if (associated(childNode%structures)) then
        deallocate(childNode%structures)
      end if
    end do
    deallocate(node%children)
  end subroutine mergeChildNodes

  subroutine removeUnnecessaryNodes(tree, newLeafsStack)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(NodeStack), intent(inout) :: newLeafsStack

    integer(i4) :: i
    type(QuadTreeNode), pointer :: node

    do i = 1, size(tree%leafs)
      if (associated(tree%leafs(i)%node)) then
        if (associated(tree%leafs(i)%node%parent)) then
          node => tree%leafs(i)%node%parent
          call removeUnnecessaryUpwardNodes(node, tree, newLeafsStack)
        end if
      end if
    end do
    ! call findNewLeafsAfterRemoval(tree%root, newLeafsStack)

  end subroutine removeUnnecessaryNodes

  recursive subroutine removeUnnecessaryNodesDownwards(node, tree, newLeafsStack)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(NodeStack), intent(inout) :: newLeafsStack

    type(QuadTreeNode), pointer :: childNode
    integer(i4) :: i, numElements
    logical :: isMergeable

    if (.not.associated(node%children)) then
      return
    end if

    isMergeable = .true.

    numElements = 0
    do i=1,4
      childNode => node%children(i)
      call removeUnnecessaryNodesDownwards(childNode, tree, newLeafsStack)
      if (.not.childNode%isCollapsable .or. associated(childNode%children)) then
        isMergeable = .false.
      else 
        if (childNode%tmpParticleCount < 0) then
          numElements = numElements + tree%particleNumbers(childNode%nodeIdx)
        else
          numElements = numElements + childNode%tmpParticleCount
        end if
      end if
    end do
    if (isMergeable) then
      if (numElements <= tree%treeParams%elementMergeThreshold) then
        call mergeChildNodes(node, tree)
      end if
    end if
    if (associated(node%children)) then
      do i=1,4
        childNode => node%children(i)
        if (.not.associated(node%children)) then
          call push(newLeafsStack, childNode)
        end if
      end do
    end if
  end subroutine removeUnnecessaryNodesDownwards

  recursive subroutine findLeafs(n, newLeafsStack)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: n
    type(NodeStack), intent(inout) :: newLeafsStack

    type(QuadTreeNode), pointer :: childNode
    integer(i1) :: i

    if (associated(n%children)) then
      !> TODO: can be parallelized
      do i = 1,4
        childNode => n%children(i)
        call findLeafs(childNode, newLeafsStack)
      end do
    else
      call push(newLeafsStack, n)
    end if
  end subroutine findLeafs

  ! recursive subroutine findNewLeafsAfterRemoval(n, newLeafsStack)
  !   implicit none
  !   type(QuadTreeNode), pointer, intent(inout) :: n
  !   type(NodeStack), intent(inout) :: newLeafsStack
  !  
  !   type(QuadTreeNode), pointer :: childNode
  !   integer(i1) :: i
  !
  !   if (n%tmpParticleCount == -2) then
  !     call markChildNodesAsProcessed(n)
  !   else if (.not.associated(n%children)) then
  !     call push(newLeafsStack, n)
  !   else
  !     do i = 1,4
  !       childNode => n%children(i)
  !       call findNewLeafsAfterRemoval(childNode, newLeafsStack)
  !     end do
  !   end if
  ! end subroutine findNewLeafsAfterRemoval
  !
  ! recursive subroutine markChildNodesAsProcessed(node)
  !   implicit none
  !     n%tmpParticleCount = -1
  !   type(QuadTreeNode), pointer, intent(inout) :: node
  !
  !   integer(i1) :: i
  !   type(QuadTreeNode), pointer :: childNode
  !
  !   if (node%tmpParticleCount == -2) then
  !     node%tmpParticleCount = -1
  !
  !     do i = 1,4
  !       childNode => node%children(i)
  !       call markChildNodesAsProcessed(childNode)
  !     end do
  !   end if
  ! end subroutine markChildNodesAsProcessed

  !> traverse from the bottom to the top and merge all possible nodes
  recursive subroutine removeUnnecessaryUpwardNodes(node, tree, newLeafsStack)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(NodeStack), intent(inout) :: newLeafsStack
   
    integer(i1) :: i
    logical :: canBeRemoved
    integer(i4) :: numElements 
    type(QuadTreeNode), pointer :: childNode
  
  
    canBeRemoved = .true.
  
    numElements = 0
    do i=1,4
      childNode => node%children(i)
      if (.not.childNode%isCollapsable .or. associated(childNode%children)) then
        canBeRemoved = .false.
        exit
      end if
      if (childNode%tmpParticleCount < 0) then
        numElements = numElements + tree%particleNumbers(childNode%nodeIdx)
      else
        numElements = numElements + childNode%tmpParticleCount
      end if
    end do
    if (canBeRemoved) then
      if (numElements <= tree%treeParams%elementMergeThreshold) then
        call mergeChildNodes(node, tree)
        ! call push(newLeafsStack, node)
        !> check if the parent node still has too few elements
        if (associated(node%parent)) then
          call removeUnnecessaryUpwardNodes(node%parent, tree, newLeafsStack)
        ! else
        !   !> if it is the root node
        !   call push(newLeafsStack, node)
        end if
        return
      end if
    end if
    ! !> otherwise check if the child nodes are new leafs
    ! do i=1,4
    !   childNode => node%children(i)
    !   if (.not.associated(childNode%children) .and. childNode%nodeIdx == -1) then
    !     call push(newLeafsStack, childNode)
    !   end if
    ! end do
  end subroutine removeUnnecessaryUpwardNodes

  ! !> find the mergeable nodes and remove them
  ! recursive subroutine removeUnnecessaryDownwardNodes(node, tree, newLeafsStack)
  !   implicit none
  !   type(QuadTreeNode), pointer, intent(inout) :: node
  !   type(QuadTree), pointer, intent(inout) :: tree
  !   type(NodeStack), intent(inout) :: newLeafsStack
  !
  !   type(QuadTreeNode), pointer :: childNode
  !
  !   integer(i4) :: i,numElements
  !
  !   logical :: isMergeable
  !
  !   isMergeable = .true.
  !
  !   if (.not.associated(node%children)) then
  !     return
  !   end if
  !   ! if (node%nodeIdx > 0) then !> if the node was a leaf in the last iteration
  !   !   return
  !   ! end if
  !
  !   numElements = 0
  !   do i=1,4
  !     childNode => node%children(i)
  !     call removeUnnecessaryNodes(childNode, tree, newLeafsStack)
  !     if (.not.childNode%isCollapsable .or. associated(childNode%children)) then
  !       isMergeable = .false.
  !     else 
  !       if (childNode%tmpParticleCount < 0) then
  !         numElements = numElements + tree%particleNumbers(childNode%nodeIdx)
  !       else
  !         numElements = numElements + childNode%tmpParticleCount
  !       end if
  !     end if
  !   end do
  !   if (isMergeable .and. numElements <= tree%treeParams%elementMergeThreshold) then
  !     call mergeChildNodes(node, tree)
  !   else
  !     !> if the children can't be merged add them to the stack of new leafs if they haven't been a leaf before
  !     do i = 1,4
  !       childNode => node%children(i)
  !       if (childNode%nodeIdx < 0) then
  !         call push(newLeafsStack, childNode)
  !       end if
  !     end do
  !   end if
  ! end subroutine removeUnnecessaryNodes



  subroutine updateTreeNodes(tree, simParams)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(in) :: simParams

    type(NodeStack) :: stack
    type(QuadTreeNode), pointer :: n

    integer(i4) :: i, idx, num

    type(QuadTreeNodePointer), pointer, dimension(:) :: tmpLeafs
    real(fp), pointer, dimension(:, :) :: tmpParticles
    integer(i1), pointer, dimension(:) :: tmpParticleTypes
    integer(i4), pointer, dimension(:) :: tmpParticleNumbers
    integer(i4), pointer, dimension(:) :: tmpParticleStartIndices


    !> store new leafs in the stack
    call initializeStack(stack)
    !> split all nodes containing more than treeParams%elementSplitThreshold elements
    !> the nodes added to the stack don't need to be considered for merging since otherwise they wouldn't have been created
    !> so they also won't be merged and merging only adds new nodes to the stack
    ! print *, "num leafs ", tree%leafNumber
    call splitNodes(tree, stack)
    ! print *, "num leafs ", tree%leafNumber

    !> remove unnecessary nodes
    !> TODO: run this only on the "old" leaf nodes as newly splitted nodes don't need to be merged
    ! call removeUnnecessaryNodes(tree, stack)
    call removeUnnecessaryNodesDownwards(tree%root, tree, stack)
    ! print *, "num leafs after node removal", tree%leafNumber
    ! print *, stack%topIndex

    call findLeafs(tree%root, stack)
    ! print *, "num leafs after searching leafs", stack%topIndex
    
    tmpLeafs => tree%leafs
    nullify(tree%leafs)
    ! call move_alloc(tree%leafs, tmpLeafs)
    tmpParticleNumbers => tree%particleNumbers
    nullify(tree%particleNumbers)
    ! call move_alloc(tree%particleNumbers, tmpParticleNumbers)
    tmpParticleStartIndices => tree%particleStartIndices
    nullify(tree%particleStartIndices)
    ! call move_alloc(tree%particleStartIndices, tmpParticleStartIndices)
    tmpParticles => tree%particles
    nullify(tree%particles)

    
    ! allocate(tmpLeafs(size(tree%leafs, 1)))
    ! tmpLeafs = tree%leafs
    
    !> store the references to the leafs in the tree%leafs array
    ! deallocate(tree%leafs)
    allocate(tree%leafs(tree%leafNumber))
    !> also update the particle-numbers for each leaf
    ! deallocate(tree%particleNumbers)
    allocate(tree%particleNumbers(tree%leafNumber))
    allocate(tree%particleStartIndices(tree%leafNumber))

    idx = 1
    ! do i = 1,size(tmpLeafs, 1)
    !   !> copy the leaf, if it is still a leaf 
    !   if (associated(tmpLeafs(i)%node)) then
    !     print *, "still leaf", i
    !     tree%leafs(idx)%node => tmpLeafs(i)%node
    !     tree%leafs(idx)%node%nodeIdx = idx
    !     tree%particleNumbers(idx) = tmpParticleNumbers(i)
    !     !> store the start index of the corresponding particles in the node
    !     tree%leafs(idx)%node%tmpParticleIndex = tmpParticleStartIndices(i)
    !     idx = idx + 1
    !   end if
    ! end do
    ! do i = 1, stack%topIndex !> stack%topIndex == size(tree%leafs,1) - idx
    !   call pop(stack,n)
    !   tree%leafs(idx)%node => n
    !   n%nodeIdx = idx
    !   tree%particleNumbers(idx) = n%tmpParticleCount
    !   idx = idx + 1
    ! end do
    do i = 1, stack%topIndex !> stack%topIndex == size(tree%leafs,1) - idx
      call pop(stack,n)
      if (n%tmpParticleCount < 0) then
        tree%particleNumbers(i) = tmpParticleNumbers(n%nodeIdx)
        n%tmpParticleIndex = tmpParticleStartIndices(n%nodeIdx)
      else
        tree%particleNumbers(i) = n%tmpParticleCount
      end if
      tree%leafs(i)%node => n
      n%nodeIdx = i
      ! n%nodeIdx = idx
      ! idx = idx + 1
    end do

    call deleteStack(stack)
    
    !> update the start-indices for the new leafs
    idx = 1
    do i = 1,tree%leafNumber
      tree%particleStartIndices(i) = idx
      idx = idx + (tree%particleNumbers(i) / tree%treeParams%elementChunkSize + 1) * tree%treeParams%elementChunkSize
    end do

    !> remove old pointers to the structures
    deallocate(tree%structures)
    allocate(tree%structures(tree%leafNumber))

    !> allocate the new particles array
    allocate(tree%particles(idx, 5))
    
    !> TODO: can be parallelized
    !$OMP PARALLEL DO private(idx, num, n) shared(tree, tmpParticles,simParams)
    do i = 1, tree%leafNumber
      idx = tree%particleStartIndices(i)
      num = tree%particleNumbers(i)
      n => tree%leafs(i)%node

      !> update the cell stats
      call addParticleCount(n%stats%particleCounter, num)
      n%stats%rho = simParams%m(1) * num / (simParams%V_c * cellWidth(n) * cellHeight(n))

      if (n%tmpParticleCount < 0) then
        tree%particles(idx:idx+num-1, :) = tmpParticles(n%tmpParticleIndex:n%tmpParticleIndex+num-1, :)
        ! tree%particleTypes(idx:idx+num) = tmpParticleTypes(n%tmpParticleIndex:n%tmpParticleIndex+num)
        n%tmpParticleIndex = -1
        n%tmpParticleCount = -1
      else
        ! if (num > 0) then
        if (associated(n%tmpParticles)) then
          tree%particles(idx:idx+num-1, :) = n%tmpParticles
          !> TODO: implement copying the particleTypes in mergeChildNodes
          ! tree%particleTypes(idx:idx+num) = n%tmpParticleTypes
          ! deallocate(n%tmpParticleTypes)
          deallocate(n%tmpParticles)
        end if
        ! deallocate(n%tmpParticles)
        n%tmpParticleCount = -1
      end if
      if (associated(n%structures)) then
        tree%structures(i)%numStructures = size(n%structures,2)
        tree%structures(i)%structures => n%structures
      end if 
    end do
    !$OMP END PARALLEL DO
    deallocate(tmpLeafs)
    deallocate(tmpParticles)
    deallocate(tmpParticleStartIndices)
    deallocate(tmpParticleNumbers)
    ! deallocate(tmpParticleTypes)
  end subroutine updateTreeNodes

  
  subroutine insertParticles(tree, particles, leafIdx)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    real(fp), pointer, dimension(:,:), intent(inout) :: particles
    integer(i4), pointer, dimension(:), intent(inout) :: leafIdx

    integer(i4) :: i, j, start, currentIndex, particleStorage, &
      nodeIdx, particleStartIdx, currentNumParticles
    integer(i4), allocatable, dimension(:) :: numParticles, startIdx
    type(QuadTreeNode), pointer :: node

    if (size(leafIdx) == 0) then
      return
    end if

    allocate(numParticles(size(leafIdx)))
    allocate(startIdx(size(leafIdx)))
    numParticles = 0

    !> sort the particles array based on the leafIdx
    ! call conditionedMergeSort(particles, leafIdx, 1, size(leafIdx))
    call sort(particles, leafIdx)

    !> count the number of occurrences of each leaf index 
    !> and store the first particle index for each leaf 
    currentIndex = leafIdx(1)
    startIdx(1) = 1
    j = 1
    do i = 1, size(leafIdx)
      if (currentIndex /= leafIdx(i)) then
        j = j + 1
        currentIndex = leafIdx(i)
        startIdx(j) = i
      end if
      numParticles(j) = numParticles(j) + 1
    end do
    
    !> ignore the invalid particles (outside the area)
    if (leafIdx(1) == -1) then
      start = 2
    else
      start = 1
    end if 

    !> TODO: can be parallelized
    !> copy the new particles to the particles array
    !!$OMP PARALLEL DO private(nodeIdx, node, currentNumParticles, particleStorage, particleStartIdx) &
    !!$OMP shared(tree, leafIdx, particles, startIdx)
    do i = start,j
      ! nodeIdx = startIdx(leafIdx(i))
      nodeIdx = leafIdx(startIdx(i))
      node => tree%leafs(nodeIdx)%node
      currentNumParticles = tree%particleNumbers(nodeIdx)
      if (nodeIdx == size(tree%leafs)) then
        particleStorage = size(tree%particles, 2)+1 - tree%particleStartIndices(nodeIdx)
      else
        particleStorage = tree%particleStartIndices(nodeIdx+1) - tree%particleStartIndices(nodeIdx)
      end if 
      ! print *, "particle storage 1", particleStorage
      ! particleStorage = (tree%particleNumbers(nodeIdx) / tree%treeParams%elementChunkSize + 1) * tree%treeParams%elementChunkSize
      ! print *, "particle storage 2", particleStorage
      !> if the space for the particles in tree%particles is not sufficient
      !> store the particles temporarily in node%tmpParticles
      if (currentNumParticles + numParticles(i) > particleStorage) then
        ! particleStorage = ((currentNumParticles + numParticles(i)) / tree%treeParams%elementChunkSize + 1) &
        !   * tree%treeParams%elementChunkSize
        particleStorage = currentNumParticles + numParticles(i)
        ! print *, particleStorage
        allocate(node%tmpParticles(particleStorage, 5))
        !> copy the original particles
        node%tmpParticles(:currentNumParticles, :) = &
          tree%particles(tree%particleStartIndices(nodeIdx):tree%particleStartIndices(nodeIdx)+currentNumParticles-1, :)
        !> copy the new particles
        node%tmpParticles(currentNumParticles+1:currentNumParticles+numParticles(i), :) = &
          particles(startIdx(i):startIdx(i)+numParticles(i)-1, :)
        node%tmpParticleCount = currentNumParticles + numParticles(i)
      !> otherwise just copy the new particles
      else
        particleStartIdx = tree%particleStartIndices(nodeIdx)+currentNumParticles
        ! print *, particleStorage, currentNumParticles, numParticles(i), particleStartIdx
        tree%particles(particleStartIdx:particleStartIdx+numParticles(i)-1,:) = &
          particles(startIdx(i):startIdx(i)+numParticles(i)-1,:)
        tree%particleNumbers(nodeIdx) = currentNumParticles + numParticles(i)
      end if
      ! node%numParticles = node%numParticles + numParticles(i)
    end do
    !!$OMP END PARALLEL DO

    deallocate(numParticles)
    deallocate(startIdx)
    deallocate(leafIdx)
    
  end subroutine insertParticles

  subroutine findParticleCells(tree, particles, leafIdx)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    real(fp), pointer, dimension(:,:), intent(in) :: particles
    integer(i4), pointer, dimension(:), intent(inout) :: leafIdx

    integer(i4) :: i, numParticles
    integer(i1) :: j
    type(QuadTreeNode), pointer :: n

    numParticles = size(particles, 1)

    allocate(leafIdx(numParticles))

    !> TODO: can be parallelized
    !$OMP PARALLEL DO private(i, n, j) shared(particles, tree, leafIdx)
    do i = 1, numParticles
      !> mark particles outside the tree as invalid
      if (&
        particles(i, 1) < 0 .or. particles(i, 1) >= tree%treeParams%width .or. &
        particles(i, 2) < 0 .or. particles(i, 2) >= tree%treeParams%height &
      ) then
        leafIdx(i) = -1
        cycle
      end if 
      
      n => tree%root
      !> this is only run while all leafs cells are stored in tree%leafs
      !> and therefore have a corresponding nodeIdx
      do while (n%nodeIdx == -1)
        if (particles(i, 1) < (cellX(n) + cellWidth(n)/2)*tree%treeParams%width) then
          j = 1_i1
        else
          j = 2_i1
        end if
        if (particles(i, 2) >= (cellY(n) + cellHeight(n)/2)*tree%treeParams%height) then
          j = j + 2_i1
        end if
        n => n%children(j)
      end do

      leafIdx(i) = n%nodeIdx
    end do
    !$OMP END PARALLEL DO

  end subroutine findParticleCells

  recursive subroutine deleteSubtree(node, tree)
    implicit none
    type(QuadTreeNode), pointer :: node
    type(QuadTree), pointer :: tree

    type(QuadTreeNode), pointer :: childNode
    integer(i1) :: i
    
    if (.not.associated(node%children)) then
      if (associated(node%structures)) then
        deallocate(node%structures)
      end if
      call deleteCellStats(node%stats)
      return
    end if

    do i=1,4
      childNode => node%children(i)
      call deleteSubtree(childNode, tree)
    end do
    deallocate(node%children)
  end subroutine deleteSubtree

  subroutine deleteTree(tree)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree

    call deleteSubtree(tree%root, tree)
    if (associated(tree%root%structures)) then
      deallocate(tree%root%structures)
    end if
    deallocate(tree%root)

    deallocate(tree%particleNumbers)
    deallocate(tree%particles)
    deallocate(tree%particleTypes)
    deallocate(tree%particleStartIndices)
    deallocate(tree%structures)

    deallocate(tree%leafs)

    deallocate(tree)

  end subroutine deleteTree

  subroutine initializeStack(stack)
    implicit none
    type(NodeStack), intent(inout) :: stack

    nullify(stack%first)
    nullify(stack%top)
    stack%topIndex = 0
  end subroutine initializeStack

  subroutine deleteStack(stack)
    implicit none
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
    implicit none
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
    implicit none
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


end module m_quadtree
