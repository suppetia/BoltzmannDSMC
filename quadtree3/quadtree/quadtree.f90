module m_quadtree
  use m_types, only: fp, i1, i2, i4, i8
  use m_datastructures, only: QuadTreeParameters, SimulationParameters, StructuresPointer, CellStats, &
    addIntegerCount, initializeCellStats, deleteCellStats, initializeStatisticsCell, deleteStatisticsCell, &
    StatisticsCell, addRealCount
  use m_util, only: lineIntersection, sort, k_B, rectContainsRect

  use omp_lib

  implicit none

  type QuadTreeNode

    !> encoding for the position of the cell
    !> 10 | 11
    !> -------
    !> 00 | 01
    !> binary representation of each level
    !> the first level is stored on the most right bits progressing to the left
    !> in the bits 58-63 (starting with 0) the level of the cell is stored 
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

    !> array of CellStats for averaging in fixed rectangles
    type(StatisticsCell), pointer, dimension(:,:) :: statisticsCells => null()

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

    integer(i2) :: i,j
    
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

    !> initialize the statistics cell array
    allocate(tree%statisticsCells(params%numStatisticsCellRows,params%numStatisticsCellColumns))
    do i = 1_i2, params%numStatisticsCellRows
      do j = 1_i2, params%numStatisticsCellColumns
        call initializeStatisticsCell(tree%statisticsCells(i,j), params)
      end do 
    end do 

  end subroutine initializeTree 


  recursive subroutine splitNode(node, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree

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

    integer(i4) :: i
    integer(i1) :: j
    integer(i2) :: k
    integer(i4), dimension(4) :: idx

    type(QuadTreeNode), pointer :: childNode

    allocate(node%children(4))

    !> update the number of leafs in the tree
    !$OMP CRITICAL
    tree%leafNumber = tree%leafNumber + 3
    !$OMP END CRITICAL

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
      !> store the information about the parent cell
      node%children(i)%parent => node

      !> split the cellstats evenly for all child nodes
      call initializeCellStats(node%children(i)%stats, tree%treeParams)
      !> split the particle stats
      node%children(i)%stats%particleCounter%history = node%stats%particleCounter%history / 4
      node%children(i)%stats%particleCounter%average = node%stats%particleCounter%average / 4
      node%children(i)%stats%particleCounter%currentPosition = node%stats%particleCounter%currentPosition

      !> split the cell stats
      node%children(i)%stats%statsCounter%history = node%stats%statsCounter%history / 4
      node%children(i)%stats%statsCounter%average = node%stats%statsCounter%average / 4
      node%children(i)%stats%statsCounter%currentPosition = node%stats%statsCounter%currentPosition

      !> split the species stats
      do k = 1, tree%treeParams%numParticleSpecies
        node%children(i)%stats%speciesStatsCounter(k)%counter%history = node%stats%speciesStatsCounter(k)%counter%history / 4
        node%children(i)%stats%speciesStatsCounter(k)%counter%average = node%stats%speciesStatsCounter(k)%counter%average / 4
        node%children(i)%stats%speciesStatsCounter(k)%counter%currentPosition = &
          node%stats%speciesStatsCounter(k)%counter%currentPosition
      end do 
    end do
    call deleteCellStats(node%stats)

    !> split structures if present
    if (associated(node%structures)) then
      ! idx(:) = 1
      ! allocate(tmpStructures(4, 6, size(node%structures,2)))
      ! do i = 1, size(node%structures,2)
      !   
      !   !> store the orientation of the structure line to find the corresponding subcells later
      !   if (node%structures(1, i) - node%structures(3, i) < 0._fp) then
      !     structureDirectionX = 1_i1
      !   else
      !     structureDirectionX = -1_i1
      !   end if
      !   if (node%structures(2, i) - node%structures(4, i) < 0._fp) then
      !     structureDirectionY = 1_i1
      !   else
      !     structureDirectionY = -1_i1
      !   end if
      !   
      !   verticalIntersection = lineIntersection([x+newWidth, y, x+newWidth, y+newHeight*2], node%structures(:4,i))
      !   horizontalIntersection = lineIntersection([x, y+newHeight, x+newWidth*2, y+newHeight], node%structures(:4,i))
      !
      !   !> determine the subcells which the line crosses
      !   if (node%structures(1,i) >= x + newWidth) then
      !     j = 2_i1
      !   else
      !     j = 1_i1
      !   end if
      !   if (node%structures(2,i) >= y + newHeight) then
      !     j = j + 2_i1
      !   end if
      !   tmpStructures(j, 1:2, idx(j)) = node%structures(1:2,i)
      !   
      !   if (any(verticalIntersection /= -1._fp)) then
      !     if (any(horizontalIntersection /= -1._fp)) then
      !       !> TODO: handle the case where the structure is exactly aligning with the middle
      !       if ((verticalIntersection(1)-horizontalIntersection(1))*structureDirectionX < 0) then
      !         !> if the vertical intersection happens first
      !         tmpStructures(j,3:4,idx(j)) = verticalIntersection
      !         tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !         idx(j) = idx(j) + 1
      !         j = j + structureDirectionX
      !         tmpStructures(j,1:2,idx(j)) = verticalIntersection
      !   
      !         tmpStructures(j,3:4,idx(j)) = horizontalIntersection
      !         tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !         idx(j) = idx(j) + 1
      !         j = j + 2_i1 * structureDirectionY
      !         tmpStructures(j,1:2,idx(j)) = horizontalIntersection
      !   
      !       else
      !         !> if the horizontalIntersection happens first
      !         tmpStructures(j,3:4,idx(j)) = horizontalIntersection
      !         tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !         idx(j) = idx(j) + 1
      !         j = j + 2_i1 * structureDirectionY
      !         tmpStructures(j,1:2,idx(j)) = horizontalIntersection
      !      
      !         tmpStructures(j,3:4,idx(j)) = verticalIntersection
      !         tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !         idx(j) = idx(j) + 1
      !         j = j + structureDirectionX
      !         tmpStructures(j,1:2,idx(j)) = verticalIntersection
      !       end if
      !     else
      !       !> only the verticalIntersection happens
      !       tmpStructures(j,3:4,idx(j)) = verticalIntersection
      !       tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !       idx(j) = idx(j) + 1
      !       j = j + structureDirectionX
      !       tmpStructures(j,1:2,idx(j)) = verticalIntersection
      !     end if
      !   elseif (any(horizontalIntersection /= -1._fp)) then
      !     !> only the horizontalIntersection happens
      !     tmpStructures(j,3:4,idx(j)) = horizontalIntersection
      !     tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !     idx(j) = idx(j) + 1
      !     j = j + 2_i1 * structureDirectionY
      !     tmpStructures(j,1:2,idx(j)) = horizontalIntersection
      !   end if
      !   tmpStructures(j,3:4,idx(j)) = node%structures(3:4,i)
      !   tmpStructures(j,5:6,idx(j)) = node%structures(5:6,i)
      !   idx(j) = idx(j) + 1
      !   
      ! end do
      !   
      ! do j = 1_i1,4_i1
      !   if (idx(j)-1 > 0) then
      !     allocate(node%children(j)%structures(6, idx(j)-1))
      !     node%children(j)%structures(:, :idx(j)-1) = tmpStructures(j, :, :idx(j)-1)
      !   end if
      ! end do

      do j = 1_i1,4_i1
        allocate(node%children(j)%structures(6, size(node%structures,2)))
        node%children(j)%structures = node%structures
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
    allocate(tmpParticleTypes(numParticles))

    ! !> copy the particles to a temporary array for fast access and easy sorting
    ! tmpParticles = tree%particles(idxStartParticles:idxStartParticles+numParticles-1, :)
    ! tmpParticleTypes = tree%particleTypes(idxStartParticles:idxStartParticles+numParticles-1)
    !> copy the particles to a temporary array for fast access and easy sorting
    !> TODO: check whether this is faster on the original arrays
    ! print *, node%tmpParticleCount, associated(node%tmpParticles)
    if (node%tmpParticleCount < 0) then
      tmpParticles = tree%particles(idxStartParticles:idxStartParticles+numParticles-1, :)
      tmpParticleTypes = tree%particleTypes(idxStartParticles:idxStartParticles+numParticles-1)
    else
      tmpParticles = node%tmpParticles
      tmpParticleTypes = node%tmpParticleTypes
      deallocate(node%tmpParticles)
      deallocate(node%tmpParticleTypes)
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
      allocate(node%children(j)%tmpParticleTypes(numParticlesPerChildCell(j)))
      node%children(j)%tmpParticleCount = numParticlesPerChildCell(j)
    end do
      

    !> actually sort the particles into the tree array
    !> TODO: can be parallelized
    do i = 1, numParticles
      ! tree%particles(:, idx(particleLocations(i))) = tmpParticles(:, i)
      ! tree%particleTypes(idx(particleLocations(i))) = tmpParticleTypes(i)

      node%children(particleLocations(i))%tmpParticles(idx(particleLocations(i)),:) = tmpParticles(i, :)
      node%children(particleLocations(i))%tmpParticleTypes(idx(particleLocations(i))) = tmpParticleTypes(i)

      idx(particleLocations(i)) = idx(particleLocations(i)) + 1
    end do
    deallocate(tmpParticles)
    deallocate(tmpParticleTypes)
    deallocate(particleLocations)

    !> check if all the child cells fulfill the threshold condition
    if (level < 63) then
      do j = 1,4
        childNode => node%children(j)
        if (childNode%tmpParticleCount > tree%treeParams%elementSplitThreshold) then
          call splitNode(childNode, tree)
        end if
      end do
    end if 
    
  end subroutine splitNode

  !> split all nodes that exceed the threshold defined in the tree parameters
  subroutine splitNodes(tree)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree

    ! type(QuadTreeNode), pointer :: node, childNode
    integer(i4) :: i

    !> TODO: can be parallelized if the leaf number is updated in a thread safe way
    !$OMP PARALLEL DO private(i) shared(tree)
    do i = 1, size(tree%leafs)
      !> use node%tmpParticleCount as this metric is updated when inserting particles
      if (max(tree%particleNumbers(i),tree%leafs(i)%node%tmpParticleCount) &
        >= tree%treeParams%elementSplitThreshold) then
        call splitNode(tree%leafs(i)%node, tree)
        nullify(tree%leafs(i)%node)
      end if
    end do
    !$OMP END PARALLEL DO
  end subroutine splitNodes

  !> merge the child nodes making node a leaf
  subroutine mergeChildNodes(node, tree)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(QuadTreeNode), pointer :: childNode

    integer(i1) :: i
    integer(i2) :: k
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
    allocate(node%tmpParticleTypes(numParticles))

    !> merge the cell stats from the child nodes
    call initializeCellStats(node%stats, tree%treeParams)
    node%stats%particleCounter%currentPosition = node%children(1)%stats%particleCounter%currentPosition
    node%stats%statsCounter%currentPosition = node%children(1)%stats%particleCounter%currentPosition
    do j = 1, tree%treeParams%numParticleSpecies
      node%stats%speciesStatsCounter(j)%counter%currentPosition = node%children(1)%stats%particleCounter%currentPosition
    end do 

    idx = 1

    do i=1,4
      childNode => node%children(i)
      !> update particle stats
      node%stats%particleCounter%history = node%stats%particleCounter%history + childNode%stats%particleCounter%history
      node%stats%particleCounter%average = node%stats%particleCounter%average + childNode%stats%particleCounter%average
      !> update cell stats
      node%stats%statsCounter%history = node%stats%statsCounter%history + childNode%stats%statsCounter%history
      node%stats%statsCounter%average = node%stats%statsCounter%average + childNode%stats%statsCounter%average
      !> update species stats
      do k = 1, tree%treeParams%numParticleSpecies
        node%stats%speciesStatsCounter(k)%counter%history = node%stats%speciesStatsCounter(k)%counter%history &
          + childNode%stats%speciesStatsCounter(k)%counter%history
        node%stats%speciesStatsCounter(k)%counter%average = node%stats%speciesStatsCounter(k)%counter%average &
          + childNode%stats%speciesStatsCounter(k)%counter%average
      end do

      if (childNode%tmpParticleCount < 0) then !> if the node was a leaf in the last iteration
        numParticles = tree%particleNumbers(childNode%nodeIdx)
        j = tree%particleStartIndices(childNode%nodeIdx)
        node%tmpParticles(idx:idx+numParticles-1, :) = tree%particles(j:j+numParticles-1, :)
        node%tmpParticleTypes(idx:idx+numParticles-1) = tree%particleTypes(j:j+numParticles-1)
      else
        numParticles = childNode%tmpParticleCount
        node%tmpParticles(idx:idx+numParticles-1, :) = childNode%tmpParticles
        node%tmpParticleTypes(idx:idx+numParticles-1) = childNode%tmpParticleTypes
        deallocate(childNode%tmpParticles)
        deallocate(childNode%tmpParticleTypes)
        ! if (associated(childNode%tmpParticles)) then
        !   if (numParticles > 0) then
        !     node%tmpParticles(:, idx:idx+numParticles-1) = childNode%tmpParticles
        !   end if 
        !   deallocate(childNode%tmpParticles)
        ! end if
      end if
      if (childNode%nodeIdx > 0) then
        !> remove the pointer to this childNode from the tree
        tree%leafs(childNode%nodeIdx)%node => null()
      end if 
      idx = idx + numParticles
      ! childNode%parent => null()
      call deleteCellStats(childNode%stats)
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
        if (.not.associated(childNode%children)) then
          call push(newLeafsStack, childNode)
        end if
      end do
    end if
  end subroutine removeUnnecessaryNodesDownwards

  recursive subroutine findLeafs(n, leafsStack)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: n
    type(NodeStack), intent(inout) :: leafsStack

    type(QuadTreeNode), pointer :: childNode
    integer(i1) :: i

    if (associated(n%children)) then
      do i = 1,4
        childNode => n%children(i)
        call findLeafs(childNode, leafsStack)
      end do
    else
      call push(leafsStack, n)
    end if
  end subroutine findLeafs

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

  subroutine updateTreeNodes(tree, simParams)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(in) :: simParams

    type(NodeStack) :: stack
    type(QuadTreeNode), pointer :: n

    integer(i4) :: i, idx, idx2, num
    integer(i2) :: j

    type(QuadTreeNodePointer), pointer, dimension(:) :: tmpLeafs
    real(fp), pointer, dimension(:, :) :: tmpParticles
    integer(i1), pointer, dimension(:) :: tmpParticleTypes
    integer(i4), pointer, dimension(:) :: tmpParticleNumbers
    integer(i4), pointer, dimension(:) :: tmpParticleStartIndices

    real(fp) :: density !> particle density
    real(fp) :: mass !> mass of particles in a cell
    real(fp) :: rho !> mass density
    real(fp) :: c0_x, c0_y, c0_z !> mean velocity components
    real(fp) :: cx_sq, cy_sq, cz_sq !> mean squared relative velocity components
    real(fp) :: c_sq !> mean squared relative velocity (relative = c - c_mean)
    real(fp) :: p !> pressure
    real(fp) :: T !> translational temperature

    !> store new leafs in the stack
    call initializeStack(stack)
    !> split all nodes containing more than treeParams%elementSplitThreshold elements
    !> the nodes added to the stack don't need to be considered for merging since otherwise they wouldn't have been created
    !> so they also won't be merged and merging only adds new nodes to the stack
    ! print *, "num leafs ", tree%leafNumber
    call splitNodes(tree)
    print *, "num leafs ", tree%leafNumber

    !> remove unnecessary nodes
    !> TODO: run this only on the "old" leaf nodes as newly splitted nodes don't need to be merged
    ! call removeUnnecessaryNodes(tree, stack)
    call removeUnnecessaryNodesDownwards(tree%root, tree, stack)
    ! do i = 1, size(tree%leafs)
    !   if (associated(tree%leafs(i)%node)) then
    !     if (associated(tree%leafs(i)%node%children)) then
    !       call removeUnnecessaryUpwardNodes(tree%leafs(i)%node%parent, tree, stack)
    !     end if
    !   end if 
    ! end do 
    ! print *, "num leafs after node removal", tree%leafNumber, stack%topIndex
    ! print *, stack%topIndex

    ! call findLeafs(tree%root, stack)
    ! print *, "num leafs after searching leafs", stack%topIndex
    ! print *, stack%topIndex - tree%leafNumber
    ! print *, stack%topIndex - size(tree%leafs)
    
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

    tmpParticleTypes => tree%particleTypes
    nullify(tree%particleTypes)

    
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
    allocate(tree%particleTypes(idx))
    
    !> TODO: can be parallelized
    !$OMP PARALLEL DO &
    !$OMP private(idx,idx2,num,n,density,mass,rho,c0_x,c0_y,c0_z,cx_sq,cy_sq,cz_sq,c_sq,p,T,j) &
    !$OMP shared(tree,tmpParticles,simParams)
    do i = 1, tree%leafNumber
      idx = tree%particleStartIndices(i)
      if (i<tree%leafNumber) then
        idx2 = tree%particleStartIndices(i+1)
      else
        idx2 = size(tree%particleTypes)+1
      end if 
      num = tree%particleNumbers(i)
      n => tree%leafs(i)%node


      if (n%tmpParticleCount < 0) then
        tree%particles(idx:idx+num-1, :) = tmpParticles(n%tmpParticleIndex:n%tmpParticleIndex+num-1, :)
        tree%particles(idx+num:idx2-1,:) = -1
        tree%particleTypes(idx:idx+num-1) = tmpParticleTypes(n%tmpParticleIndex:n%tmpParticleIndex+num-1)
        tree%particleTypes(idx+num:idx2-1) = -1
        n%tmpParticleIndex = -1
        n%tmpParticleCount = -1
      else
        ! if (num > 0) then
        if (associated(n%tmpParticles)) then
          tree%particles(idx:idx+num-1, :) = n%tmpParticles
          tree%particleTypes(idx:idx+num-1) = n%tmpParticleTypes
          tree%particles(idx+num:idx2-1,:) = -1
          tree%particleTypes(idx+num:idx2-1) = -1
          deallocate(n%tmpParticleTypes)
          deallocate(n%tmpParticles)
          n%tmpParticleCount = -1
        end if
        ! deallocate(n%tmpParticles)
        ! n%tmpParticleCount = -1
      end if

      ! print *, i, num, tree%particles(idx:idx2-1, 1)
      !> update the cell stats
      ! call addParticleCount(n%stats%particleCounter, num)
      call addIntegerCount(n%stats%particleCounter, [num])

      ! !> TODO: STORE THE STATS FOR EACH PARTICLE TYPE ?
      ! !> calculate the node stats
      ! density = num / (simParams%V_c * cellWidth(n) * cellHeight(n)) * simParams%F_N
      ! !> TODO: update using the particle types
      ! rho = simParams%m(1) * density
      !
      ! !> calculate the mean velocity components
      ! c0_x = sum(tree%particles(idx:idx+num-1, 3)) / num
      ! c0_y = sum(tree%particles(idx:idx+num-1, 4)) / num
      ! c0_z = sum(tree%particles(idx:idx+num-1, 5)) / num
      ! !> calculate the mean relative velocities
      ! cx_sq = sum((tree%particles(idx:idx+num-1, 3)-c0_x) ** 2)/num
      ! cy_sq = sum((tree%particles(idx:idx+num-1, 4)-c0_y) ** 2)/num
      ! cz_sq = sum((tree%particles(idx:idx+num-1, 5)-c0_z) ** 2)/num
      !
      ! c_sq = cx_sq + cy_sq + cz_sq
      ! p = rho / 3 * c_sq
      ! !> translational temperature 3/2 k * T_tr = 1/2 m * \overbar{c^2}
      ! T = simParams%m(1) / k_B * c_sq / 3 
      !
      ! n%stats%n = num / (simParams%V_c * cellWidth(n) * cellHeight(n)) * simParams%F_N
      ! !> TODO: update using the particle types
      ! n%stats%rho = simParams%m(1) * n%stats%n
      !
      ! !> calculate the mean velocity components
      ! n%stats%c0_x = sum(tree%particles(idx:idx+num-1, 3)) / num
      ! n%stats%c0_y = sum(tree%particles(idx:idx+num-1, 4)) / num
      ! n%stats%c0_z = sum(tree%particles(idx:idx+num-1, 5)) / num
      ! !> calculate the mean relative velocities
      ! n%stats%cx_sq = sum((tree%particles(idx:idx+num-1, 3)-n%stats%c0_x) ** 2)/num
      ! n%stats%cy_sq = sum((tree%particles(idx:idx+num-1, 4)-n%stats%c0_y) ** 2)/num
      ! n%stats%cz_sq = sum((tree%particles(idx:idx+num-1, 5)-n%stats%c0_z) ** 2)/num
      !
      ! n%stats%c_sq = n%stats%cx_sq + n%stats%cy_sq + n%stats%cz_sq
      ! n%stats%p = n%stats%rho / 3 * n%stats%c_sq
      ! !> translational temperature 3/2 k * T_tr = 1/2 m * \overbar{c^2}
      ! n%stats%T = simParams%m(1) / k_B * n%stats%c_sq / 3 
      !
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
    deallocate(tmpParticleTypes)
  end subroutine updateTreeNodes

  
  subroutine insertParticles(tree, particles, particleTypes, leafIdx)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    real(fp), pointer, dimension(:,:), intent(inout) :: particles
    integer(i1), pointer, dimension(:), intent(inout) :: particleTypes
    integer(i4), pointer, dimension(:), intent(inout) :: leafIdx

    integer(i4) :: i, j, start, currentIndex, particleStorage, &
      nodeIdx, particleStartIdx, currentNumParticles
    integer(i4), allocatable, dimension(:) :: numParticles, startIdx
    type(QuadTreeNode), pointer :: node

    if (size(leafIdx) == 0) then
      deallocate(leafIdx)
      return
    end if

    allocate(numParticles(size(leafIdx)))
    allocate(startIdx(size(leafIdx)))
    numParticles = 0

    !> sort the particles array based on the leafIdx
    ! call conditionedMergeSort(particles, leafIdx, 1, size(leafIdx))
    call sort(particles, particleTypes, leafIdx)

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
        particleStorage = size(tree%particles, 1)+1 - tree%particleStartIndices(nodeIdx)
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
        allocate(node%tmpParticleTypes(particleStorage))
        !> copy the original particles
        node%tmpParticles(:currentNumParticles, :) = &
          tree%particles(tree%particleStartIndices(nodeIdx):tree%particleStartIndices(nodeIdx)+currentNumParticles-1, :)
        node%tmpParticleTypes(:currentNumParticles) = &
          tree%particleTypes(tree%particleStartIndices(nodeIdx):tree%particleStartIndices(nodeIdx)+currentNumParticles-1)
        !> copy the new particles
        node%tmpParticles(currentNumParticles+1:currentNumParticles+numParticles(i), :) = &
          particles(startIdx(i):startIdx(i)+numParticles(i)-1, :)
        node%tmpParticleTypes(currentNumParticles+1:currentNumParticles+numParticles(i)) = &
          particleTypes(startIdx(i):startIdx(i)+numParticles(i)-1)
        node%tmpParticleCount = currentNumParticles + numParticles(i)
      !> otherwise just copy the new particles
      else
        particleStartIdx = tree%particleStartIndices(nodeIdx)+currentNumParticles
        ! print *, particleStorage, currentNumParticles, numParticles(i), particleStartIdx
        tree%particles(particleStartIdx:particleStartIdx+numParticles(i)-1,:) = &
          particles(startIdx(i):startIdx(i)+numParticles(i)-1,:)
        tree%particleTypes(particleStartIdx:particleStartIdx+numParticles(i)-1) = &
          particleTypes(startIdx(i):startIdx(i)+numParticles(i)-1)
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


  !> find the node that fully contains a rectangle
  subroutine getNodeThatFullyContainsRect(rect, tree, node)
    implicit none
    real(fp), dimension(4), intent(in) :: rect
    type(QuadTree), pointer, intent(in) :: tree
    type(QuadTreeNode), pointer, intent(inout) :: node

    logical :: fitsInChild
    integer(i4) :: i
    real(fp), dimension(4) :: childNodeRect

    node => tree%root
    do while (associated(node%children))
      fitsInChild = .false.
      do i = 1, 4
        childNodeRect = (/cellX(node%children(i)) * tree%treeParams%width, cellY(node%children(i)) * tree%treeParams%height, &
          cellWidth(node%children(i)) * tree%treeParams%width, cellHeight(node%children(i)) * tree%treeParams%height/)
        if (rectContainsRect(childNodeRect, rect)) then
          node => node%children(i)
          fitsInChild = .true.
          exit
        end if 
      end do 
      if (.not. fitsInChild) then
        exit
      end if 
    end do
  end subroutine getNodeThatFullyContainsRect

  subroutine findCellsContainedInRect(rect, tree, stack)
    implicit none
    real(fp), dimension(4), intent(in) :: rect
    type(QuadTree), pointer, intent(in) :: tree
    type(NodeStack), intent(inout) :: stack

    type(QuadTreeNode), pointer :: nodeFullyContainingRect

    call getNodeThatFullyContainsRect(rect, tree, nodeFullyContainingRect)
    call initializeStack(stack)
    call findLeafs(nodeFullyContainingRect, stack)
    
  end subroutine findCellsContainedInRect


  ! subroutine calculateAveragesInStatisticsCells(tree, simParams)
  !   implicit none
  !   type(QuadTree), pointer, intent(inout) :: tree
  !   type(SimulationParameters), pointer, intent(in) :: simParams
  !
  !   type(StatisticsCell), pointer :: stats
  !   integer(i2) :: i,j
  !   logical, dimension(:), pointer :: maskX, maskY
  !   real(fp) :: x,y,width,height
  !
  ! end subroutine calculateAveragesInStatisticsCells 

  subroutine calculateAveragesInStatisticsCells(tree, simParams)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(in) :: simParams

    type(QuadTreeNode), pointer :: n
    type(StatisticsCell), pointer :: stats
    type(NodeStack) :: stack
    integer(i1) :: m
    integer(i2) :: i,j
    integer(i4) :: k,l,idx, num
    real(fp) :: x,y,width,height
    real(fp), dimension(:), allocatable :: rho,T,p,density
    real(fp), dimension(:), allocatable :: cx,cy,cz,cx_sq,cy_sq,cz_sq
    integer(i4), dimension(:), allocatable :: numInRect

    integer(i2), dimension(:,:), allocatable :: statCellID
    logical, pointer, dimension(:) :: maskLoc, mask
    logical, pointer, dimension(:,:) :: maskType

    width = tree%treeParams%width / tree%treeParams%numStatisticsCellColumns
    height = tree%treeParams%height / tree%treeParams%numStatisticsCellRows

    allocate(statCellID(size(tree%particleTypes),2))
    allocate(maskType(size(tree%particleTypes), simParams%numParticleSpecies))
    statCellID = 0


    !$OMP PARALLEL shared(tree, maskType, statCellID,width,height) &
    !$OMP private(x,y,stats,i,j,l,m,numInRect,mask,maskLoc,cx,cy,cz,cx_sq,cy_sq,cz_sq,rho,T,p,density)

    !$OMP sections
    !$OMP section
    ! where (tree%particles(:,1) >= 0.0_fp)
    !   statCellID(:, 1) = int(tree%particles(:, 1) / width) + 1
    !   statCellID(:, 2) = int(tree%particles(:, 2) / height) + 1
    ! end where
    statCellID(:, 1) = int(floor(tree%particles(:, 1) / width), i2) + 1_i2
    !$OMP section
    statCellID(:, 2) = int(floor(tree%particles(:, 2) / height), i2) + 1_i2
    ! !$OMP section
    ! do i = 1_i2, tree%treeParams%numStatisticsCellColumns
    !   x = (i-1)*width
    !   where (tree%particles(:, 1) >= x .and. tree%particles(:, 1) < x+width)
    !     statCellID(:, 1) = i
    !   end where
    ! end do 
    ! !$OMP section
    ! do i = 1_i2, tree%treeParams%numStatisticsCellRows
    !   y = (i-1)*height
    !   where (tree%particles(:, 2) >= y .and. tree%particles(:, 2) < y+height)
    !     statCellID(:, 2) = i
    !   end where
    ! end do 
    !$OMP section
    do m = 1_i1, simParams%numParticleSpecies
      maskType(:, m) = tree%particleTypes == m
    end do 
    !$OMP end sections

    j = simParams%numParticleSpecies+1
    allocate(numInRect(j))
    allocate(cx(j),cy(j),cz(j))
    allocate(cx_sq(j),cy_sq(j),cz_sq(j))
    allocate(rho(j),T(j),p(j),density(j))
    allocate(maskLoc(size(statCellID, 1)))
    allocate(mask(size(statCellID, 1)))
    !$OMP DO
    do i = 1_i2, tree%treeParams%numStatisticsCellColumns
      do j = 1_i2, tree%treeParams%numStatisticsCellRows
        stats => tree%statisticsCells(j,i)

        numInRect = 0
        maskLoc = statCellID(:, 1) == i .and. statCellID(:, 2) == j

        do l = 1, tree%treeParams%numParticleSpecies
          mask = maskLoc .and. maskType(:, l)
          numInRect(l+1) = count(mask)
          cx(l+1) = sum(tree%particles(:,3), mask)
          cy(l+1) = sum(tree%particles(:,4), mask)
          cz(l+1) = sum(tree%particles(:,5), mask)
          cx_sq(l+1) = sum(tree%particles(:,3)**2, mask)
          cy_sq(l+1) = sum(tree%particles(:,4)**2, mask)
          cz_sq(l+1) = sum(tree%particles(:,5)**2, mask)
        end do 

        numInRect(1) = sum(numInRect(2:))
        if (numInRect(1) > 0) then
          !> mean velocity components
          !> cX_0 = (sum_{species p} m_p * mean velocity of species p)/(sum_{species p} numParticles of species p)
          cx(1) = sum(simParams%m*cx(2:))/sum(simParams%m*numInRect(2:))
          cy(1) = sum(simParams%m*cy(2:))/sum(simParams%m*numInRect(2:))
          cz(1) = sum(simParams%m*cz(2:))/sum(simParams%m*numInRect(2:))
        end if 

        !> number density
        density = numInRect / (simParams%V_c/tree%treeParams%width/tree%treeParams%height * width * height) * simParams%F_N
        !> mass density
        rho(2:) = density(2:)*simParams%m
        rho(1) = sum(rho(2:))

        !> translational temperature for each species (see eq.4.40 Bird2012)
        where (numInRect(2:) > 0)
          T(2:) = simParams%m * ((cx_sq(2:)+cy_sq(2:)+cz_sq(2:)) / numInRect(2:) - cx(1)**2-cy(1)**2-cz(1)**2)/(3*k_B)
        elsewhere
          T(2:) = 0.0
        end where
        !> total translational temperature (see eq.4.39 Bird2012)
        if (numInRect(1) > 0) then
          T(1) = (sum(simParams%m * (cx_sq(2:)+cy_sq(2:)+cz_sq(2:))) &
            - sum(simParams%m * numInRect(2:))*(cx(1)**2+cy(1)**2+cz(1)**2)&
            )/(3*k_B*numInRect(1))
        else
          T(1) = 0.0
        end if

        !> pressure
        p = density * k_B * T
        !> calculate the averages velocities for each particle species
        where (numInRect(2:) /= 0)
          cx(2:) = cx(2:) / numInRect(2:)
          cy(2:) = cy(2:) / numInRect(2:)
          cz(2:) = cz(2:) / numInRect(2:)
          cx_sq(2:) = cx_sq(2:) / numInRect(2:)
          cy_sq(2:) = cy_sq(2:) / numInRect(2:)
          cz_sq(2:) = cz_sq(2:) / numInRect(2:)
        end where

        call deleteStack(stack)

        call addIntegerCount(stats%numParticles, numInRect)
        call addRealCount(stats%n, density)
        call addRealCount(stats%rho, rho)
        call addRealCount(stats%cx_0, cx)
        call addRealCount(stats%cy_0, cy)
        call addRealCount(stats%cz_0, cz)
        call addRealCount(stats%p, p)
        call addRealCount(stats%T, T)

      end do 
    end do 
    !$OMP END DO
    deallocate(numInRect)
    deallocate(cx,cy,cz)
    deallocate(cx_sq,cy_sq,cz_sq)
    deallocate(rho,T,p,density)
    deallocate(maskLoc, mask)
    !$OMP END PARALLEL

    deallocate(maskType, statCellID)


  end subroutine calculateAveragesInStatisticsCells 



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

    integer(i2) :: i,j

    do j = 1_i2, tree%treeParams%numStatisticsCellColumns
      do i = 1_i2, tree%treeParams%numStatisticsCellRows
        call deleteStatisticsCell(tree%statisticsCells(i,j))
      end do 
    end do 
    deallocate(tree%statisticsCells)

    deallocate(tree%treeParams)
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
