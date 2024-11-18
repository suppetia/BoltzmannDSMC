module m_quadtree_io
  use m_types, only: fp, i1, i4, i8
  use m_quadtree, only: QuadTree, QuadTreeNode, cellX, cellY, cellWidth, cellHeight, &
    NodeStack, initializeStack, deleteStack, push, pop
  use m_matrix_io, only: writeRealMatrixToH5Dataset
  use m_datastructures, only: QuadTreeParameters, initializeCellStats
  implicit none
contains

  subroutine writeQuadTreeToHDF5(tree, filename, datasetID, writeParticles)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: datasetID
    logical, intent(in) :: writeParticles

    real(fp), pointer, dimension(:,:) :: matrix
    integer(i4), pointer, dimension(:) :: particleStartPosition
    type(QuadTreeNode), pointer :: n
    integer(i4) :: i, numParticles, idx1, idx2, num
    real(fp) :: width, height
    character(len=5) :: datasetName

    integer(i4) :: error

    width = tree%treeParams%width
    height = tree%treeParams%height
    allocate(matrix(5, tree%leafNumber))
    do i = 1, tree%leafNumber
      n => tree%leafs(i)%node
      matrix(:, i) = [cellX(n)*width, cellY(n)*height, cellWidth(n)*width, cellHeight(n)*height, n%stats%rho]
    end do 


    write(datasetName, "(i5.5)") datasetID

    call writeRealMatrixToH5Dataset(filename, datasetName//"_tree", matrix,tree%leafNumber, 5, error)
    deallocate(matrix)

    if (writeParticles) then
      !> write the particle matrix to the file
      numParticles = 0
      allocate(particleStartPosition(tree%leafNumber))
      particleStartPosition = 1
      do i = 1, tree%leafNumber
        num = tree%particleNumbers(i)
        numParticles = numParticles + num
        particleStartPosition(i+1:) = particleStartPosition(i+1:) + num
      end do 
      allocate(matrix(5,numParticles))
      !> TODO: can be parallelized
      do i = 1, tree%leafNumber
        num = tree%particleNumbers(i)
        idx1 = tree%particleStartIndices(i)
        idx2 = particleStartPosition(i)
        matrix(:, idx2:idx2+num-1) = tree%particles(:, idx1:idx1+num-1)
      end do
      
      call writeRealMatrixToH5Dataset(filename, datasetName//"_particles", matrix, numParticles, 5, error)
      deallocate(matrix)
      deallocate(particleStartPosition)
    end if 
  end subroutine writeQuadTreeToHDF5 


  subroutine createTreeFromFile(tree, params, filename)
    implicit none
    type(QuadTree), pointer, intent(out) :: tree
    type(QuadTreeParameters), pointer, intent(inout) :: params
    character(len=*), intent(in) :: filename

    type(NodeStack) :: stack
    type(QuadTreeNode), pointer :: n, parent
    integer(i4) :: io, ioStatus, level, maxLevel
    integer(i4) :: i,j,k, nRows, nLeafs

    integer(i8) :: cellID
    integer(i1) :: numberOfContours, isLeaf
    real(fp) :: width, height

    io = 10
    !> open the file for reading
    open(unit=io, file=filename, status='old', action="read", iostat=ioStatus)

    read(io, *, iostat=ioStatus) nRows, width, height

    params%width = width
    params%height = height

    !> initialize the first part of the tree
    allocate(tree)
    allocate(tree%root)
    tree%treeParams => params

    call initializeStack(stack)

    do k = 1, nRows
      read(io, *) cellID, isLeaf, numberOfContours

      !> if the file end is reached
      if (ioStatus /= 0) then
        exit
      end if 

      maxLevel = shiftr(cellID, 58) - 1
      if (maxLevel == -1) then !> if the line describes the root
        n => tree%root
      else
        parent => tree%root
        do level = 1, maxLevel
          parent => parent%children(ibits(cellID, 2*(level-1), 2) + 1)
        end do
        n => parent%children(ibits(cellID, 2*(maxLevel), 2) + 1)
        n%parent => parent
      end if 
      n%cellID = cellID
      n%isCollapsable = .false.
      if (isLeaf == 0) then
        allocate(n%children(4))
      else
        call push(stack, n)
        if (numberOfContours > 0) then
          allocate(n%structures(6, numberOfContours))
          read(io, *) ((n%structures(i,j), i=1,6), j=1, numberOfContours)
        end if 
        call initializeCellStats(n%stats, tree%treeParams)
      end if 
    end do  
    close(io)

    !> allocate the tree array belonging to the leafs
    nLeafs = stack%topIndex
    tree%leafNumber = nLeafs
    allocate(tree%leafs(nLeafs))
    allocate(tree%structures(nLeafs))
    allocate(tree%particleNumbers(nLeafs))
    tree%particleNumbers = 0
    allocate(tree%particleStartIndices(nLeafs))
    allocate(tree%particles(5, params%elementChunkSize * nLeafs))
    allocate(tree%particleTypes(params%elementChunkSize * nLeafs))
    do k = 1,nLeafs 
      call pop(stack, n)
      tree%leafs(k)%node => n
      n%nodeIdx = k
      if (associated(n%structures)) then
        tree%structures(k)%structures => n%structures
        tree%structures(k)%numStructures = size(n%structures,2)
      end if 
      tree%particleStartIndices(k) = 1 + (k-1)*params%elementChunkSize
    end do 
    call deleteStack(stack)

    
  end subroutine createTreeFromFile 


end module m_quadtree_io 
