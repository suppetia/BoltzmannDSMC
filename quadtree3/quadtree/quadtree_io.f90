module m_quadtree_io
  use m_types, only: fp, i1, i2, i4, i8
  use m_quadtree, only: QuadTree, QuadTreeNode, cellX, cellY, cellWidth, cellHeight, &
    NodeStack, initializeStack, deleteStack, push, pop
  use m_matrix_io, only: writeRealMatrixToH5Dataset, writeReal4TensorToH5Dataset
  use m_datastructures, only: QuadTreeParameters, initializeCellStats, SimulationParameters, &
    initializeStatisticsCell, deleteStatisticsCell, CellStats, StatisticsCell
  use m_util, only: k_B
  implicit none
contains

  subroutine writeQuadTreeToHDF5(tree, simParams, filename, datasetID, writeParticles)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(in) :: simParams
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: datasetID
    logical, intent(in) :: writeParticles

    real(fp), pointer, dimension(:,:) :: matrix
    real(fp), pointer, dimension(:,:,:,:) :: tensor
    integer(i4), pointer, dimension(:) :: particleStartPosition
    type(QuadTreeNode), pointer :: node
    integer(i4) :: i, numParticles, idx1, idx2, num
    real(fp) :: width, height
    character(len=5) :: datasetName

    integer(i4) :: error

    real(fp) :: n !> particle density
    real(fp) :: rho !> mass density
    real(fp) :: c0_x, c0_y, c0_z !> mean velocity components
    real(fp) :: cx_sq, cy_sq, cz_sq !> mean squared relative velocity components
    real(fp) :: c_sq !> mean squared relative velocity (relative = c - c_mean)
    real(fp) :: p !> pressure
    real(fp) :: T !> translational temperature


    width = tree%treeParams%width
    height = tree%treeParams%height
    allocate(matrix(12, tree%leafNumber))
    numParticles = 0
    allocate(particleStartPosition(tree%leafNumber))
    particleStartPosition = 1
    do i = 1, tree%leafNumber
      num = tree%particleNumbers(i)
      idx1 = tree%particleStartIndices(i)
      numParticles = numParticles + num
      particleStartPosition(i+1:) = particleStartPosition(i+1:) + num
      node => tree%leafs(i)%node

      !> calculate the node stats
      n = num / (simParams%V_c * cellWidth(node) * cellHeight(node)) * simParams%F_N
      !> TODO: update using the particle types
      rho = simParams%m(1) * n
      
      !> calculate the mean velocity components
      c0_x = sum(tree%particles(idx1:idx1+num-1, 3)) / num
      c0_y = sum(tree%particles(idx1:idx1+num-1, 4)) / num
      c0_z = sum(tree%particles(idx1:idx1+num-1, 5)) / num
      !> calculate the mean relative velocities
      cx_sq = sum((tree%particles(idx1:idx1+num-1, 3)-c0_x) ** 2)/num
      cy_sq = sum((tree%particles(idx1:idx1+num-1, 4)-c0_y) ** 2)/num
      cz_sq = sum((tree%particles(idx1:idx1+num-1, 5)-c0_z) ** 2)/num
      
      c_sq = cx_sq + cy_sq + cz_sq
      p = rho / 3 * c_sq
      !> translational temperature 3/2 k * T_tr = 1/2 m * \overbar{c^2}
      T = simParams%m(1) / k_B * c_sq / 3 
      

      matrix(:, i) = [cellX(node)*width, cellY(node)*height, cellWidth(node)*width, cellHeight(node)*height, &
        n, rho, p, T, c_sq, c0_x, c0_y, c0_z]
    end do 

    write(datasetName, "(i5.5)") datasetID

    call writeRealMatrixToH5Dataset(filename, datasetName//"_tree", matrix,tree%leafNumber, 12, error)
    deallocate(matrix)

    !> store the data from the statistics cells
    call createTensorFromStatisticsCells(tree, tensor)
    call writeReal4TensorToH5Dataset(filename, datasetName//"_stats", tensor, &
      [tree%treeParams%numParticleSpecies+1, 8, int(tree%treeParams%numStatisticsCellRows), &
       int(tree%treeParams%numStatisticsCellColumns)], &
      error)
    deallocate(tensor)

    if (writeParticles) then
      !> write the particle matrix to the file
      allocate(matrix(numParticles,5))
      !> TODO: can be parallelized
      !$OMP PARALLEL DO private(i,num,idx1,idx2) shared(matrix, tree)
      do i = 1, tree%leafNumber
        num = tree%particleNumbers(i)
        idx1 = tree%particleStartIndices(i)
        idx2 = particleStartPosition(i)
        matrix(idx2:idx2+num-1, :) = tree%particles(idx1:idx1+num-1, :)
      end do
      !$OMP END PARALLEL DO
      
      call writeRealMatrixToH5Dataset(filename, datasetName//"_particles", matrix, 5, numParticles, error)
      deallocate(matrix)
    end if 
    deallocate(particleStartPosition)
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
    type(CellStats), pointer :: stats

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
    allocate(tree%particles(params%elementChunkSize * nLeafs, 5))
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

    !> initialize the statistics cell array
    allocate(tree%statisticsCells(params%numStatisticsCellRows,params%numStatisticsCellColumns))
    do i = 1_i2, params%numStatisticsCellRows
      do j = 1_i2, params%numStatisticsCellColumns
        call initializeStatisticsCell(tree%statisticsCells(i,j), params)
      end do 
    end do 
    
  end subroutine createTreeFromFile 


  subroutine createTensorFromStatisticsCells(tree, tensor)
    implicit none
    type(QuadTree), pointer, intent(in) :: tree
    real(fp), dimension(:,:,:,:), pointer, intent(inout) :: tensor

    type(StatisticsCell), pointer :: stats
    integer(i2) :: i,j,nrows,ncols

    nrows = tree%treeParams%numStatisticsCellRows
    ncols = tree%treeParams%numStatisticsCellColumns

    allocate(tensor(tree%treeParams%numParticleSpecies+1, 8, nrows, ncols))
    
    do i = 1_i2, ncols
      do j = 1_i2, nrows
        stats => tree%statisticsCells(j,i)
        tensor(:,1,j,i) = real(stats%numParticles%average) 
        tensor(:,2,j,i) = stats%n%average
        tensor(:,3,j,i) = stats%rho%average
        tensor(:,4,j,i) = stats%cx_0%average
        tensor(:,5,j,i) = stats%cy_0%average
        tensor(:,6,j,i) = stats%cz_0%average
        tensor(:,7,j,i) = stats%p%average
        tensor(:,8,j,i) = stats%T%average
      end do 
    end do

    
  end subroutine createTensorFromStatisticsCells 


end module m_quadtree_io 
