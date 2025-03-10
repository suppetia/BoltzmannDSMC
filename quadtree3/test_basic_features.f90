program test_basic_features
  use m_types, only: fp, i4, i1
  use m_quadtree, only: QuadTreeNode, QuadTree, insertParticles, findParticleCells, splitNode, &
    NodeStack, splitNodes, updateTreeNodes, deleteTree
  use m_datastructures, only: QuadTreeParameters, SimulationParameters
  use m_quadtree_io, only: writeQuadTreeToHDF5, createTreeFromFile
  use m_simulation, only: moveParticles, step
  use m_util, only: gauss_random_2d

  use omp_lib
  implicit none

  type(QuadTree), pointer :: tree
  type(QuadTreeParameters), pointer :: params
  type(SimulationParameters), pointer :: simParams

  real(fp), dimension(:,:), pointer :: particles
  integer(i1), dimension(:), pointer :: particleTypes
  integer(i4), dimension(:), pointer :: leafIdx

  integer :: i,it, status
  integer(i4) :: tStart, tStop, numParticlesPerStep, numSimulationSteps, writeFrequency

  character(len=100) :: filebasename
  logical :: storeParticles

  ! filebasename = "data/test_18"
  filebasename = "data/wing1"
  numParticlesPerStep = 100
  writeFrequency = 10
  numSimulationSteps = 20
  storeParticles = .false.

  !$OMP PARALLEL
  !$OMP SINGLE
  print *, "number of threads:", omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL

  tStart = time()

  allocate(particles(50,2))
  call gauss_random_2d(particles)
  print *, particles
  deallocate(particles)


  !> initialize the quadtree-parameters
  allocate(params)
  params%width = 1._fp
  params%height = 1._fp
  params%elementSplitThreshold = 3
  params%elementMergeThreshold = 2
  params%elementChunkSize = 1
  params%cellHistoryLength = 30

  allocate(simParams)
  simParams%dt = 1e-4_fp
  simParams%numTimeSteps = numSimulationSteps
  simParams%writeFrequency = writeFrequency
  simParams%numParticleSpecies = 2
  allocate(simParams%m(simParams%numParticleSpecies))
  simParams%m(1) = 1e-25
  simParams%m(2) = 5e-25
  allocate(simParams%d_ref(simParams%numParticleSpecies))
  simParams%d_ref(1) = 5e-10
  simParams%d_ref(2) = 5e-9
  simParams%F_N = 1e14
  simParams%collisionModel = 1
  simParams%V_c = 1._fp

  !> copy this number before creating the tree so the initial cells have the correct size
  params%numParticleSpecies = simParams%numParticleSpecies 

  call createTreeFromFile(tree, params, trim(filebasename)//".tree")
  print *, tree%particleStartIndices(tree%leafNumber)
  simParams%width = tree%treeParams%width
  simParams%height = tree%treeParams%height
  print *, "num particle species:", tree%treeParams%numParticleSpecies 
  print *, tree%leafNumber
  ! simParams%V_c = simParams%width * simParams%height * 1._fp

  ! !> initializing the tree by hand
  ! allocate(tree)
  ! tree%treeParams => params
  ! allocate(tree%root)
  ! tree%root%nodeIdx = 1
  ! tree%leafNumber = 1
  ! allocate(tree%leafs(tree%leafNumber))
  ! tree%leafs(1)%node => tree%root
  ! allocate(tree%particles(params%elementChunkSize, 5))
  ! allocate(tree%particleTypes(params%elementChunkSize))
  ! allocate(tree%particleNumbers(tree%leafNumber))
  ! allocate(tree%particleStartIndices(tree%leafNumber))
  ! allocate(tree%structures(1))
  ! tree%particleNumbers(1) = 0
  ! tree%particleStartIndices(1) = 1

  !
  ! allocate(particles(1000000, 5))
  ! call random_number(particles)
  ! particles(:, 1) = particles(:, 1) * tree%treeParams%width
  ! particles(:, 2) = particles(:, 2) * tree%treeParams%height
  ! ! allocate(particles(2, 5))
  ! ! particles = reshape([[.1, .1, 0.,0.,0.],[.4,.6,0.,0.,0.]], [5,2])
  ! ! allocate(particles(1, 5))
  ! ! particles = reshape([[.1, .1, 0.,0.,0.]],[5,1])
  !
  !
  ! call findParticleCells(tree, particles, leafIdx)
  ! ! print *, leafIdx
  ! ! print *, tree%particles(:, 1)
  ! call insertParticles(tree, particles, leafIdx)
  ! print *, tree%particles(:, 1)
  ! if (associated(tree%leafs(1)%node%tmpParticles)) then
  !   print *, tree%leafs(1)%node%tmpParticles(:, 1)
  !   print *, size(tree%leafs(1)%node%tmpParticles, 2)
  !   print *, tree%leafs(1)%node%tmpParticleCount
  ! end if
  !
  ! ! call splitNode(tree%root, tree)
  ! call splitNodes(tree, stack)
  ! do i = 1, 4
  !   print *, tree%root%children(i)%tmpParticleCount
  !   ! print *, tree%root%children(1)%children(i)%tmpParticleCount
  ! end do
  !
  ! do i=1,size(tree%leafs)
  !   print *, associated(tree%leafs(i)%node)
  ! end do
  
  call updateTreeNodes(tree, simParams)
  ! do i=1,size(tree%leafs)
  !   print *, associated(tree%leafs(i)%node%tmpParticles)
  !   print *,tree%leafs(i)%node%tmpParticleCount
  ! end do
  ! do i =1,4
  !   print *, tree%root%children(i)%nodeIdx
  ! end do
  
  ! do i = 1, size(tree%leafs) /2
  !   ! tree%leafs(i)%node%tmpParticleCount = 0
  !   tree%particleNumbers(i) = 0
  ! end do
  
  ! call updateTreeNodes(tree)

  !> delete the last file
  open(unit=1234, iostat=status, file=trim(filebasename)//".h5", status='old')
  if (status == 0) close(1234, status='delete')

  call writeQuadTreeToHDF5(tree, simParams, trim(filebasename)//".h5", 0, storeParticles)
  ! print *, particles(:, 1)


  allocate(particles(numParticlesPerStep,5))
  allocate(particleTypes(numParticlesPerStep))
  particleTypes(:numParticlesPerStep/2) = 1
  particleTypes(numParticlesPerStep/2:) = 2
  do it = 1, simParams%numTimeSteps
    print *, it
    call random_number(particles)
    particles(:, 1) = 0._fp
    particles(:, 2) = particles(:, 2) * tree%treeParams%height
    particles(:, 3) = (10._fp + particles(:, 3) * 10._fp) * (tree%treeParams%height+tree%treeParams%width)/2
    particles(:, 4:) = (particles(:, 4:) -.5_fp) * 10._fp * (tree%treeParams%height+tree%treeParams%width)/2
    ! !> for the "rocket engine"
    ! particles(:, 1) = 400._fp + (particles(:,1) - .5_fp) * 30
    ! particles(:, 2) = 1550._fp + (particles(:,2) - .5_fp) * 30
    ! particles(:, 3) = (particles(:, 3) -.5_fp)* 10._fp *3000!* (tree%treeParams%height+tree%treeParams%width)/2
    ! particles(:, 4:) = (particles(:, 4:) -.5_fp) * 10._fp *3000!* (tree%treeParams%height+tree%treeParams%width)/2
    ! particles(:,1) = 0._fp
    ! particles(:,2) = particles(:, 2)* 100
    ! particles(:,3) = 1e4_fp
    ! particles(:,4) = 0._fp
    call findParticleCells(tree, particles, leafIdx)
    ! print *, leafIdx, tree%leafNumber
    if (it < 100) then
      call insertParticles(tree, particles, particleTypes, leafIdx)
      call updateTreeNodes(tree, simParams)
    end if 
    ! call insertParticles(tree, particles, particleTypes, leafIdx)
    ! call updateTreeNodes(tree, simParams)

    ! call random_number(particles)
    ! particles(:, 1) = particles(:, 1) * tree%treeParams%width
    ! particles(:, 2) = particles(:, 2) * tree%treeParams%height
    ! call findParticleCells(tree, particles, leafIdx)
    ! call insertParticles(tree, particles, leafIdx)
    ! call moveParticles(tree, simParams)
    ! call updateTreeNodes(tree, simParams)
    call step(tree, simParams)
    if (mod(it, simParams%writeFrequency) == 0) then
      call writeQuadTreeToHDF5(tree, simParams, trim(filebasename)//".h5", it, storeParticles)
    end if 
  end do 

  call deleteTree(tree)

  deallocate(particles)
  deallocate(params)
  deallocate(simParams%m)
  deallocate(simParams%d_ref)
  deallocate(simParams)

  tStop = time()
  open(unit=1234, iostat=status, file=trim(filebasename)//".sum", status='old')
  write(1234, *) tStop-tStart
  
end program test_basic_features
