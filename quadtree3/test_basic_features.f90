program test_basic_features
  use m_types, only: fp, i4
  use m_quadtree, only: QuadTreeNode, QuadTree, insertParticles, findParticleCells, splitNode, &
    NodeStack, splitNodes, updateTreeNodes, deleteTree
  use m_datastructures, only: QuadTreeParameters, SimulationParameters
  use m_quadtree_io, only: writeQuadTreeToHDF5, createTreeFromFile
  use m_simulation, only: moveParticles, step

  use omp_lib
  implicit none

  type(QuadTree), pointer :: tree
  type(QuadTreeParameters), pointer :: params
  type(SimulationParameters), pointer :: simParams

  real(fp), dimension(:,:), pointer :: particles
  integer(i4), dimension(:), pointer :: leafIdx

  integer :: i,it, status

  character(len=100) :: filebasename
  logical :: storeParticles

  filebasename = "data/test_18"
  storeParticles = .true.

  !$OMP PARALLEL
  print *, "number of threads:", omp_get_num_threads()
  !$OMP END PARALLEL


  !> initialize the quadtree-parameters
  allocate(params)
  params%width = 1._fp
  params%height = 1._fp
  params%elementSplitThreshold = 50
  params%elementMergeThreshold = 40
  params%elementChunkSize = 10
  params%cellHistoryLength = 30

  allocate(simParams)
  simParams%dt = 1e-4_fp
  simParams%numTimeSteps = 2000
  simParams%writeFrequency = 40
  allocate(simParams%m(1))
  simParams%m(1) = 1e-25
  allocate(simParams%d_ref(1))
  simParams%d_ref(1) = 5e-10
  simParams%F_N = 1e14
  simParams%collisionModel = 1
  simParams%V_c = 1._fp


  call createTreeFromFile(tree, params, trim(filebasename)//".tree")
  simParams%width = tree%treeParams%width
  simParams%height = tree%treeParams%height
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

  call writeQuadTreeToHDF5(tree, trim(filebasename)//".h5", 0, storeParticles)
  ! print *, particles(:, 1)


  allocate(particles(1000,5))
  do it = 1, simParams%numTimeSteps
    print *, it
    call random_number(particles)
    particles(:, 1) = 0._fp
    particles(:, 2) = particles(:, 2) * tree%treeParams%height
    particles(:, 3) = (10._fp + particles(:, 3) * 10._fp) * (tree%treeParams%height+tree%treeParams%width)/2
    particles(:, 4:) = (particles(:, 4:) -.5_fp) * 10._fp * (tree%treeParams%height+tree%treeParams%width)/2
    call findParticleCells(tree, particles, leafIdx)
    call insertParticles(tree, particles, leafIdx)

    ! call random_number(particles)
    ! particles(:, 1) = particles(:, 1) * tree%treeParams%width
    ! particles(:, 2) = particles(:, 2) * tree%treeParams%height
    ! call findParticleCells(tree, particles, leafIdx)
    ! call insertParticles(tree, particles, leafIdx)
    ! call moveParticles(tree, simParams)
    ! call updateTreeNodes(tree, simParams)
    call step(tree, simParams)
    if (mod(it, simParams%writeFrequency) == 0) then
      call writeQuadTreeToHDF5(tree, trim(filebasename)//".h5", it, storeParticles)
    end if 
  end do 

  call deleteTree(tree)

  deallocate(particles)
  deallocate(params)
  deallocate(simParams)

end program test_basic_features
