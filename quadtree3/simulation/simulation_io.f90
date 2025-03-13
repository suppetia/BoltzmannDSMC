module m_simulation_io  !
  use m_types, only: fp, i4, i1
  use m_datastructures, only: SimulationParameters, SimulationSequence, QuadTreeParameters
  use m_quadtree, only: QuadTree, updateTreeNodes, findParticleCells, insertParticles
  use m_quadtree_io, only: createTreeFromFile
  use m_matrix_io, only: readParticleMatrixFromFile
  implicit none
contains

  subroutine parseSimulationSequence(filename, sim)
    !> the simulation is stored in a file in the following format
    !> #step #particles spawnArea               meanVelocity  stdVelocity
    !> i     N1 N2 ...  (x y w h) * numSpecies  (v_x v_y v_z  s_x s_y s_z)*numSpecies
    !> j     N1 N2 ...  (x y w h) * numSpecies  (v_x v_y v_z  s_x s_y s_z)*numSpecies
    !> 
    !> which means that between step i and step j each step are Nk particles of species k
    !> spawned in the corresponding area with the corresponding velocity distribution
    implicit none
    character(len=*), intent(in) :: filename
    type(SimulationSequence), pointer :: sim

    integer(i4) :: numRows, numParticleSpecies, i

    allocate(sim)

    open(unit=10, file=filename, status='old', action='read')

    !> read the number of rows in the first line
    read(10, *) numRows, numParticleSpecies

    !> allocate the arrays in sim based on that info
    allocate(sim%stepNumbers(numRows+1))
    allocate(sim%numNewParticles(numParticleSpecies, numRows+1))
    allocate(sim%spawnArea(4*numParticleSpecies, numRows+1))
    allocate(sim%velocityDistribution(6*numParticleSpecies, numRows+1))
    !> store dummy values in the additional line where only the stepNumbers entry is used as a comparison value
    sim%stepNumbers(numRows+1) = 2147483647 !> max int4 value
    sim%numNewParticles(:, numRows+1) = 0
    sim%spawnArea(:, numRows+1) = 0._fp 
    sim%velocityDistribution(:, numRows+1) = 0._fp

    do i = 1, numRows
      read(10, *) sim%stepNumbers(i), sim%numNewParticles(:, i), &
        sim%spawnArea(:, i), sim%velocityDistribution(:,i)
    end do 

    close(10)
    
  end subroutine parseSimulationSequence

  subroutine parseSimulationConfiguration(filename, simParams, treeParams)
    use cfgio_mod, only: cfg_t, parse_cfg
    implicit none
    character(len=*), intent(in) :: filename
    type(SimulationParameters), pointer, intent(inout) :: simParams
    type(QuadTreeParameters), pointer, intent(inout) :: treeParams

    type(cfg_t) :: cfg
    character(len=100) :: tmpFilebasename
    integer(i4) :: tmpInt4
    real(fp), allocatable, dimension(:) :: tmpDoubleArray

    allocate(simParams)
    allocate(treeParams)

    cfg = parse_cfg(filename)
    call cfg%get("", "filebasename", tmpFilebasename)
    simParams%filebasename = trim(tmpFilebasename)
    call cfg%get("Simulation", "numTimeSteps", simParams%numTimeSteps)
    call cfg%get("Simulation", "writeStatsFrequency", simParams%writeStatsFrequency)
    if (cfg%has_key("Simulation", "writeStatsOffset")) then
      call cfg%get("Simulation", "writeStatsOffset", simParams%writeStatsOffset)
    else
      simParams%writeStatsOffset = 0
    end if 
    call cfg%get("Simulation", "writeParticles", simParams%writeParticles)
    if (simParams%writeParticles) then
      call cfg%get("Simulation", "writeParticlesOffset", simParams%writeParticlesOffset)
      call cfg%get("Simulation", "writeParticlesFrequency", simParams%writeParticlesFrequency)
      if (mod(simParams%writeParticlesFrequency, simParams%writeStatsFrequency) /= 0 &
        .or. mod(simParams%writeParticlesOffset-simParams%writeStatsOffset, simParams%writeStatsFrequency) /= 0) then
        
        print *, "Note that the particles will only be written to the file if also the stats are written. &
          &Make sure to align the offsets and the write frequencies."
      end if 
    else
      simParams%writeParticlesOffset = 0
      simParams%writeParticlesFrequency = 0
    end if 
    call cfg%get("Simulation", "particleSpawnDistribution", tmpInt4)
    simParams%particleSpawnDistribution = tmpInt4
    call cfg%get("Simulation", "dt", simParams%dt)
    call cfg%get("Simulation", "F_N", simParams%F_N)
    call cfg%get("Simulation", "V_c", simParams%V_c)

    !> quadtree parameters
    call cfg%get("Quadtree", "splitThreshold", treeParams%elementSplitThreshold)
    call cfg%get("Quadtree", "mergeThreshold", treeParams%elementMergeThreshold)
    call cfg%get("Quadtree", "elementsChunkSize", treeParams%elementChunkSize)
    call cfg%get("Quadtree", "cellHistoryLength", tmpInt4)
    treeParams%cellHistoryLength = tmpInt4
    call cfg%get("MolecularModel", "numParticleSpecies", tmpInt4)
    treeParams%numParticleSpecies = tmpInt4
    simParams%numParticleSpecies = tmpInt4
    call cfg%get("Quadtree", "numStatisticsCellRows", tmpInt4)
    treeParams%numStatisticsCellRows = tmpInt4
    call cfg%get("Quadtree", "numStatisticsCellColumns", tmpInt4)
    treeParams%numStatisticsCellColumns = tmpInt4

    !> parameters of the MolecularModel stored in simParams
    call cfg%get("MolecularModel", "nu", simParams%nu)
    call cfg%get("MolecularModel", "collisionModel", tmpInt4)
    simParams%collisionModel = tmpInt4
    allocate(simParams%m(treeParams%numParticleSpecies))
    allocate(simParams%d_ref(treeParams%numParticleSpecies))
    allocate(simParams%T_ref(treeParams%numParticleSpecies))
    call cfg%get("MolecularModel", "m", tmpDoubleArray)
    simParams%m = tmpDoubleArray
    call cfg%get("MolecularModel", "d_ref", tmpDoubleArray)
    simParams%d_ref = tmpDoubleArray
    call cfg%get("MolecularModel", "T_ref", tmpDoubleArray)
    simParams%T_ref = tmpDoubleArray
    deallocate(tmpDoubleArray)

  end subroutine parseSimulationConfiguration 

  subroutine setupSimulation(filename, tree, simParams)
    implicit none
    character(len=*), intent(in) :: filename
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(inout) :: simParams

    type(SimulationSequence), pointer :: simSeq
    type(QuadTreeParameters), pointer :: treeParams
    real(fp), pointer, dimension(:,:) :: particles
    integer(i4), pointer, dimension(:):: leafIdx
    integer(i1), pointer, dimension(:):: particleTypes

    allocate(treeParams)
    allocate(simParams)

    call parseSimulationConfiguration(filename, simParams, treeParams)
    call parseSimulationSequence(simParams%filebasename//".sim", simSeq)

    call createTreeFromFile(tree, treeParams, simParams%filebasename//".tree")
    simParams%width = tree%treeParams%width
    simParams%height = tree%treeParams%height
    simParams%simSeq => simSeq

    !> read in the initial particles and insert them into the tree
    call readParticleMatrixFromFile(simParams%filebasename//".particles", particles, particleTypes)
    call findParticleCells(tree, particles, leafIdx)
    call insertParticles(tree, particles, particleTypes, leafIdx)
    
    call updateTreeNodes(tree, simParams)
    
  end subroutine setupSimulation 
  
end module m_simulation_io 
