program run_simulation
  use m_types, only: fp, i4, i1
  use m_quadtree, only: QuadTree, deleteTree, insertParticles, findParticleCells, updateTreeNodes
  use m_datastructures, only: QuadTreeParameters, SimulationParameters, SimulationSequence
  use m_quadtree_io, only: writeQuadTreeToHDF5, createTreeFromFile, createTensorFromStatisticsCells
  use m_simulation, only: step
  use m_simulation_io, only: parseSimulationSequence, parseSimulationConfiguration, setupSimulation
  use m_util, only: gauss_random_2d

  use omp_lib
  implicit none

  type(QuadTree), pointer :: tree
  type(SimulationParameters), pointer :: simParams


  real(fp), pointer, dimension(:,:) :: particles, ptrVelocities
  integer(i1), pointer, dimension(:) :: particleTypes
  integer(i4), dimension(:), pointer :: leafIdx
  integer(i4) :: it, currentSimState,j,k,n,m,numNewParticles,status, numX, numY
  logical :: writeParticles
  real(fp) :: spawnArea

  character(len=200) :: simulationFilename

  call get_command_argument(1, simulationFilename)
  if (len_trim(simulationFilename) == 0) then
    print *, "No simulation configuration file passed. Stopping the simulation."
    call exit(-1)
  end if 

  call setupSimulation(trim(simulationFilename), tree, simParams)
  nullify(particles)
  nullify(particleTypes)

  !> delete the last file
  open(unit=1234, iostat=status, file=simParams%filebasename//".h5", status='old')
  if (status == 0) close(1234, status='delete')
  call writeQuadTreeToHDF5(tree, simParams, simParams%filebasename//".h5", 0, &
    simParams%writeParticles .and. simParams%writeParticlesOffset == 0)

  print *, "store particles:", simParams%writeParticles
  print *, "frequency to write the stats:", simParams%writeStatsFrequency

  !> perform the simulation sequence
  currentSimState = 0
  numNewParticles = 0
  do it = 1, simParams%numTimeSteps
    print *, "step:", it
    if (it == simParams%simSeq%stepNumbers(currentSimState+1)) then
      currentSimState = currentSimState + 1
      if (associated(particles)) then
        deallocate(particles)
        deallocate(particleTypes)
      end if 
      numNewParticles = sum(simParams%simSeq%numNewParticles(:,currentSimState)) 
      if (numNewParticles > 0) then
        allocate(particles(numNewParticles, 5))
        allocate(particleTypes(numNewParticles))
      end if 
      if (simParams%particleSpawnDistribution == 1) then !> equidistant particle spawn location
        n = 1
        do j = 1, simParams%numParticleSpecies
          !> calculate the positions of the particles
          numNewParticles = simParams%simSeq%numNewParticles(j,currentSimState)
          print *, numNewParticles
          if (simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState) == 0) then
            particles(n:n+numNewParticles-1, 1) = simParams%simSeq%spawnArea(4*(j-1)+1, currentSimState)
            particles(n:n+numNewParticles-1, 2) = (/ (simParams%simSeq%spawnArea(4*(j-1)+2, currentSimState) &
              + real(k,fp)/numNewParticles * simParams%simSeq%spawnArea(4*(j-1)+4, currentSimState)&
              , k=0,numNewParticles-1) /) 
          else if (simParams%simSeq%spawnArea(4*(j-1)+4, currentSimState) == 0) then
            particles(n:n+numNewParticles-1, 1) = (/ (simParams%simSeq%spawnArea(4*(j-1)+1, currentSimState) &
              + real(k,fp)/numNewParticles * simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState) &
              , k=0,numNewParticles-1) /)
            particles(n:n+numNewParticles-1, 2) = simParams%simSeq%spawnArea(4*(j-1)+2, currentSimState)
          else
            spawnArea = simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState) &
            * simParams%simSeq%spawnArea(4*(j-1)+4, currentSimState)
            numX = simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState) * sqrt(numNewParticles / spawnArea)
            numY = simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState) * sqrt(numNewParticles / spawnArea)

            do m = 1, numY
              particles(n+(m-1)*numX:n+m*numX-1, 1) = (/ (simParams%simSeq%spawnArea(4*(j-1)+1, currentSimState) &
                + real(k,fp)/numX* simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState) &
                , k=0,numX-1) /)
              particles(n+(m-1)*numX:n+m*numX-1, 2) = simParams%simSeq%spawnArea(4*(j-1)+2, currentSimState) &
                + (m-1)/numY * simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState)
            end do 
          end if 
          n = n + numNewParticles
        end do
        numNewParticles = n-1 
        ! print *, particles(:, 2)
        ! print *, particles(:, 3)
      end if
    end if 
    
    if (numNewParticles > 0) then
      if (simParams%particleSpawnDistribution == 2) then
        !> spawn location is uniformly distributed
        call random_number(particles(:, 1:2))
      end if 
      !> velocity distribution is a normal distribution
      ptrVelocities => particles(:, 3:5)
      call gauss_random_2d(ptrVelocities)
      n = 1
      do j = 1, simParams%numParticleSpecies
        m = simParams%simSeq%numNewParticles(j, currentSimState)
        !> assign the particleTypes
        particleTypes(n:n+m-1) = j
        if (simParams%particleSpawnDistribution == 2) then
          !> rescale the positions
          particles(n:n+m-1, 1) = simParams%simSeq%spawnArea(4*(j-1)+1, currentSimState) &
            + particles(n:n+m-1,1) * simParams%simSeq%spawnArea(4*(j-1)+3, currentSimState)
          particles(n:n+m-1, 2) = simParams%simSeq%spawnArea(4*(j-1)+2, currentSimState) &
            + particles(n:n+m-1, 2) * simParams%simSeq%spawnArea(4*(j-1)+4, currentSimState)
        end if 
        !> rescale the velocities
        particles(n:n+m-1, 3) = simParams%simSeq%velocityDistribution(6*(j-1)+1, currentSimState) &
          + particles(n:n+m-1, 3) * simParams%simSeq%velocityDistribution(6*(j-1)+4, currentSimState)
        particles(n:n+m-1, 4) = simParams%simSeq%velocityDistribution(6*(j-1)+2, currentSimState) &
          + particles(n:n+m-1, 4) * simParams%simSeq%velocityDistribution(6*(j-1)+5, currentSimState)
        particles(n:n+m-1, 5) = simParams%simSeq%velocityDistribution(6*(j-1)+3, currentSimState) &
          + particles(n:n+m-1, 5) * simParams%simSeq%velocityDistribution(6*(j-1)+6, currentSimState)
        ! print *, "hi", j, particles(n, :)
        ! print *, m, simParams%simSeq%spawnArea(4*(j-1)+1:4*j, currentSimState)
        ! print *, m, simParams%simSeq%velocityDistribution(6*(j-1)+1:6*j, currentSimState)
  
        n = n+m
      end do 
      ! print *, particles(:, 1)
      ! print *, particles(:, 2)
      ! print *, particles(:, 3)
  
      call findParticleCells(tree, particles, leafIdx)
      call insertParticles(tree, particles, particleTypes, leafIdx)
      call updateTreeNodes(tree, simParams)
    end if 

    call step(tree, simParams)
    if (it >= simParams%writeStatsOffset .and. &
      mod(it-simParams%writeStatsOffset, simParams%writeStatsFrequency) == 0) then
      
      !> only write the particles if the offset is crossed and at the given frequency
      if (simParams%writeParticles .and. it >= simParams%writeParticlesOffset .and. &
        mod(it-simParams%writeParticlesOffset, simParams%writeParticlesFrequency) == 0) then

        writeParticles = .true.
      else
        writeParticles = .false.
      end if 
      call writeQuadTreeToHDF5(tree, simParams, simParams%filebasename//".h5", it, &
        writeParticles)
    end if 
  end do 

  if (associated(particles)) then
    deallocate(particles)
    deallocate(particleTypes)
  end if 

  call deleteTree(tree)

  

end program run_simulation

