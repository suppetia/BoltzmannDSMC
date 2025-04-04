module m_datastructures
  use m_types, only: fp, i4, i1, i2
  implicit none

  !> just stores the pointer to structures for a single node
  type :: StructuresPointer
    integer(i2) :: numStructures = 0
    real(fp), pointer, dimension(:,:) :: structures => null()
  end type StructuresPointer

  !> store the sequence of spawned particles during the simulation
  type :: SimulationSequence
    integer(i4), pointer, dimension(:) :: stepNumbers
    integer(i4), pointer, dimension(:,:) :: numNewParticles
    
    real(fp), pointer, dimension(:,:) :: spawnArea
    real(fp), pointer, dimension(:,:) :: velocityDistribution
  end type SimulationSequence

  !> store metadata about the simulation
  type :: SimulationParameters
    character(len=:), allocatable :: filebasename !> name of all files
    real(fp) :: width, height !> dimensions
    real(fp) :: V_c !> cell volume
    real(fp) :: F_N !> number of real particles per simulated particles
    real(fp) :: dt  !> time step per iteration
    real(fp) :: nu  !> parameter related to the temperature exponent omega proportional to viscosity mu ~ T^omega, nu = omega - 1/2
    integer(i1) :: numParticleSpecies !> number of used particle species
    real(fp), pointer, dimension(:) :: d_ref !> diameter of the particles (index represents the particle type)
    real(fp), pointer, dimension(:) :: T_ref !> reference temperature for the diameter of the particles (index represents the particle type)
    real(fp), pointer, dimension(:) :: m !> masses of the particles
    integer(i1) :: collisionModel !> 1: hard sphere, 2: variable hard sphere (VHS), 3: variable soft sphere (VSS)
    integer(i1) :: surfaceCollisionModel !> 1: specular, 2: diffuse
    real(fp) :: surfaceTemperature
    integer(i1) :: particleSpawnDistribution !> 1: equidistant, 2: uniform
    integer(i4) :: writeStatsFrequency, writeStatsOffset
    logical :: writeParticles
    integer(i4) :: writeParticlesFrequency, writeParticlesOffset
    integer(i4) :: numTimeSteps
    type(SimulationSequence), pointer :: simSeq
  end type SimulationParameters

  !> store all metadata about a quadtree
  type :: QuadTreeParameters
    real(fp) :: width, height !> (copy from SimulationParameters)
    integer(i4) :: elementSplitThreshold !> number of elements at which the node is split
    integer(i4) :: elementMergeThreshold !> number of elements at which the node is merged
    integer(i4) :: elementChunkSize !> for each node the memory for the elements is allocated in chunks

    integer(i2) :: cellHistoryLength

    integer(i1) :: numParticleSpecies !> number of used particle species (copy from SimulationParameters, required for the cell-stats)

    !> number of rows and columns for the statistics cells
    integer(i2) :: numStatisticsCellRows, numStatisticsCellColumns
  end type

  type :: ParticleList
    real(fp), pointer, dimension(:, :) :: particles
    integer(i1), pointer, dimension(:) :: particleTypes
    integer(i4) :: chunksize
    integer(i4) :: length, numParticles
  end type
  
  type :: ParticleAverageCounter
    integer(i2) :: historyLength
    integer(i4), dimension(:), allocatable :: history

    integer(i2) :: currentPosition
    real(fp) :: average
    logical :: historyIsFull = .false.
  end type

  type :: IntegerAverageCounter
    integer(i2) :: historyLength
    integer(i4), dimension(:,:), allocatable :: history

    integer(i2) :: currentPosition
    real(fp), dimension(:), allocatable :: average
    logical :: historyIsFull = .false.
  end type

  type :: RealAverageCounter
    integer(i2) :: historyLength
    real(fp), dimension(:,:), allocatable :: history

    integer(i2) :: currentPosition
    real(fp), dimension(:), allocatable :: average
    logical :: historyIsFull = .false.
  end type

  type :: RealAverageCounterArrayItem
    !> use this type as a workaround to store RealAverageCounter objects in an array
    type(RealAverageCounter), pointer :: counter
  end type

  type :: CellStats
    type(IntegerAverageCounter), pointer :: particleCounter
    type(RealAverageCounter), pointer :: statsCounter
    type(RealAverageCounterArrayItem), dimension(:), pointer :: speciesStatsCounter

    integer(i4), pointer :: N_avg !> average number of particles in the cell

    real(fp) :: maxSigmaC !> maxValue of sigma_T * c_r

    real(fp) :: n !> particle density
    real(fp) :: rho !> mass density
    real(fp) :: c0_x, c0_y, c0_z !> mean velocity components
    real(fp) :: cx_sq, cy_sq, cz_sq !> mean squared relative velocity components
    real(fp) :: c_sq !> mean squared relative velocity (relative = c - c_mean)
    real(fp) :: p !> pressure
    real(fp) :: T !> translational temperature

  end type CellStats


  type :: StatisticsCell
    !> in each of the counters the values are stored in the following layout
    !> [total value, value for species 1, value for species 2, ...]
    type(IntegerAverageCounter), pointer :: numParticles
    type(RealAverageCounter), pointer :: n !> number density
    type(RealAverageCounter), pointer :: rho !> mass density
    type(RealAverageCounter), pointer :: cx_0, cy_0, cz_0 !> mean velocity components
    type(RealAverageCounter), pointer :: p !> pressure
    type(RealAverageCounter), pointer :: T !> translational temperature
  end type StatisticsCell


contains

  subroutine initializeParticleList(list, chunksize)
    implicit none
    type(ParticleList), intent(out) :: list
    integer(i4), intent(in), optional :: chunksize

    if (present(chunksize)) then
      list%chunksize = chunksize
    else
      list%chunksize = 10000
    end if

    list%numParticles = 0
    list%length = list%chunksize
    allocate(list%particles(list%length, 5))
    allocate(list%particleTypes(list%length))
  end subroutine initializeParticleList

  subroutine extendParticleListLength(list)
    implicit none
    type(ParticleList), intent(inout) :: list
    
    real(fp), pointer, dimension(:,:) :: tmp
    integer(i1), pointer, dimension(:) :: tmpTypes

    !> temporarily copy the current particles to tmp
    tmp => list%particles
    nullify(list%particles)
    allocate(list%particles(list%length + list%chunksize, 5))

    tmpTypes => list%particleTypes
    nullify(list%particleTypes)
    allocate(list%particleTypes(list%length + list%chunksize))

    ! !> temporarily copy the current particles to tmp
    ! allocate(tmp(5,list%length))
    ! tmp(:, :) = list%particles(:, :)
    ! !> deallocate the particles array and reallocate with new size
    ! deallocate(list%particles)
    ! allocate(list%particles(5,list%length + list%chunksize))

    !> copy the particles back to the particles array
    list%particles(1:list%length,:) = tmp(:, :)
    list%particleTypes(1:list%length) = tmpTypes(:)

    list%length = list%length + list%chunksize

    deallocate(tmp)
    deallocate(tmpTypes)
  end subroutine extendParticleListLength

  subroutine deleteParticleList(list)
    implicit none
    type(ParticleList), intent(inout) :: list

    deallocate(list%particles)
    deallocate(list%particleTypes)
    list%numParticles = 0
    list%length = 0
    list%chunksize = 0
  end subroutine deleteParticleList

  
  subroutine append(list, particle, particleType)
    implicit none
    type(ParticleList), intent(inout) :: list
    real(fp), dimension(5), intent(in) :: particle
    integer(i1), intent(in) :: particleType

    if (list%numParticles == list%length) then
      call extendParticleListLength(list)
    end if

    list%numParticles = list%numParticles + 1
    list%particles(list%numParticles,:) = particle
    list%particleTypes(list%numParticles) = particleType
  end subroutine

  subroutine appendMultiple(list, particles, particleTypes)
    implicit none
    type(ParticleList), intent(inout) :: list
    real(fp), dimension(:,:), pointer, intent(in) :: particles
    integer(i1), dimension(:), pointer, intent(in) :: particleTypes

    if (.not.(size(particleTypes,1) == size(particles,1 ))) then
      print *, "particles array must have same length as particleTypes"
      return
    end if 

    if (.not.associated(particles)) then
      return
    end if
    if (list%numParticles + size(particles,1) >= list%length) then
      call extendParticleListLength(list)
    end if

    list%particles(list%numParticles:list%numParticles+size(particles,1),:) = particles(:,:)
    list%particleTypes(list%numParticles:list%numParticles+size(particles,1)) = particleTypes
    list%numParticles = list%numParticles + size(particles,1) 
  end subroutine

  subroutine initializeIntegerAverageCounter(counter, numValues, historyLength)
    implicit none
    type(IntegerAverageCounter), pointer, intent(out) :: counter
    integer(i2), intent(in) :: numValues, historyLength

    allocate(counter)
    allocate(counter%history(numValues, historyLength))
    allocate(counter%average(numValues))
    counter%history = 0
    counter%historyLength = historyLength
    counter%currentPosition = 0
    counter%average = 0
  end subroutine initializeIntegerAverageCounter

  subroutine initializeRealAverageCounter(counter, numValues, historyLength)
    implicit none
    type(RealAverageCounter), pointer, intent(out) :: counter
    integer(i2), intent(in) :: numValues, historyLength

    allocate(counter)
    allocate(counter%history(numValues, historyLength))
    allocate(counter%average(numValues))
    counter%history = 0
    counter%historyLength = historyLength
    counter%currentPosition = 0
    counter%average = 0
  end subroutine initializeRealAverageCounter

  subroutine initializeParticleAverageCounter(counter, historyLength)
    implicit none
    type(ParticleAverageCounter), pointer, intent(out) :: counter
    integer(i2), intent(in) :: historyLength

    allocate(counter)
    allocate(counter%history(historyLength))
    counter%history = 0
    counter%historyLength = historyLength
    counter%currentPosition = 0
    counter%average = 0
  end subroutine initializeParticleAverageCounter

  subroutine deleteParticleAverageCounter(counter)
    implicit none
    type(ParticleAverageCounter), pointer, intent(inout) :: counter

    deallocate(counter%history)
    deallocate(counter)
  end subroutine deleteParticleAverageCounter

  subroutine deleteIntegerAverageCounter(counter)
    implicit none
    type(IntegerAverageCounter), pointer, intent(inout) :: counter

    deallocate(counter%history)
    deallocate(counter%average)
    deallocate(counter)
  end subroutine deleteIntegerAverageCounter

  subroutine deleteRealAverageCounter(counter)
    implicit none
    type(RealAverageCounter), pointer, intent(inout) :: counter

    deallocate(counter%history)
    deallocate(counter%average)
    deallocate(counter)
  end subroutine deleteRealAverageCounter

  subroutine addIntegerCount(counter, val)
    implicit none
    type(IntegerAverageCounter), pointer, intent(inout) :: counter
    integer(i4), dimension(:), intent(in) :: val

    if (counter%currentPosition == counter%historyLength) then
      counter%historyIsFull = .true. !> if the position is reset to the beginning the history is filled at least once
      counter%currentPosition = 1_i2
    else
      counter%currentPosition = counter%currentPosition + 1_i2
    end if
    if (.not.counter%historyIsFull) then
      counter%history(:, counter%currentPosition) = val
      counter%average = real(sum(counter%history,2), fp)/counter%currentPosition
    else
      ! counter%average = counter%average - real(counter%history(:, counter%currentPosition), fp)/counter%historyLength
      counter%history(:, counter%currentPosition) = val
      ! counter%average = counter%average + real(val, fp)/counter%historyLength
      counter%average = real(sum(counter%history,2), fp)/counter%historyLength
    end if 
  end subroutine addIntegerCount

  subroutine addRealCount(counter, val)
    implicit none
    type(RealAverageCounter), pointer, intent(inout) :: counter
    real(fp), dimension(:), intent(in) :: val

    if (counter%currentPosition == counter%historyLength) then
      counter%historyIsFull = .true. !> if the position is reset to the beginning the history is filled at least once
      counter%currentPosition = 1_i2
    else
      counter%currentPosition = counter%currentPosition + 1_i2
    end if
    if (.not.counter%historyIsFull) then
      counter%history(:, counter%currentPosition) = val
      counter%average = sum(counter%history,2)/counter%currentPosition
    else
      ! counter%average = counter%average - counter%history(:, counter%currentPosition)/counter%historyLength
      counter%history(:, counter%currentPosition) = val
      ! counter%average = counter%average + val/counter%historyLength
      counter%average = sum(counter%history,2)/counter%historyLength
    end if 
  end subroutine addRealCount

  subroutine initializeCellStats(stats, treeParams)
    implicit none
    type(CellStats), pointer, intent(inout) :: stats
    type(QuadTreeParameters), intent(in) :: treeParams

    integer(i1) :: i

    allocate(stats)
    ! call initializeParticleAverageCounter(stats%particleCounter, treeParams%cellHistoryLength)
    call initializeIntegerAverageCounter(stats%particleCounter, 1_i2, treeParams%cellHistoryLength)
    call initializeRealAverageCounter(stats%statsCounter, 11_i2, treeParams%cellHistoryLength)
    allocate(stats%speciesStatsCounter(treeParams%numParticleSpecies))
    do i = 1, treeParams%numParticleSpecies
      call initializeRealAverageCounter(stats%speciesStatsCounter(i)%counter, 5_i2, treeParams%cellHistoryLength)
    end do

    stats%maxSigmaC = 1e-16_fp
  end subroutine initializeCellStats 

  subroutine deleteCellStats(stats)
    implicit none
    type(CellStats), pointer, intent(inout) :: stats

    integer(i1) :: i

    call deleteIntegerAverageCounter(stats%particleCounter)
    call deleteRealAverageCounter(stats%statsCounter)

    do i = 1, size(stats%speciesStatsCounter)
      call deleteRealAverageCounter(stats%speciesStatsCounter(i)%counter)
    end do
    deallocate(stats%speciesStatsCounter)
    deallocate(stats)
  end subroutine deleteCellStats 

  subroutine initializeStatisticsCell(stats, treeParams)
    implicit none
    type(StatisticsCell), intent(inout) :: stats
    type(QuadTreeParameters), pointer, intent(in) :: treeParams
    
    integer(i2) :: historyLength, numValues

    historyLength = treeParams%cellHistoryLength
    numValues = treeParams%numParticleSpecies + 1_i2
    call initializeIntegerAverageCounter(stats%numParticles, numValues, historyLength)
    call initializeRealAverageCounter(stats%n, numValues, historyLength)
    call initializeRealAverageCounter(stats%rho, numValues, historyLength)
    call initializeRealAverageCounter(stats%cx_0, numValues, historyLength)
    call initializeRealAverageCounter(stats%cy_0, numValues, historyLength)
    call initializeRealAverageCounter(stats%cz_0, numValues, historyLength)
    call initializeRealAverageCounter(stats%p, numValues, historyLength)
    call initializeRealAverageCounter(stats%T, numValues, historyLength)
  end subroutine initializeStatisticsCell 

  subroutine deleteStatisticsCell(stats)
    implicit none
    type(StatisticsCell), intent(inout) :: stats

    call deleteIntegerAverageCounter(stats%numParticles)
    call deleteRealAverageCounter(stats%n)
    call deleteRealAverageCounter(stats%rho)
    call deleteRealAverageCounter(stats%cx_0)
    call deleteRealAverageCounter(stats%cy_0)
    call deleteRealAverageCounter(stats%cz_0)
    call deleteRealAverageCounter(stats%p)
    call deleteRealAverageCounter(stats%T)
  end subroutine deleteStatisticsCell 

  subroutine deleteSimulationSequence(simSeq)
    implicit none
    type(SimulationSequence), pointer, intent(inout) :: simSeq

    deallocate(simSeq%stepNumbers)
    deallocate(simSeq%numNewParticles)
    deallocate(simSeq%stepNumbers)
    deallocate(simSeq%velocityDistribution)
  end subroutine deleteSimulationSequence 

  subroutine deleteSimulationParameters(params)
    implicit none
    type(SimulationParameters), pointer, intent(inout) :: params

    deallocate(params%d_ref)
    deallocate(params%m)
    deallocate(params%filebasename)
    call deleteSimulationSequence(params%simSeq)
  end subroutine deleteSimulationParameters 

end module m_datastructures
