module m_datastructures
  use m_types, only: fp, i4, i1, i2
  implicit none

  !> just stores the pointer to structures for a single node
  type :: StructuresPointer
    integer(i2) :: numStructures = 0
    real(fp), pointer, dimension(:,:) :: structures => null()
  end type StructuresPointer

  !> store metadata about the simulation
  type :: SimulationParameters
    real(fp) :: width, height !> dimensions
    real(fp) :: V_c !> cell volume
    real(fp) :: F_N !> number of real particles per simulated particles
    real(fp) :: dt  !> time step per iteration
    real(fp), pointer, dimension(:) :: d_ref !> diameter of the particles (index represents the particle type)
    real(fp), pointer, dimension(:) :: m !> masses of the particles
    integer(i1) :: collisionModel !> 1: hard sphere, 2: variable hard sphere (VHS), 3: variable soft sphere (VSS)
    integer(i4) :: writeFrequency
    integer(i4) :: numTimeSteps
  end type SimulationParameters

  !> store all metadata about a quadtree
  type :: QuadTreeParameters
    real(fp) :: width, height
    integer(i4) :: elementSplitThreshold !> number of elements at which the node is split
    integer(i4) :: elementMergeThreshold !> number of elements at which the node is merged
    integer(i4) :: elementChunkSize !> for each node the memory for the elements is allocated in chunks

    integer(i2) :: cellHistoryLength
  end type

  type :: ParticleList
    real(fp), pointer, dimension(:, :) :: particles
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

  type :: CellStats
    type(ParticleAverageCounter), pointer :: particleCounter
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
  end subroutine initializeParticleList

  subroutine extendParticleListLength(list)
    implicit none
    type(ParticleList), intent(inout) :: list
    
    real(fp), pointer, dimension(:,:) :: tmp

    !> temporarily copy the current particles to tmp
    tmp => list%particles
    nullify(list%particles)
    allocate(list%particles(list%length + list%chunksize, 5))


    ! !> temporarily copy the current particles to tmp
    ! allocate(tmp(5,list%length))
    ! tmp(:, :) = list%particles(:, :)
    ! !> deallocate the particles array and reallocate with new size
    ! deallocate(list%particles)
    ! allocate(list%particles(5,list%length + list%chunksize))

    !> copy the particles back to the particles array
    list%particles(1:list%length,:) = tmp(:, :)
    list%length = list%length + list%chunksize

    deallocate(tmp)
  end subroutine extendParticleListLength

  subroutine deleteParticleList(list)
    implicit none
    type(ParticleList), intent(inout) :: list

    deallocate(list%particles)
    list%numParticles = 0
    list%length = 0
    list%chunksize = 0
  end subroutine deleteParticleList

  
  subroutine append(list, particle)
    implicit none
    type(ParticleList), intent(inout) :: list
    real(fp), dimension(5), intent(in) :: particle

    if (list%numParticles == list%length) then
      call extendParticleListLength(list)
    end if

    list%numParticles = list%numParticles + 1
    list%particles(list%numParticles,:) = particle(:)
  end subroutine

  subroutine appendMultiple(list, particles)
    implicit none
    type(ParticleList), intent(inout) :: list
    real(fp), dimension(:,:), pointer, intent(in) :: particles

    if (.not.associated(particles)) then
      return
    end if
    if (list%numParticles + size(particles,1) >= list%length) then
      call extendParticleListLength(list)
    end if

    list%particles(list%numParticles:list%numParticles+size(particles,1),:) = particles(:,:)
    list%numParticles = list%numParticles + size(particles,1) 
  end subroutine


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

  subroutine addParticleCount(counter, val)
    implicit none
    type(ParticleAverageCounter), pointer, intent(inout) :: counter
    integer(i4), intent(in) :: val

    if (counter%currentPosition == counter%historyLength) then
      counter%historyIsFull = .true. !> if the position is reset to the beginning the history is filled at least once
      counter%currentPosition = 1_i2
    else
      counter%currentPosition = counter%currentPosition + 1_i2
    end if
    if (.not.counter%historyIsFull) then
      counter%history(counter%currentPosition) = val
      counter%average =  real(sum(counter%history), fp)/counter%currentPosition
    else
      counter%average = counter%average - real(counter%history(counter%currentPosition), fp)/counter%historyLength
      counter%history(counter%currentPosition) = val
      counter%average = counter%average + real(val, fp)/counter%historyLength
    end if 
  end subroutine addParticleCount

  subroutine initializeCellStats(stats, treeParams)
    implicit none
    type(CellStats), pointer, intent(inout) :: stats
    type(QuadTreeParameters), intent(in) :: treeParams

    allocate(stats)
    call initializeParticleAverageCounter(stats%particleCounter, treeParams%cellHistoryLength)

    stats%maxSigmaC = 1e-16_fp
  end subroutine initializeCellStats 

  subroutine deleteCellStats(stats)
    implicit none
    type(CellStats), pointer, intent(inout) :: stats

    call deleteParticleAverageCounter(stats%particleCounter)
    deallocate(stats)
  end subroutine deleteCellStats 

end module m_datastructures
