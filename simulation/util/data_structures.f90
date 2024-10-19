module m_data_structures
  use m_types, only: dp, i4, i2
  implicit none

  type :: ParticleList
    real(dp), allocatable, dimension(:) :: particles
    integer(i4) :: chunksize
    integer(i4) :: length, numParticles
  end type

  type :: ParticleAverageCounter
    integer(i2) :: historyLength
    integer(i4), dimension(:), allocatable :: history

    integer(i2) :: currentPosition
    real(dp) :: average
    logical :: historyIsFull = .false.
  end type

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
    allocate(list%particles(list%length*5))
  end subroutine initializeParticleList

  subroutine extendParticleListLength(list)
    implicit none
    type(ParticleList), intent(inout) :: list
    
    real(dp), allocatable, dimension(:) :: tmp

    !> temporarily copy the current particles to tmp
    allocate(tmp(list%length*5))
    tmp(:) = list%particles(:)
    !> deallocate the particles array and reallocate with new size
    deallocate(list%particles)
    allocate(list%particles((list%length + list%chunksize) * 5))

    !> copy the particles back to the particles array
    list%particles(1:list%length*5) = tmp(:)
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
    real(dp), dimension(5), intent(in) :: particle

    if (list%numParticles == list%length) then
      call extendParticleListLength(list)
    end if

    list%particles(5*list%numParticles+1:5*list%numParticles+5) = particle(:)
    list%numParticles = list%numParticles + 1
  end subroutine

  subroutine getByIndex(list, idx, particle)
    implicit none
    type(ParticleList), intent(in) :: list
    integer(i4), intent(in) :: idx
    real(dp), dimension(5), intent(out) :: particle

    if (idx > list%numParticles) then
      write(*, fmt="(a,1x,i0)", advance="no") "invalid index:", idx
      return
    end if
    particle(:) = list%particles(5*idx-4:5*idx)

  end subroutine getByIndex

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
      counter%currentPosition = 1
    else
      counter%currentPosition = counter%currentPosition + 1
    end if
    if (.not.counter%historyIsFull) then
      counter%history(counter%currentPosition) = val
      counter%average =  real(sum(counter%history), dp)/counter%currentPosition
    else
      counter%average = counter%average - real(counter%history(counter%currentPosition), dp)/counter%historyLength
      counter%history(counter%currentPosition) = val
      counter%average = counter%average + real(val, dp)/counter%historyLength
    end if 

  end subroutine addParticleCount

end module m_data_structures
