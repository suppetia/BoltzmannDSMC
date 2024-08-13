module m_data_structures
  use m_types, only: dp, i4
  implicit none

  type :: ParticleList
    real(dp), allocatable, dimension(:) :: particles
    integer(i4) :: chunksize
    integer(i4) :: length, numParticles
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
    allocate(list%particles(list%length*4))
  end subroutine initializeParticleList

  subroutine extendParticleListLength(list)
    implicit none
    type(ParticleList), intent(inout) :: list
    
    real(dp), allocatable, dimension(:) :: tmp

    !> temporarily copy the current particles to tmp
    allocate(tmp(list%length*4))
    tmp(:) = list%particles(:)
    !> deallocate the particles array and reallocate with new size
    deallocate(list%particles)
    allocate(list%particles((list%length + list%chunksize) * 4))

    !> copy the particles back to the particles array
    list%particles(1:list%length*4) = tmp(:)
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
    real(dp), dimension(4), intent(in) :: particle

    if (list%numParticles == list%length) then
      call extendParticleListLength(list)
    end if

    list%particles(4*list%numParticles+1:4*list%numParticles+4) = particle(:)
    list%numParticles = list%numParticles + 1
  end subroutine

  subroutine getByIndex(list, idx, particle)
    implicit none
    type(ParticleList), intent(in) :: list
    integer(i4), intent(in) :: idx
    real(dp), dimension(4), intent(out) :: particle

    if (idx > list%numParticles) then
      write(*, fmt="(a,1x,i0)", advance="no") "invalid index:", idx
      return
    end if
    particle(:) = list%particles(4*idx-3:4*idx)

  end subroutine getByIndex
end module m_data_structures
