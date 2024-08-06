module m_simulation
  use m_quadtree, only: QuadTree, QuadTreeNode, removeElementFromNode, insertElement
  use m_quadtree_io, only: NodeStack, initializeStack, findLeaves, pop
  use m_data_structures, only: ParticleList, initializeParticleList, append, getByIndex, deleteParticleList
  use m_types, only: dp, i4

  implicit none
contains
  
  subroutine moveParticles(tree, dt)
    implicit none
    type(QuadTree), intent(inout) :: tree
    real(dp), intent(in) :: dt

    type(NodeStack) :: stack
    type(QuadTreeNode), pointer :: n

    integer(i4) :: numLeaves, i, j

    type(ParticleList) :: list
    real(dp), dimension(4) :: particle

    call initializeStack(stack)
    call initializeParticleList(list)

    call findLeaves(tree%root, stack)

    numLeaves = stack%topIndex

    do i = 1, numLeaves
      call pop(stack, n)

      j = n%numberOfElements
      !> update x components: x_new = x + dt * v_x
      n%elements(1:4*j:4) = n%elements(1:4*j:4) + dt * n%elements(3:4*j:4)
      !> update y components: y_new = y + dt*v_y
      n%elements(2:4*j:4) = n%elements(2:4*j:4) + dt * n%elements(4:4*j:4)

      !> check if the particle is still in the same cell after the location update
      do j = n%numberOfElements, 1, -1
        !> iterate backwards because of how the particles are removed:
        !> if the particle at index i is removed, it is replaced with the last particle
        if ( &
          n%elements(4*j-3) < n%x .or. &
          n%elements(4*j-3) >= n%x+n%width .or. &
          n%elements(4*j-2) < n%y .or. &
          n%elements(4*j-2) >= n%y+n%height &
        ) then
          !> if the particle is outside the cell, store it in list and remove from the cell
          particle(:) = n%elements(4*j-3:4*j)
          call append(list, particle)
          call removeElementFromNode(n, j)
        end if
      end do
    end do

    !> reinsert all particles from list in the now correct cells
    do i = 1, list%numParticles
      call getByIndex(list, i, particle)
      if (&
        particle(1) < tree%root%x .or. &
        particle(1) >= tree%root%x+tree%root%width .or. &
        particle(2) < tree%root%y .or. &
        particle(2) >= tree%root%y+tree%root%height &
      ) then
        cycle
      end if
      call insertElement(tree%root, particle)
      
    end do

    call deleteParticleList(list)
  end subroutine moveParticles
end module m_simulation
