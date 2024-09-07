module m_simulation
  use m_quadtree, only: QuadTree, QuadTreeNode, removeElementFromNode, insertElement, &
    NodeStack, initializeStack, pop, cellHeight, cellWidth, cellX, cellY
  use m_quadtree_io, only: findLeaves
  use m_data_structures, only: ParticleList, initializeParticleList, append, getByIndex, deleteParticleList, addParticleCount
  use m_types, only: dp, i4, i2, sp, i1
  use m_util, only: SimulationParameters, sigma_T, PI

  implicit none

contains

  !> select collision pairs using the NTC method
  subroutine selectCollisionPairs(cell, params, pairs, nCollisions)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: cell
    type(SimulationParameters), intent(in) :: params
    integer(i4), intent(out) :: nCollisions

    integer(i2), dimension(:), pointer, intent(out) :: pairs

    integer(i4) :: nCurrentCollisions
    integer(i4) :: p1, p2, i

    real(dp), dimension(5) :: tmpParticle
    real(dp) :: c_r, sigma
    real(sp) :: rnd

    !> number of collisions in the cell = 1/2 * N N_avg F_N (sigma_T c_r)_max dt/V_c
    ! nCollisions = .5_dp * cell%numberOfElements * cell%stats%particleCounter%average &
    !   * params%F_N * cell%stats%maxSigmaC * params%dt &
    !   / (params%V_c * cellWidth(cell) * cellHeight(cell)) 
    nCollisions = .5_dp * cell%numberOfElements * cell%numberOfElements &
      * params%F_N * cell%stats%maxSigmaC * params%dt &
      / (params%V_c * cellWidth(cell) * cellHeight(cell)) 

    ! print *, nCollisions
    ! print *, params%F_N
    ! print *, cell%stats%maxSigmaC
    ! print *, params%dt
    ! print *,params%V_c * cellWidth(cell) * cellHeight(cell)
    ! if (nCollisions == 0) nCollisions = 1

    allocate(pairs(nCollisions))
    pairs(:) = 0
    nCurrentCollisions = 0

    !> start from the first particle to find a collision partner
    p1 = 1

    do while (nCurrentCollisions < nCollisions)
      !> to avoid statistical imbalance from always starting with the first particle
      !> swap a random particle to the first place
      call random_number(rnd)
      p2 = int(rnd*(cell%numberOfElements-p1))+p1

      tmpParticle(:) = cell%elements(5*p2-4:5*p2)
      cell%elements(5*p2-4:5*p2) = cell%elements(5*p1-4:5*p1)
      cell%elements(5*p1-4:5*p1) = tmpParticle(:)

      do p2 = p1+1, cell%numberOfElements
        c_r = norm2(abs(cell%elements(5*p1-2:5*p1)-cell%elements(5*p2-2:5*p2)))
        !> TODO: change to actual types
        sigma = sigma_T(cell%elements(5*p1-2:5*p1), 1_i1, cell%elements(5*p2-2:5*p2), 1_i1, params)

        call random_number(rnd)
        !> the collision occurs with probability c_r*sigma_T/(c_r*sigma_T)_max
        if (rnd < c_r * sigma / cell%stats%maxSigmaC) then

          if (p2 .ne. p1+1) then
            !> swap the particle p2 to position p1+1 to avoid that it is respected in other collisions
            tmpParticle(:) = cell%elements(5*p2-4:5*p2)
            cell%elements(5*p2-4:5*p2) = cell%elements(5*p1+1:5*p1+5)
            cell%elements(5*p1+1:5*p1+5) = tmpParticle(:)
          end if
          !> store the index of the first collision particle
          !> the collision partner is at the next position
          pairs(nCurrentCollisions+1) = p1
          ! print *, p1
          ! print *, "hi"
          p1 = p1+1

          if (cell%stats%maxSigmaC < c_r * sigma) then
            cell%stats%maxSigmaC = c_r * sigma
          end if

          nCurrentCollisions = nCurrentCollisions + 1
          exit
        end if
      end do
      p1 = p1+1

      if (p1 >= cell%numberOfElements-1) exit
    end do
    nCollisions = nCurrentCollisions

  end subroutine selectCollisionPairs

  subroutine collide(cell, pairs, numCollisions, params)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: cell
    integer(i2), dimension(:), pointer, intent(inout) :: pairs
    integer(i4), intent(in) :: numCollisions
    type(SimulationParameters), intent(in) :: params

    integer(i4) :: p1, p2, i
    real(dp), dimension(3) :: vecC_r, vecNewC_r, vecC_m
    real(dp) :: c_r, m1, m2
    real(dp) :: c_chi, s_chi, eps
    real(dp), dimension(2) :: rnd

    do i = 1,numCollisions
      p1 = pairs(i)
      p2 = pairs(i) + 1
      vecC_r(:) = cell%elements(5*p1-2:5*p1)-cell%elements(5*p2-2:5*p2)
      c_r = norm2(vecC_r(:))
      m1 = params%m(1) !> TODO: look up the type of a particle
      m2 = params%m(1)
      
      vecC_m(:) = (m1*cell%elements(5*p1-2:5*p1) + m2*cell%elements(5*p2-2:5*p2))/(m1+m2)

      !> use (V)HS logic (isotropic deflection)
      !> generate random deflection angles
      call random_number(rnd)
      c_chi = 2.*rnd(1)-1. !> cosine of a random elevation angle
      s_chi = sqrt(1.-c_chi*c_chi) !> sine of that angle
      eps = 2.*PI*rnd(2) !> random azimuth angle

      !> calculate the new relativ velocity components
      vecNewC_r(1) = c_r * c_chi
      vecNewC_r(2) = c_r * s_chi * cos(eps)
      vecNewC_r(3) = c_r * s_chi * sin(eps)

      !> update the velocity components after the collision
      cell%elements(5*p1-2:5*p1) = vecC_m(:) + m2/(m1+m2) * vecNewC_r(:)
      cell%elements(5*p2-2:5*p2) = vecC_m(:) - m1/(m1+m2) * vecNewC_r(:)

      ! print *, i
      ! print *, size(pairs, i)
      ! if (i == size(pairs, i)) exit

    end do

    deallocate(pairs)
  end subroutine collide

  
  subroutine moveParticles(dt, tree)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    real(dp), intent(in) :: dt

    type(NodeStack) :: stack
    type(QuadTreeNode), pointer :: n

    integer(i4) :: numLeaves, i, j

    type(ParticleList) :: list
    real(dp), dimension(5) :: particle
    real(dp) :: x, y, width, height

    call initializeStack(stack)
    call initializeParticleList(list)

    call findLeaves(tree%root, stack)

    numLeaves = stack%topIndex

    do i = 1, numLeaves
      call pop(stack, n)

      j = n%numberOfElements
      !> update x components: x_new = x + dt * v_x
      n%elements(1:5*j:5) = n%elements(1:5*j:5) + dt * n%elements(3:5*j:5)
      !> update y components: y_new = y + dt*v_y
      n%elements(2:5*j:5) = n%elements(2:5*j:5) + dt * n%elements(4:5*j:5)

      !> check if the particle is still in the same cell after the location update
      x = cellX(n) * tree%width
      y = cellY(n) * tree%height
      width = cellWidth(n) * tree%width
      height = cellHeight(n) * tree%height
      do j = n%numberOfElements, 1, -1
        !> iterate backwards because of how the particles are removed:
        !> if the particle at index i is removed, it is replaced with the last particle
        if ( &
          n%elements(5*j-4) < x .or. &
          n%elements(5*j-4) >= x+width .or. &
          n%elements(5*j-3) < y .or. &
          n%elements(5*j-3) >= y+height &
        ) then
          !> if the particle is outside the cell, store it in list and remove from the cell
          particle(:) = n%elements(5*j-4:5*j)
          call append(list, particle)
          call removeElementFromNode(n, j)
        end if
      end do
    end do

    !> reinsert all particles from list in the now correct cells
    x = cellX(tree%root) * tree%width
    y = cellY(tree%root) * tree%height
    width = cellWidth(tree%root) * tree%width
    height = cellHeight(tree%root) * tree%height

    do i = 1, list%numParticles
      call getByIndex(list, i, particle)
      if (&
        particle(1) < x .or. &
        particle(1) >= x+width .or. &
        particle(2) < y .or. &
        particle(2) >= y+height &
      ) then
        cycle
      end if
      call insertElement(tree%root, particle, tree)
    end do

    call deleteParticleList(list)
  end subroutine moveParticles

end module m_simulation
