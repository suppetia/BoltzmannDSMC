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
    nCollisions = .5_dp * cell%numberOfElements * cell%stats%particleCounter%average &
      * params%F_N * cell%stats%maxSigmaC * params%dt &
      / (params%V_c * cellWidth(cell) * cellHeight(cell)) 
    ! nCollisions = .5_dp * cell%numberOfElements * cell%numberOfElements &
    !   * params%F_N * cell%stats%maxSigmaC * params%dt &
    !   / (params%V_c * cellWidth(cell) * cellHeight(cell)) 

    if (nCollisions < 0) then
      return
    end if
    ! print *, "N_col", nCollisions
    ! print *, cell%numberOfElements
    ! print *, cell%stats%particleCounter%average
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
    ! print *, "N_col", nCollisions, nCurrentCollisions
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

    if (associated(pairs)) then
      deallocate(pairs)
    end if
  end subroutine collide

  
  subroutine moveParticles(dt, tree)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    real(dp), intent(in) :: dt

    type(NodeStack) :: stack
    type(QuadTreeNode), pointer :: n

    integer(i4) :: numLeaves, i, j, k, numCollisions

    type(ParticleList) :: list
    real(dp), dimension(5) :: particle
    real(dp) :: x, y, width, height

    ! integer, dimension(tree%params%maxElementsPerCell) :: collidedParticles

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
      
      ! collidedParticles(:) = -1
      call collideWithStructures(n, tree%params)!, collidedParticles)
      ! numCollisions = 0
      ! do j = 1,tree%params%maxElementsPerCell
      !   if (collidedParticles(j) == -1) then
      !     exit
      !   end if
      !   numCollisions = numCollisions + 1
      ! end do
      ! !> check if the particle collided with a structure
      ! !> in that case move it again in the new direction
      ! !> TODO: think whether this is fine
      ! !> idea: it's okay because the mean distance should be have the distance from the structure
      ! !> averaged over many particles
      ! do j = 1, numCollisions
      !   k = collidedParticles(j)
      !   print *, k
      !   n%elements(5*k-4) = n%elements(5*k-4) + dt * n%elements(5*k-2)
      !   n%elements(5*k-3) = n%elements(5*k-3) + dt * n%elements(5*k-1)
      ! end do

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

  subroutine collideWithStructures(cell, params)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: cell
    type(SimulationParameters), intent(in) :: params
    integer :: k

    integer :: i, j
    real(dp) :: v_x,v_y, s_x,s_y, s, t
    real(dp), dimension(:), pointer :: collisionDistance !> store the distance to the nearest structure

    k = 0

    if (.not.associated(cell%structures)) then
      return
    end if

    allocate(collisionDistance(size(cell%structures, 1)))

    do i = 1, cell%numberOfElements
      do j = 1, size(cell%structures, 1)
        !> find the distance to the structure as the projection on the normal vector
        ! print *, cell%elements(5*i-4:5*i)
        ! print *, cell%structures(j, :)
        collisionDistance(j) = (cell%elements(5*i-4) - cell%structures(j,1)) * cell%structures(j,5) &
          + (cell%elements(5*i-3)-cell%structures(j,2)) * cell%structures(j,6)
      end do

      if (all(collisionDistance <= 0)) then
        j = maxloc(collisionDistance, 1)
        v_x = cell%elements(5*i-2)
        v_y = cell%elements(5*i-1)
        !> update the velocity of the particle as v <- v - 2(v*n)n
        cell%elements(5*i-2) = v_x - 2*(v_x*cell%structures(j, 5) + v_y*cell%structures(j,6)) * cell%structures(j,5)
        cell%elements(5*i-1) = v_y - 2*(v_x*cell%structures(j, 5) + v_y*cell%structures(j,6)) * cell%structures(j,6)

        !> find the intersection of the structure and the particle trajectory
        s_x = cell%structures(j,3) - cell%structures(j,1)
        s_y = cell%structures(j,4) - cell%structures(j,2)
        ! !> denominator
        ! if (abs(d) < 1e-10) then
        !   t = (cell%elements(5*i-3)-cell%structures(j,2)-(cell%elements(5*i-4)-cell%structures(j,1))/s_x*s_y)/d
        ! else
        !   !> in this case the particle might be trapped
        !   t = params%dt * 1.0001_dp
        ! end if

        if (abs(s_x) < 1e-10_dp) then
          !> use a slight offset to move the particle into the right cell
          t = (cell%structures(j,1) - cell%elements(5*i-4))/v_x! * 1.00000001_dp
        else
          t = (cell%elements(5*i-3)-cell%structures(j,2)-(cell%elements(5*i-4)-cell%structures(j,1))/s_x*s_y)/(v_x*s_y/s_x-v_y)
        end if

        ! t = (cell%elements(5*i-3)-cell%structures(j,2)-(cell%elements(5*i-4)-cell%structures(j,1))/s_x*s_y)/(v_x*s_y/s_x-v_y)
        ! print *, "s", collisionDistance(j)
        ! print *, "t", t
        
        cell%elements(5*i-4) = cell%elements(5*i-4) + t*v_x - t*cell%elements(5*i-2)
        cell%elements(5*i-3) = cell%elements(5*i-3) + t*v_y - t*cell%elements(5*i-1)
        
        ! print *, "s_after", (cell%elements(5*i-4) - cell%structures(j,1)) * cell%structures(j,5) &
        !   + (cell%elements(5*i-3)-cell%structures(j,2)) * cell%structures(j,6)

      end if
    end do

    deallocate(collisionDistance)


    ! if (associated(cell%structures)) then
    !   do j = 1, cell%numberOfElements
    !     do i = 1, size(cell%structures, 1)
    !       !> check whether the particle collided with the structure using the cross-product
    !       ! s_x = (cell%elements(5*j-4) - cell%structures(i,1))
    !       ! s_y = (cell%elements(5*j-3) - cell%structures(i,2))
    !       ! print *, s_x, s_y
    !       ! s = sqrt(s_x**2+s_y**2)
    !       s = (cell%elements(5*j-4)-cell%structures(i,1)) * cell%structures(i,5) &
    !         + (cell%elements(5*j-3)-cell%structures(i,2)) * cell%structures(i,6)
    !       print *, "s", s
    !       if (s == 0._dp) then
    !         collisionDistance(i) = 0._dp
    !       else
    !         collisionDistance(i) = s
    !         ! collisionDistance(i) = -(s_x * (cell%structures(i,4)-cell%structures(i,2) &
    !         !   - s_y * (cell%structures(i,3)-cell%structures(i,1))))/s
    !
    !       end if
    !       ! print *, v_x,v_y,y_, cell%elements(5*j-2),cell%elements(5*j-1)
    !     end do
    !
    !     print *, collisionDistance(:)
    !     print *, collisionDistance <= 0
    !     if (all(collisionDistance <= 0)) then
    !       i = maxloc(collisionDistance, 1)
    !       print *, cell%elements(5*j-4:5*j)
    !       v_x = cell%elements(5*j-2)
    !       v_y = cell%elements(5*j-1)
    !       cell%elements(5*j-2) = v_x - 2*(cell%structures(i,5) * v_x + cell%structures(i,6)*v_y) * cell%structures(i,5)
    !       cell%elements(5*j-1) = v_y - 2*(cell%structures(i,5) * v_x + cell%structures(i,6)*v_y) * cell%structures(i,6)
    !
    !       !> calculate the intersection point of the motion of the particle and the structure
    !       !> based on https://paulbourke.net/geometry/pointlineplane/
    !       t = ((cell%structures(i, 3)-cell%structures(i,1))*(cell%elements(5*j-3)-cell%structures(i,2)) &
    !         - (cell%structures(i, 4)-cell%structures(i,2))*(cell%elements(5*j-4)-cell%structures(i,1))) &
    !         / ((cell%structures(i, 4)-cell%structures(i,2)) * v_x - (cell%structures(i, 3)-cell%structures(i,1)) * v_y)
    !       print *, "t", t
    !
    !       !> x_x = cell%elements(5*j-4), x_y = cell%elements(5*j-3)
    !       !> v_x = cell%elements(5*j-2), v_y = cell%elements(5*j-1)
    !       !> a_x = cell%structures(i,1), a_y = cell%structures(i,2)
    !       !> b_x = cell%structures(i,3), b_y = cell%structures(i,4)
    !       !> s_x = b_x-a_x, s_y = b_y - a_y
    !       s_x = cell%structures(i,3) - cell%structures(i,1)
    !       s_y = cell%structures(i,4) - cell%structures(i,2)
    !       t = (cell%elements(5*j-3)-cell%structures(i,2)-(cell%elements(5*j-4)-cell%structures(i,1))/s_x*s_y)/(v_x*s_y/s_x-v_y)
    !       print *, "t2", t
    !       t = abs(t)
    !
    !       if (t < 0) then
    !         cell%elements(5*j-4) = cell%elements(5*j-4) - t*v_x !+ (params%dt-t) * cell%elements(5*j-2)
    !         cell%elements(5*j-3) = cell%elements(5*j-3) - t*v_y !+ (params%dt-t) * cell%elements(5*j-1)
    !       end if
    !    
    !       print *, cell%elements(5*j-4:5*j)
    !       
    !       k = k+1
    !       ! collidedParticles(k) = j
    !     end if
    !   end do
    !
    !   deallocate(collisionDistance)
    ! end if
    

  end subroutine collideWithStructures

end module m_simulation
