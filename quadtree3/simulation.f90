module m_simulation  !
  use m_types, only: fp, i1, i4, i8
  use m_quadtree, only: QuadTreeNode, QuadTree, cellX, cellY, cellWidth, cellHeight, &
    insertParticles, findParticleCells, updateTreeNodes
  use m_datastructures, only: SimulationParameters, ParticleList, append, initializeParticleList, deleteParticleList
  use m_util, only: PI, sigma_T
  implicit none

contains

  subroutine step(tree, simParams)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(in) :: simParams

    integer(i4) :: i, nCollisions
    type(QuadTreeNode), pointer :: n
    integer(i4), dimension(:), pointer :: collisionPairs

    call moveParticles(tree, simParams)
    call updateTreeNodes(tree, simParams)
    do i = 1, tree%leafNumber
      n => tree%leafs(i)%node
      call selectCollisionPairs(n, tree, simParams, collisionPairs, nCollisions)
      call collide(n, tree, collisionPairs, nCollisions, simParams)
    end do
  end subroutine step 

  subroutine moveParticles(tree, params)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), intent(in) :: params

    type(QuadTreeNode), pointer :: n
    type(ParticleList) :: list
    integer(i4), dimension(:), pointer :: leafIdx
    integer(i4) :: i,j, numParticles, idx
    real(fp) :: x,y,width,height
    real(fp), dimension(:,:), pointer :: particles

    integer :: counter

    !> TODO: write a thread-safe version for this
    call initializeParticleList(list, chunksize=1000)


    counter = 0
    !> TODO: can be parallelized if at the end the particle list is merged
    do i = 1, tree%leafNumber
      n => tree%leafs(i)%node
      numParticles = tree%particleNumbers(i)
      idx = tree%particleStartIndices(i)
      !> update x and y position components
      tree%particles(1,idx:idx+numParticles-1) = tree%particles(1,idx:idx+numParticles-1)&
        + params%dt * tree%particles(3,idx:idx+numParticles-1) 
      tree%particles(2,idx:idx+numParticles-1) = tree%particles(2,idx:idx+numParticles-1)&
        + params%dt * tree%particles(4,idx:idx+numParticles-1) 

      call collideWithStructures(tree, i)

      x = cellX(n) * tree%treeParams%width
      y = cellY(n) * tree%treeParams%height
      width = cellWidth(n) * tree%treeParams%width
      height = cellWidth(n) * tree%treeParams%height


      do j = numParticles-1, 0, -1
      !> iterate backwards because if a particle is remove the last valid particle copied inplace
        !> TODO: add target cell hints?
        if (&
          tree%particles(1, idx+j) < x .or. tree%particles(1, idx+j) >= x+width .or. &
          tree%particles(2, idx+j) < y .or. tree%particles(2, idx+j) >= y+height &
        ) then
          call append(list, tree%particles(:, idx+j))
          numParticles = numParticles - 1
          tree%particles(:, idx+j) = tree%particles(:, idx+numParticles)
          counter = counter + 1
        end if 
      end do 
      tree%particleNumbers(i) = numParticles
    end do 

    !> TODO: determine the actual cells using hints
    ! allocate(leafIdx(list%numParticles))

    particles => list%particles(:, :list%numParticles)
    call findParticleCells(tree, particles, leafIdx)
    call insertParticles(tree, particles, leafIdx)

    nullify(particles)
    ! deallocate(leafIdx)
    call deleteParticleList(list)

  end subroutine moveParticles 

  subroutine collideWithStructures(tree, nodeIdx)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    integer(i4), intent(in) :: nodeIdx

    real(fp), dimension(:,:), pointer :: collisionDistance, structures, particles
    integer(i4) :: i, numParticles
    integer(i1) :: numStructures, j
    real(fp) :: v_x, v_y, s_x, s_y, vdotn, t

    if (tree%structures(nodeIdx)%numStructures == 0) then
      return
    end if 

    structures => tree%structures(nodeIdx)%structures
    numStructures = tree%structures(nodeIdx)%numStructures

    numParticles = tree%particleNumbers(nodeIdx)
    particles => tree%particles(:, tree%particleStartIndices(nodeIdx):tree%particleStartIndices(nodeIdx)+numParticles-1)

    allocate(collisionDistance(numParticles, numStructures))

    do i = 1,numStructures
      collisionDistance(:, i) = (particles(1,:) - structures(1,i)) * structures(5,i) &
        + (particles(2,:) - structures(2,i)) * structures(6,i)
    end do 

    do i = 1, numParticles
      if (all(collisionDistance(i, :) <= 0)) then
        j = maxloc(collisionDistance(i, :), 1) !> find index of closest structure
        v_x = particles(3,i)
        v_y = particles(4,i)
        !> update the velocity of the particle as v <- v - 2(v*n)n
        vdotn = 2*(v_x*structures(5,j) + v_y*structures(6,j))
        particles(3,i) = v_x - vdotn * structures(5,j)
        particles(4,i) = v_y - vdotn * structures(6,j)
        !> find the intersection of the structure and the particle trajectory
        s_x = structures(3,j) - structures(1,j)
        s_y = structures(4,j) - structures(2,j)

        if (abs(s_x) < 1e-10_fp) then
          t = (structures(1,j)-particles(1,i))/v_x
        else
          t = (particles(2, i) - structures(2,j) - (particles(1,i)-structures(1,j))/s_x*s_y)/(v_x*s_y/s_x-v_y)
        end if 

        particles(1,i) = particles(1,i) + t*(v_x - particles(3,i))
        particles(2,i) = particles(2,i) + t*(v_y - particles(4,i))
      end if 
    end do 
    deallocate(collisionDistance)
    
  end subroutine collideWithStructures 

  subroutine selectCollisionPairs(node, tree, simParams, pairs, nCollisions)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), intent(inout) :: simParams
    integer(i4), dimension(:), pointer, intent(out) :: pairs
    integer(i4), intent(out) :: nCollisions

    integer(i4) :: nCurrentCollisions, numParticles
    integer(i4) :: p1, p2, i, idx

    real(fp), dimension(:,:), pointer :: particles
    real(fp), dimension(5) :: tmpParticle
    real(fp) :: c_r, sigma
    real(fp) :: rnd

    idx = node%nodeIdx
    numParticles = tree%particleNumbers(idx)
    particles => tree%particles(:, tree%particleStartIndices(idx):tree%particleStartIndices(idx)+numParticles-1)

    !> number of collisions in the cell = 1/2 * N N_avg F_N (sigma_T, c_r)_max dt/V_c
    nCollisions = .5_fp * numParticles * node%stats%particleCounter%average &
      * simParams%F_N * node%stats%maxSigmaC * simParams%dt &
      / (simParams%V_c * cellWidth(node) * cellHeight(node))

    if (nCollisions < 0) then
      return
    end if 

    allocate(pairs(nCollisions))

    nCurrentCollisions = 0

    !> start from the first particle to find a collision partner
    p1 = 1

    do while (nCurrentCollisions < nCollisions)
      !> to avoid statistical imbalance from always starting with the first particle
      !> swap a random particle to the first place
      call random_number(rnd)
      p2 = int(rnd*(numParticles-p1))+p1

      tmpParticle = particles(:,p2)
      particles(:, p2) = particles(:, p1)
      particles(:, p1) = tmpParticle

      do p2 = p1+1, numParticles
        c_r = norm2(abs(particles(3:, p1)-particles(3:, p2)))
        !> TODO: change to actual types
        sigma = sigma_T(particles(:, p1), 1_i1, particles(:, p2), 1_i1, simParams)

        call random_number(rnd)
        !> the collision occurs with probability c_r*sigma_T/(c_r*sigma_T)_max
        if (rnd < c_r * sigma / node%stats%maxSigmaC) then

          if (p2 .ne. p1+1) then
            !> swap the particle p2 to position p1+1 to avoid that it is respected in other collisions
            tmpParticle = particles(:,p2)
            particles(:, p2) = particles(:, p1)
            particles(:, p1) = tmpParticle
          end if
          !> store the index of the first collision particle
          !> the collision partner is at the next position
          pairs(nCurrentCollisions+1) = p1
          p1 = p1+1

          if (node%stats%maxSigmaC < c_r * sigma) then
            node%stats%maxSigmaC = c_r * sigma
          end if

          nCurrentCollisions = nCurrentCollisions + 1
          exit
        end if
      end do
      p1 = p1+1

      if (p1 >= numParticles-1) exit
    end do
    nCollisions = nCurrentCollisions
  end subroutine selectCollisionPairs

  subroutine collide(node, tree, pairs, numCollisions, simParams)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    integer(i4), dimension(:), pointer, intent(inout) :: pairs
    integer(i4), intent(in) :: numCollisions
    type(SimulationParameters), intent(in) :: simParams

    integer(i4) :: p1, p2, i, idx, numParticles
    real(fp), dimension(:,:), pointer :: particles
    real(fp), dimension(3) :: vecC_r, vecNewC_r, vecC_m
    real(fp) :: c_r, m1, m2
    real(fp) :: c_chi, s_chi, eps
    real(fp), dimension(2) :: rnd

    idx = node%nodeIdx
    numParticles = tree%particleNumbers(idx)
    particles => tree%particles(:, tree%particleStartIndices(idx):tree%particleStartIndices(idx)+numParticles-1)

    do i = 1,numCollisions
      p1 = pairs(i)
      p2 = pairs(i) + 1
      vecC_r(:) = particles(3:, p1) - particles(3:, p2)
      c_r = norm2(vecC_r(:))
      m1 = simParams%m(1) !> TODO: look up the type of a particle
      m2 = simParams%m(1)
      
      vecC_m(:) = (m1*particles(3:, p1) + m2*particles(3:, p2))/(m1+m2)

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
      particles(3:, p1) = vecC_m(:) + m2/(m1+m2) * vecNewC_r(:)
      particles(3:, p2) = vecC_m(:) - m1/(m1+m2) * vecNewC_r(:)
    end do

    if (associated(pairs)) then
      deallocate(pairs)
    end if
  end subroutine collide

  
end module m_simulation 
