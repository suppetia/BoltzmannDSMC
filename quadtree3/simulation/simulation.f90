module m_simulation  !
  use m_types, only: fp, i1, i4, i8
  use m_quadtree, only: QuadTreeNode, QuadTree, cellX, cellY, cellWidth, cellHeight, &
    insertParticles, findParticleCells, updateTreeNodes
  use m_datastructures, only: SimulationParameters, ParticleList, append, initializeParticleList, deleteParticleList
  use m_util, only: PI, sigma_T
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
  implicit none

contains

  subroutine step(tree, simParams)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), pointer, intent(in) :: simParams

    integer(i4) :: i, nCollisions, totalCollisions
    type(QuadTreeNode), pointer :: n
    integer(i4), dimension(:), pointer :: collisionPairs

    call moveParticles(tree, simParams)
    call updateTreeNodes(tree, simParams)
    totalCollisions = 0
    !$OMP PARALLEL DO private(n, nCollisions, collisionPairs) shared(tree, simParams) &
    !$OMP reduction(+: totalCollisions)
    do i = 1, tree%leafNumber
      n => tree%leafs(i)%node
      nullify(collisionPairs)
      call selectCollisionPairs(n, tree, simParams, collisionPairs, nCollisions)
      ! print *, "num col:", nCollisions, ASSOCIATED(collisionPairs)
      totalCollisions = totalCollisions + nCollisions
      call collide(n, tree, collisionPairs, nCollisions, simParams)
    end do
    !$OMP END PARALLEL DO
    print*, "number of collisions:", totalCollisions
  end subroutine step 

  subroutine moveParticles(tree, params)
    implicit none
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), intent(in) :: params

    type(QuadTreeNode), pointer :: n
    type(ParticleList) :: list
    integer(i4), dimension(:), pointer :: leafIdx, particleStartIdx, numParticlesThread
    integer(i4) :: i,j, numParticles, idx, threadID
    real(fp) :: x,y,width,height
    real(fp), dimension(:,:), pointer :: particles
    integer(i1), dimension(:), pointer :: particleTypes



    !$OMP PARALLEL private(list, n, numParticles, idx, x,y,width,height, threadID) &
    !$OMP shared(tree, particles, particleStartIdx, numParticlesThread)

    !$OMP SINGLE
    allocate(particleStartIdx(omp_get_num_threads()))
    allocate(numParticlesThread(omp_get_num_threads()))
    numParticlesThread = 0
    !$OMP END SINGLE

    threadID = omp_get_thread_num() + 1

    !> TODO: write a thread-safe version for this
    call initializeParticleList(list, chunksize=10000)


    !> TODO: can be parallelized if at the end the particle list is merged
    !$OMP DO
    do i = 1, tree%leafNumber
      n => tree%leafs(i)%node
      numParticles = tree%particleNumbers(i)
      idx = tree%particleStartIndices(i)
      !> update x and y position components
      tree%particles(idx:idx+numParticles-1, 1) = tree%particles(idx:idx+numParticles-1, 1)&
        + params%dt * tree%particles(idx:idx+numParticles-1, 3) 
      tree%particles(idx:idx+numParticles-1, 2) = tree%particles(idx:idx+numParticles-1, 2)&
        + params%dt * tree%particles(idx:idx+numParticles-1, 4) 

      call collideWithStructures(tree, i)

      x = cellX(n) * tree%treeParams%width
      y = cellY(n) * tree%treeParams%height
      width = cellWidth(n) * tree%treeParams%width
      height = cellWidth(n) * tree%treeParams%height


      do j = numParticles-1, 0, -1
      !> iterate backwards because if a particle is remove the last valid particle copied inplace
        !> TODO: add target cell hints?
        if (&
          tree%particles(idx+j, 1) < x .or. tree%particles(idx+j, 1) >= x+width .or. &
          tree%particles(idx+j, 2) < y .or. tree%particles(idx+j, 2) >= y+height &
        ) then
          call append(list, tree%particles(idx+j, :), tree%particleTypes(idx+j))
          numParticles = numParticles - 1
          tree%particles(idx+j, :) = tree%particles(idx+numParticles, :)
          tree%particleTypes(idx+j) = tree%particleTypes(idx+numParticles)
        end if 
      end do 
      tree%particleNumbers(i) = numParticles
    end do 
    !$OMP END DO
    numParticlesThread(threadID) = list%numParticles

    
    !> TODO: determine the actual cells using hints
    ! allocate(leafIdx(list%numParticles))

    !$OMP BARRIER
    !$OMP SINGLE
    particleStartIdx = 1
    do j = 1, omp_get_num_threads()-1
      particleStartIdx(j+1:) = particleStartIdx(j+1:) + numParticlesThread(j)
    end do 
    allocate(particles(sum(numParticlesThread),5))
    allocate(particleTypes(sum(numParticlesThread)))
    !$OMP END SINGLE
    
    particles(particleStartIdx(threadID):particleStartIdx(threadID)+numParticlesThread(threadID)-1,:) &
      = list%particles(:list%numParticles,:)
    particleTypes(particleStartIdx(threadID):particleStartIdx(threadID)+numParticlesThread(threadID)-1) &
      = list%particleTypes(:list%numParticles)

    call deleteParticleList(list)
    !$OMP END PARALLEL

    ! particles => list%particles(:list%numParticles, :)
    call findParticleCells(tree, particles, leafIdx)
    call insertParticles(tree, particles, particleTypes, leafIdx)
    deallocate(particles)
    deallocate(particleTypes)

    ! nullify(particles)
    ! deallocate(leafIdx)
    ! call deleteParticleList(list)

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
    particles => tree%particles(tree%particleStartIndices(nodeIdx):tree%particleStartIndices(nodeIdx)+numParticles-1, :)

    allocate(collisionDistance(numParticles, numStructures))

    do i = 1,numStructures
      collisionDistance(:, i) = (particles(:, 1) - structures(1,i)) * structures(5,i) &
        + (particles(:, 2) - structures(2,i)) * structures(6,i)
    end do 

    do i = 1, numParticles
      j = 0
      if (numStructures == 1 .and. collisionDistance(i,1) <= 0) then
        j = 1
      else if (numStructures == 2) then
        !> check if the inside area is concave or convex in this region
        !> convex => n_1 x n_2 < 0
        !> in this case both must be <= 0 for a collision
        ! if (all(collisionDistance(i, :) <= 0)) then
        !   j = maxloc(collisionDistance(i, :), 1) !> find index of closest structure
        ! end if
        if (structures(5, 1)*structures(6,2) - structures(6,1)*structures(5,2) < 0) then
          !> in this case only one dotp can be < 0 => this means collision
          if (any(collisionDistance(i,:) <= 0)) then
            j = minloc(collisionDistance(i,:), 1)
          end if 
        else
          !> in this case both must be <= 0 for a collision
          if (all(collisionDistance(i, :) <= 0)) then
            j = maxloc(collisionDistance(i, :), 1) !> find index of closest structure
          end if
        end if 
      else if (size(structures,2) > 2) then
        print *, "more than two structures per cell not implemented"
        return
      end if 
      if (j > 0) then
      ! if (all(collisionDistance(i, :) <= 0)) then
      !   j = maxloc(collisionDistance(i, :), 1) !> find index of closest structure
        v_x = particles(i, 3)
        v_y = particles(i, 4)
        !> update the velocity of the particle as v <- v - 2(v*n)n
        vdotn = 2*(v_x*structures(5,j) + v_y*structures(6,j))
        particles(i, 3) = v_x - vdotn * structures(5,j)
        particles(i, 4) = v_y - vdotn * structures(6,j)
        !> find the intersection of the structure and the particle trajectory
        s_x = structures(3,j) - structures(1,j)
        s_y = structures(4,j) - structures(2,j)
  
        if (abs(s_x) < 1e-10_fp) then
          t = (structures(1,j)-particles(i, 1))/v_x
        else
          t = (particles(i, 2) - structures(2,j) - (particles(i, 1)-structures(1,j))/s_x*s_y)/(v_x*s_y/s_x-v_y)
        end if 
  
        particles(i, 1) = particles(i, 1) + t*(v_x - particles(i, 3))
        particles(i, 2) = particles(i, 2) + t*(v_y - particles(i, 4))
      end if 
    end do 
    deallocate(collisionDistance)
    
  end subroutine collideWithStructures 

  subroutine selectCollisionPairs(node, tree, simParams, pairs, nCollisions)
    implicit none
    type(QuadTreeNode), pointer, intent(inout) :: node
    type(QuadTree), pointer, intent(inout) :: tree
    type(SimulationParameters), intent(inout) :: simParams
    integer(i4), dimension(:), pointer, intent(inout) :: pairs
    integer(i4), intent(inout) :: nCollisions

    integer(i4) :: nCurrentCollisions, numParticles
    integer(i4) :: p1, p2, i, idx

    real(fp), dimension(:,:), pointer :: particles
    integer(i1), dimension(:), pointer :: particleTypes
    real(fp), dimension(5) :: tmpParticle
    integer(i1) :: tmpParticleType
    real(fp) :: c_r, sigma
    real(fp) :: rnd

    idx = node%nodeIdx
    numParticles = tree%particleNumbers(idx)
    particles => tree%particles(tree%particleStartIndices(idx):tree%particleStartIndices(idx)+numParticles-1, :)
    particleTypes => tree%particleTypes(tree%particleStartIndices(idx):tree%particleStartIndices(idx)+numParticles-1)

    !> number of collisions in the cell = 1/2 * N N_avg F_N (sigma_T, c_r)_max dt/V_c
    ! print *, "hi"
    ! print *,.5_fp * numParticles * node%stats%particleCounter%average(1) 
    ! print *, simParams%F_N * node%stats%maxSigmaC
    ! print *, simParams%V_c * cellWidth(node) * cellHeight(node)
    ! print *, .5_fp * numParticles * node%stats%particleCounter%average(1) &
    !   * simParams%F_N * node%stats%maxSigmaC * simParams%dt &
    !   / (simParams%V_c * cellWidth(node) * cellHeight(node))
    nCollisions = .5_fp * numParticles * node%stats%particleCounter%average(1) &
      * simParams%F_N * node%stats%maxSigmaC * simParams%dt &
      / (simParams%V_c * cellWidth(node) * cellHeight(node))

    ! print *, nCollisions
    if (nCollisions <= 0) then
      nCollisions = 0
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

      tmpParticle = particles(p2, :)
      tmpParticleType = particleTypes(p2)
      particles(p2, :) = particles(p1, :)
      particleTypes(p2) = particleTypes(p1)
      particles(p1, :) = tmpParticle
      particleTypes(p1) = tmpParticleType

      do p2 = p1+1, numParticles
        c_r = norm2(abs(particles(p1, 3:)-particles(p2, 3:)))
        sigma = sigma_T(particles(p1, :), particleTypes(p1), particles(p2, :), particleTypes(p2), simParams)

        call random_number(rnd)
        !> the collision occurs with probability c_r*sigma_T/(c_r*sigma_T)_max
        if (rnd < c_r * sigma / node%stats%maxSigmaC) then

          if (p2 .ne. p1+1) then
            !> swap the particle p2 to position p1+1 to avoid that it is respected in other collisions
            tmpParticle = particles(p2, :)
            particles(p2, :) = particles(p1, :)
            particles(p1, :) = tmpParticle
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
    integer(i4), intent(inout) :: numCollisions
    type(SimulationParameters), intent(in) :: simParams

    integer(i4) :: p1, p2, i, idx, numParticles
    real(fp), dimension(:,:), pointer :: particles
    integer(i1), dimension(:), pointer :: particleTypes
    real(fp), dimension(3) :: vecC_r, vecNewC_r, vecC_m
    real(fp) :: c_r, m1, m2, tmp
    real(fp) :: c_chi, s_chi, chi, c_eps, s_eps, eps
    real(fp), dimension(2) :: rnd

    idx = node%nodeIdx
    numParticles = tree%particleNumbers(idx)
    particles => tree%particles(tree%particleStartIndices(idx):tree%particleStartIndices(idx)+numParticles-1, :)
    particleTypes => tree%particleTypes(tree%particleStartIndices(idx):tree%particleStartIndices(idx)+numParticles-1)

    do i = 1,numCollisions
      p1 = pairs(i)
      p2 = pairs(i) + 1
      vecC_r(:) = particles(p1, 3:) - particles(p2, 3:)
      c_r = norm2(vecC_r(:))
      m1 = simParams%m(particleTypes(p1)) !> TODO: look up the type of a particle
      m2 = simParams%m(particleTypes(p2))
      
      vecC_m(:) = (m1*particles(p1, 3:) + m2*particles(p2, 3:))/(m1+m2)

      !> use (V)HS logic (isotropic deflection)
      !> generate random deflection angles
      call random_number(rnd)
      chi = 2.*PI*rnd(1) !> random deflection angle
      c_chi = cos(chi)
      ! c_chi = 2.*rnd(1)-1. !> cosine of a random elevation angle
      s_chi = sqrt(1.-c_chi*c_chi) !> sine of that angle
      eps = 2.*PI*rnd(2) !> random azimuth angle
      c_eps = cos(eps)
      s_eps = sqrt(1.-c_eps*c_eps)

      !> (y_r^2+w_r^2)^{1/2}
      tmp = sqrt(vecC_r(2)*vecC_r(2)+vecC_r(3)*vecC_r(3))

      !> calculate the new relativ velocity components
      vecNewC_r(1) = vecC_r(1) * c_chi + s_chi * s_eps * tmp
      vecNewC_r(2) = vecC_r(2) * c_chi + s_chi * (c_r*vecC_r(3)*c_eps - vecC_r(1)*vecC_r(2)*s_eps) / tmp
      vecNewC_r(2) = vecC_r(3) * c_chi - s_chi * (c_r*vecC_r(2)*c_eps + vecC_r(1)*vecC_r(3)*s_eps) / tmp

      !> update the velocity components after the collision
      particles(p1, 3:) = vecC_m(:) + m2/(m1+m2) * vecNewC_r(:)
      particles(p2, 3:) = vecC_m(:) - m1/(m1+m2) * vecNewC_r(:)
    end do

    if (associated(pairs)) then
      deallocate(pairs)
    end if
  end subroutine collide

  
end module m_simulation 
