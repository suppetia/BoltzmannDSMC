module m_util
  use m_types, only: dp, i1, i2, i4
  implicit none

  real(dp) :: PI = 4._dp * datan(1._dp)

  type :: SimulationParameters
    real(dp) :: width, height !> dimensions
    real(dp) :: V_c !> cell volume
    real(dp) :: F_N !> number of real particles per simulated particles
    real(dp) :: dt  !> time step per iteration
    real(dp), pointer, dimension(:) :: d_ref !> diameter of the particles (index represents the particle type)
    real(dp), pointer, dimension(:) :: m !> masses of the particles
    integer(i1) :: collisionModel !> 1: hard sphere, 2: variable hard sphere (VHS), 3: variable soft sphere (VSS)
    integer(i2) :: cellHistoryLength

    integer(i4) :: maxElementsPerCell
  end type SimulationParameters

contains

  function sigma_T(particle1, particleType1, particle2, particleType2, params) result (sigma)
    implicit none
    real(dp), dimension(5), intent(in) :: particle1, particle2
    integer(i1), intent(in) :: particleType1, particleType2
    type(SimulationParameters), intent(in) :: params

    real(dp) :: sigma, d

    d = .5_dp * (params%d_ref(particleType1)+params%d_ref(particleType2))

    if (params%collisionModel == 1) then !> hard sphere
      sigma = PI * d**2
    else
      print *, "collision_model not implemented yet: ", params%collisionModel
      sigma = -1
    end if 

  end function sigma_T


  !> calculate the intersection of two lines, where the line is given by two points: [x1,y1,x2,y2]
  !> if there is no intersection or the lines are identical, return [-1,-1]
  function lineIntersection(line1, line2) result (intersection)
    implicit none
    real(dp), dimension(4), intent(in) :: line1, line2
    real(dp), dimension(2) :: intersection
    
    real(dp) :: d, u1, u2

    !> denominator
    d = (line2(4)-line2(2))*(line1(3)-line1(1)) - (line2(3)-line2(1))*(line1(4)-line1(2))

    !> if the lines are parallel
    if (abs(d) < 1e-10_dp) then
      intersection(:) = -1
      return
    end if

    u1 = ((line2(3)-line2(1)) * (line1(2)-line2(2)) - (line2(4)-line2(2))*(line1(1)-line2(1))) / d
    u2 = ((line1(3)-line1(1)) * (line1(2)-line2(2)) - (line1(4)-line1(2))*(line1(1)-line2(1))) / d

    ! print *, "u1, u2", u1, u2
    !> if the intersection is not in the line segment defined by the points for at least one of the lines
    if (u1 < 0 .or. u1 > 1 .or. u2 < 0 .or. u2 > 1) then
      intersection(:) = -1
      return
    end if

    intersection = line1(:2) + u1 * (line1(3:4)-line1(:2))
  end function lineIntersection


end module m_util
