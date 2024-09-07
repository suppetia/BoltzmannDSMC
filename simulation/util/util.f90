module m_util
  use m_types, only: dp, i1, i2
  implicit none

  real(dp) :: PI = 4._dp * datan(1._dp)

  type :: SimulationParameters
    real(dp) :: V_c !> cell volume
    real(dp) :: F_N !> number of real particles per simulated particles
    real(dp) :: dt  !> time step per iteration
    real(dp), pointer, dimension(:) :: d_ref !> diameter of the particles (index represents the particle type)
    real(dp), pointer, dimension(:) :: m !> masses of the particles
    integer(i1) :: collisionModel !> 1: hard sphere, 2: variable hard sphere (VHS), 3: variable soft sphere (VSS)
    integer(i2) :: cellHistoryLength
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



end module m_util
