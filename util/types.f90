module m_types
  use iso_c_binding

  implicit none
  !> Char length for integers, range -2⁷ to 2⁷-1; 8 bits
  integer, parameter :: i1 = selected_int_kind(2)
  ! integer, parameter :: i1 = c_int8_t
  !> Short length for integers, range -2¹⁵ to 2¹⁵-1; 16 bits
  integer, parameter :: i2 = selected_int_kind(4)
  ! integer, parameter :: i2 = c_int16_t
  !> Length of default integers, range -2³¹ to 2³¹-1; 32 bits
  integer, parameter :: i4 = selected_int_kind(9)
  ! integer, parameter :: i4 = c_int32_t
  !> Long length for integers, range -2⁶³ to 2⁶³-1; 64 bits
  integer, parameter :: i8 = selected_int_kind(18)
  ! integer, parameter :: i8 = c_int64_t

  !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
  integer, parameter :: sp = selected_real_kind(6, 37)
  ! integer, parameter :: sp = c_float
  !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
  integer, parameter :: dp = selected_real_kind(15, 307)
  ! integer, parameter :: dp = c_double
  !> Quadruple precision real numbers, 33 digits, range 10⁻⁴⁹³¹ to 10⁴⁹³¹-1; 128 bits
  integer, parameter :: qp = selected_real_kind(33, 4931)
  ! integer, parameter :: qp = c_float128


  !> set the floating point precision
  integer, parameter :: fp = dp
  !> set the payload precision
  integer, parameter :: pp = i4


end module m_types

