module mod_types
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: dp, PI, I_C

  real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
  complex(dp), parameter :: I_C = (0.0_dp, 1.0_dp)
end module mod_types
