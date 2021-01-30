! *** Module: nao2gto_common ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Common NAO2GTO constants and low-level routines
!!
!! This module stores constants and low-level routines which are used
!! by nearly all the NAO2GTO components of HONPAS.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2018 Created [Yann Pouillon]
! *****************************************************************************
module nao2gto_common

  implicit none

  private

                    ! ------------------------------------ !

  !
  ! Number kinds
  !

  !> Long integer kind, for Libint
  integer, parameter, public :: int_8 = selected_int_kind(10)

  !> Double-precision real kind
  integer, parameter, public :: dp = selected_real_kind(14, 200)

  !> Single-precision real kind
  integer, parameter, public :: sp = selected_real_kind(6, 30)

                    ! ------------------------------------ !

  !
  ! HFX potential types
  !

  !> Use the full Coulomb potential to compute Hartree-Fock XC contributions
  integer, parameter, public :: do_hfx_potential_coulomb = 1

  !> Use the short-range potential to compute Hartree-Fock XC contributions
  integer, parameter, public :: do_hfx_potential_short   = 2

  !> Truncate the Coulomb potential to compute Hartree-Fock XC contributions
  integer, parameter, public :: do_hfx_potential_truncated = 3

                    ! ------------------------------------ !

  !
  ! Boundaries
  !

  !> Limit to approximate log(0) when evaluating Hartree-Fock gradients
  real(dp), parameter, public :: log_zero = -1000.0_dp

  !> Minimum log value for the Powell optimization method
  real(dp), parameter, public :: powell_min_log = -20.0_dp

                    ! ------------------------------------ !

  !
  ! Gaussian-Type Orbitals
  !

  ! NAO2GTO: NAO = sum RGTOs = sum CGTOs
  ! RGTOs  : Real Spherical Harmonic GTOs or slater-type GTOs
  ! CGTOs  : Cartesian GTOs

  ! Number of CGTOs larger than maxnorbs (PAOs of an atom)
  integer, parameter, public :: maxnorbs_cphi = 120

  ! Max numbers of CGTOs and RGTOs
  integer, parameter, public :: maxn_contract = 10
  integer, parameter, public :: ncon_max = 60

  ! Support Max 6 GTOs for f orbitals, 10 GTOs for d orbitals
  integer, parameter, public :: l_max=3

end module nao2gto_common
