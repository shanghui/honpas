! *** Module: nao2gto_data ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of NAO2GTO data
!!
!! This module defines global variables required by the NAO2GTO routines.
!! All of them should ideally be substituted by a design ensuring a
!! seamless data flow throughout the whole program.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2018 Created [Yann Pouillon]
! *****************************************************************************
module nao2gto_data

  use nao2gto_common
  use nao2gto_libint, only: Libderiv_t, Libint_t
  use nao2gto_types

  implicit none

  private

  integer, public :: hfx_call_counter = 0

  real(dp), public :: coeffs_kind_max0 = 0.0_dp
  real(dp), public :: log10_eps_schwarz = 0.0_dp

  type(eri_link_type), pointer, public :: eri_head => null()
  type(eri_link_type), pointer, public :: eri_ptr => null()
  type(eri_link_type), pointer, public :: eri_tail => null()

  real(dp), pointer, public :: eri_prescreen(:,:) => null()

  type(Libint_t)        , public :: hfx_libint
  type(Libderiv_t)      , public :: hfx_libderiv
  type(hfx_options_type), public :: hfx_options
  type(hfx_system_type) , public :: hfx_system

  type(pair_list_type), public :: list_ij, list_kl

  integer, pointer, public :: subshell(:) => null()

  type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer, public :: &
    pair_dist_radii_pgf => null(), sfc_pgf => null()
  type(hfx_screen_coeff_type), dimension(:,:,:,:), pointer, public :: &
    sfc_shell => null()
  type(hfx_screen_coeff_type), dimension(:,:), pointer, public :: &
    sfc_kind => null()

  logical, pointer, public :: um_cut(:,:) => null()

                    ! ------------------------------------ !

  ! CO: Cartesian GTOS, SO: Spherical GTOs
  !   1s+3p+6d     = 10
  !   1s+3p+6d+10f = 20
  !   ncosum(lmax)=sum_l (l+1)(l+2)/2 ; l=0,l_max
  !
  ! Consider magnetic quantum number m for RGTOs
  integer, save, target, public :: &
    nco(0:l_max), &
    nso(0:l_max), &
    ncosum(-1:l_max), &
    co(0:l_max,0:l_max,0:l_max), &
    coset(0:l_max,0:l_max,0:l_max), &
    indco(3,20)

                    ! ------------------------------------ !

  ! Information about GTOs, to use along species_info
  type(gto_info_type), target, allocatable, save, public :: hfx_gtos(:)

end module nao2gto_data
