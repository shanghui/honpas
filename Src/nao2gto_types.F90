! *** Module: nao2gto_types ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of NAO2GTO data
!!
!! This module defines data structures to handle NAO2GTO options and
!! prescreening tolerances of ERIs in HONPAS. The specifications are read
!! from the FDF input file.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2016 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_types

  use atm_types, only: maxn_orbnl, maxnorbs
  use nao2gto_common

  implicit none

  private

  !> \brief Data type to store a chained list of 4-center integral parameters
  type, public :: eri_link_type

    real(dp) :: gto_eri
    integer  :: four_index(4)
    type(eri_link_type), pointer :: next

  end type eri_link_type

  !> \brief Data type to store relevant cell parameters
  type, public :: hfx_cell_type
    integer :: nsc(3)
    real(dp) :: cell(3)
    real(dp) :: cell_r(3)
  end type hfx_cell_type


  !> \brief Data type to store Hartree-Fock exchange options that can be
  !!        read from SIESTA input files
  !!
  !! \par Default values from FDF:
  !!      - DM_trunc       = .true.  (use sparse DM to screen ERIs)
  !!      - dump_fit_data  = .true.  (dump fit data to nao2gto_fit.yml)
  !!      - farfield       = .true.  (far-near field screening)
  !!      - npts_fit       = NTBMAX  (number of data points to fit orbitals)
  !!      - potential_type = 1 (1/r/erfc(wr)/r ...)
  !!      - omega          = 0.11
  !!      - cutoff_radius  = ???
  !!      - eps_far        = 1.0e-6  (far-field screening tolerance)
  !!      - eps_pairlist   = 1.0e-6  (build shell pair-list tolerance)
  !!      - eps_schwarz    = 1.0e-6  (Schwarz tolerance)
  !!      - eps_stored     = 1.0e-6  (stored ERIs tolerance)
  type, public :: hfx_options_type

    logical   ::  DM_trunc = .false.
    logical   ::  dump_fit_data = .false.
    logical   ::  farfield = .false.
    integer   ::  npts_fit = -1
    integer   ::  potential_type = -1
    real(dp)  ::  omega = -1.0_dp
    real(dp)  ::  cutoff_radius = -1.0_dp
    real(dp)  ::  eps_farfield = -1.0_dp
    real(dp)  ::  eps_pairlist = -1.0_dp
    real(dp)  ::  eps_schwarz = -1.0_dp
    real(dp)  ::  eps_stored = -1.0_dp

  end type hfx_options_type

  !> \brief Data type to store information about orbital pairs
  type, public :: pair_list_element_type

    integer  :: pair(2)
    integer  :: nl_index
    real(dp) :: r1(3), r2(3)
    real(dp) :: dist2

  end type

  !> \brief Data type to store a list of orbital-pair information
  type, public :: pair_list_type

    type(pair_list_element_type), dimension(:), pointer :: element
    integer :: nelement

  end type pair_list_type

  !> Data type to store information about screening coefficients
  type, public :: hfx_screen_coeff_type

    real(dp) :: x(2)

  end type hfx_screen_coeff_type

  !> \brief Data type to point NAO2GTO routines to relevant SIESTA data
  type, public :: hfx_system_type

    real(dp) :: cell(3,3)
    real(dp) :: cell_r(3,3)

    integer , pointer :: maxnh
    integer , pointer :: na
    integer , pointer :: norb
    integer , pointer :: nspin
    integer , pointer :: nua
    integer , pointer :: nuo
    integer , pointer :: nuotot
    integer , pointer :: iaorb(:)
    integer , pointer :: indxua(:)
    integer , pointer :: iphorb(:)
    integer , pointer :: isa(:)
    integer , pointer :: listh(:)
    integer , pointer :: listhptr(:)
    integer , pointer :: nsc(:)
    integer , pointer :: numh(:)
    real(dp), pointer :: xa(:,:)

  end type hfx_system_type

                    ! ------------------------------------ !

  ! sphi: Slater-type orbital (spherical)
  ! cphi : Cartesian orbital
  type, public :: gto_info_type

    ! Added by Honghui Shang
    integer, dimension(maxn_orbnl)  ::  orbnl_contract
    integer, dimension(maxn_orbnl)  ::  orbnl_index_cphi
    integer, dimension(maxn_orbnl)  ::  orbnl_index_sphi
    real(dp) :: kind_radius
    real(dp),dimension(maxn_contract,maxn_orbnl) :: orbnl_zeta
    real(dp), dimension(maxn_contract,maxn_orbnl) :: orbnl_coefficient
    real(dp), dimension(maxn_contract,maxn_orbnl) :: pgf_radius
    real(dp),dimension(maxn_orbnl) ::   shell_radius

    ! Added by Xinming Qin
    ! zeta and coeff of the most diffuse GTO. Normalized adjoint GTO
    ! NAO = sum Cm*GTOm

    ! min zeta, the most diffuse gto         
    real(dp), dimension(maxn_orbnl) :: orbnl_adjoined_zeta

    ! the corresponding coefficient
    !real(dp), dimension(maxn_orbnl) :: orbnl_diff_coeff

    ! Added by Xinming Qin
    ! Max contraction coefficient: C=cartesian, R=spherical (radial)
    ! shell : same n and l,  
    ! CGTO (nco*k) ---> RGTO (nso*k) ---> PAO (nso)
    ! c2s(nco, nso)   Fit: Sum_k D_k*RGTO_k,  k = 0, M 
    ! RGTO (nso) = Transpose[c2s(nco, nso)] * CGTO (nco), matrix * vector 
    ! RGTO (io)   = sum c2s(1:nco, io) * CGTO(1:nco)
    !
    ! ERI (nsoa*nsob*nsoc*nsod),  nso PAO once !
    ! Contribution of each CGTO shell (nco) to PAO shell (nso) ERIs : 
    !
    !     Sum_k D_k* c2s_k(nco,nso)*GTO_k(nco) ), k=1, M
    !
    ! Actually, both D_k and normorlized coeff have been involved in c2s_k (nco, nso), 
    ! We calc the trans matrix sphi(nco*M, nso) !
    !
    ! So max contribution coeff of echo CGTO shell is:
    ! max( sum [sphi(1: nco, i)],  i = 1, nso ) !
    real(dp), dimension(maxn_orbnl) ::  orbnl_contraction_coeff

    !----------------cphi's data------------------------!
    integer :: norbs_cphi
    integer,  dimension(maxnorbs_cphi) :: orb_cphi_lx
    integer,  dimension(maxnorbs_cphi) :: orb_cphi_ly
    integer,  dimension(maxnorbs_cphi) :: orb_cphi_lz
    real(dp), dimension(maxnorbs_cphi) :: norm_cphi
    integer,  dimension(maxnorbs_cphi) :: orb_index_cphi
    integer,  dimension(maxnorbs_cphi) :: orb_n_cphi
    integer,  dimension(maxnorbs_cphi) :: orb_l_cphi
    real(dp), dimension(ncon_max,maxnorbs_cphi) :: cphi
    !---------------end cphi's data---------------------!

    !---------------sphi's data-------------------------!
    ! In SIESTA, all orbitals are actually spherical
    real(dp),dimension(ncon_max,maxnorbs) :: sphi
    !---------------end sphi's data---------------------!

  end type gto_info_type

end module nao2gto_types
