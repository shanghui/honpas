! *** Module: nao2gto_io ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief I/O module for Gaussian-based Hartree-Fock exchange
!!
!!  This module bridges the input file of SIESTA with the NAO2GTO routines,
!!  which calculate the Hartree-Fock exchange interaction using Gaussians.
!!
!! \note
!!      This file currently works with a version of Libint configured for
!!      LIBINT_MAX_AM=5 and LIBINT_MAX_AM1=4 only.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2010-2018 SIESTA Developers Group
!!
!! \par History
!!      - 11.2017 Reviewed for inclusion in SIESTA [Xinming Qin]
!!      - 01.2018 Brought together from separate files [Yann Pouillon]
! *****************************************************************************
module nao2gto_io

  implicit none

  private

  public :: &
    nao2gto_dump_system, &
    nao2gto_transfer

contains

  ! ***************************************************************************
  ! *** Public routines                                                     ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Displays a summary of a hfx_system_type data structure
  !!
  !! \param[in] hfx_sys: data structure to display
  ! ***************************************************************************
  subroutine nao2gto_dump_system(hfx_sys)

    use nao2gto_types, only: hfx_system_type

    implicit none

    ! Arguments
    type(hfx_system_type), intent(in) :: hfx_sys

    ! -------------------------------------------------------------------------

    if ( associated(hfx_sys%maxnh) ) then
      write(*, fmt='(2X,A,": ",I8)') "maxnh", hfx_sys%maxnh
    else
      write(*, fmt='(2X,A,": ",A8)') "maxnh", "null"
    end if
    if ( associated(hfx_sys%na) ) then
      write(*, fmt='(2X,A,": ",I8)') "na", hfx_sys%na
    else
      write(*, fmt='(2X,A,": ",A8)') "na", "null"
    end if
    if ( associated(hfx_sys%norb) ) then
      write(*, fmt='(2X,A,": ",I8)') "norb", hfx_sys%norb
    else
      write(*, fmt='(2X,A,": ",A8)') "norb", "null"
    end if
    if ( associated(hfx_sys%nspin) ) then
      write(*, fmt='(2X,A,": ",I8)') "nspin", hfx_sys%nspin
    else
      write(*, fmt='(2X,A,": ",A8)') "nspin", "null"
    end if
    if ( associated(hfx_sys%nua) ) then
      write(*, fmt='(2X,A,": ",I8)') "nua", hfx_sys%nua
    else
      write(*, fmt='(2X,A,": ",A8)') "nua", "null"
    end if
    if ( associated(hfx_sys%nuo) ) then
      write(*, fmt='(2X,A,": ",I8)') "nuo", hfx_sys%nuo
    else
      write(*, fmt='(2X,A,": ",A8)') "nuo", "null"
    end if
    if ( associated(hfx_sys%nuotot) ) then
      write(*, fmt='(2X,A,": ",I8)') "nuotot", hfx_sys%nuotot
    else
      write(*, fmt='(2X,A,": ",A8)') "nuotot", "null"
    end if
    if ( associated(hfx_sys%iaorb) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "iaorb", size(hfx_sys%iaorb, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "iaorb", "null"
    end if
    if ( associated(hfx_sys%indxua) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "indxua", size(hfx_sys%indxua, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "indxua", "null"
    end if
    if ( associated(hfx_sys%iphorb) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "iphorb", size(hfx_sys%iphorb, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "iphorb", "null"
    end if
    if ( associated(hfx_sys%isa) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "isa", size(hfx_sys%isa, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "isa", "null"
    end if
    if ( associated(hfx_sys%listh) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "listh", size(hfx_sys%listh, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "listh", "null"
    end if
    if ( associated(hfx_sys%listhptr) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "listhptr", size(hfx_sys%listhptr, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "listhptr", "null"
    end if
    if ( associated(hfx_sys%nsc) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "nsc", size(hfx_sys%nsc, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "nsc", "null"
    end if
    if ( associated(hfx_sys%numh) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "numh", size(hfx_sys%numh, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "numh", "null"
    end if
    if ( associated(hfx_sys%xa) ) then
      write(*, fmt='(2X,A,": array(",I8,",",I8")")') "xa", &
        size(hfx_sys%xa, 1), size(hfx_sys%xa, 2)
    else
      write(*, fmt='(2X,A,": ", A8)') "xa", "null"
    end if

  end subroutine nao2gto_dump_system

! *****************************************************************************
!> \brief Reads coefficients and builds the corresponding spherical Gaussians
!!
!! This routine reads the coefficients of Natural Atomic Orbitals (NAO) onto
!! Gaussian-Type Orbitals (GTO) from the FDF input of SIESTA and builds the
!! corresponding spherical Gaussians.
!!
!! \author Honghui Shang
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2010 Imported and connected to atm_transfer [Honghui Shang]
!!      - 12.2013 Modified for input from FDF format [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] l_max: ...
!! \param[in] co: ...
!! \param[out] orbtramat: ...
! *****************************************************************************
  subroutine nao2gto_transfer(gtos, hfx_opts)

    use units,     only: pi
    use alloc,     only: re_alloc, de_alloc
    use atm_types, only: maxnorbs, nspecies
    use atm_types, only: maxn_orbnl, species, species_info
    use atmfuncs,  only: lofio,mofio,labelfis
    use atomlist,  only: rmaxo
    use basis_specs, only: label2species
    use chemical,  only: species_label
    use m_io,      only: io_assign
    use radial,    only: rad_get
    use sys,       only: die
    use fdf
    use parsing
    use nao2gto_common
    !FIXME: restore after debugging
    !use nao2gto_data, only: co, coset, indco, nco, ncosum, nso
    use nao2gto_data
    use nao2gto_types, only: gto_info_type, hfx_options_type
    use nao2gto_utils, only: dfac, exp_radius, exp_radius_very_extended
    use nao2gto_transform, only: calc_c2s_matrix, cphi2sphi, orbtramat_type
    use nao2gto_wrappers, only: nao2gto_libint_dump, nao2gto_libderiv_dump
    use gaufre

    implicit none

    ! Arguments
    type(gto_info_type), intent(inout) :: gtos(nspecies)
    type(hfx_options_type), intent(out) :: hfx_opts

    ! Local variables
    character(len=132) :: msg
    character(len=GAUFRE_ERRMSG_LEN) :: errmsg
    logical :: do_fit, do_stop
    integer  :: errno, fit_fd, isp, io, iorbrd, inlz, l, m, n, jnlz, &
&     ipgf, irad, itry, jpgf, npts, lx, ly, lz, nco_max, nso_max, &
&     ico, ifdfblk, i_cphi, ren, rel, rezeta, recoeffs
    integer  :: num_cphi(maxnorbs)
    real(dp) :: expzet, fnorm, zeta, zetb, prefac, gcca, gccb, &
&     gauss_resid, dummy, fit_rad, fit_orb, fit_amin, fit_amax, fit_sig3, &
&     new_sig3
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    type(species_info), pointer :: spp => null()

    logical, dimension(:), pointer    :: chk_coeffs => null()
    integer, dimension(:,:), pointer :: hfx_contract => null()
    real(dp), dimension(:,:,:), pointer :: nao2gto_zeta => null()
    real(dp), dimension(:,:,:), pointer :: nao2gto_coefficient => null()
    real(dp), dimension(:), pointer   :: rad_pts => null()
    real(dp), dimension(:), pointer   :: orb_pts => null()
    type(orbtramat_type), dimension(:), pointer :: orbtramat => null()
    type(gaufre_autopilot_t) :: fit_pilot
    type(gaufre_orbital_t) :: orb_data

    ! -------------------------------------------------------------------------

    write(*,'(A,/)') "************************ Begin: HYBRID XC INITIALIZATION **********************"

    ! -------------------------------------------------------------------------
    !> Step 1: Initialization of nco, nso, co, coset, orbtramat, etc
    ! -------------------------------------------------------------------------
    ncosum(-1:l_max)=0
    do l=0,l_max
      nco(l) = (l + 1)*(l + 2)/2
      nso(l) = 2*l + 1
      ncosum(l) = ncosum(l-1) + nco(l)
    end do

    do lx=0,l_max
      do ly=0,l_max
        do lz=0,l_max
          l = lx + ly + lz
          if ( l > l_max ) cycle
          co(lx,ly,lz) = 1 + (l - lx)*(l - lx + 1)/2 + lz
          coset(lx,ly,lz) = ncosum(l-1) + co(lx,ly,lz)
        end do
      end do
    end do

    indco(:,:) = 0
    do l=0,l_max
      do lx=0,l
        do ly=0,l-lx
          lz = l - lx - ly
          indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
        end do
      end do
    end do

    ! One cannot use re_alloc with orbtramat, because it is a vector
    ! of structured types
    allocate(orbtramat(0:l_max))
    do l=0,l_max
      nco_max = nco(l)
      nso_max = nso(l)
      write(msg,'("orbtramat(",I1,")%c2s")') l
      nullify(orbtramat(l)%c2s)
      call re_alloc(orbtramat(l)%c2s, 1, nso_max, 1, nco_max, &
&       name=trim(msg), routine='nao2gto_transfer')
    enddo
    call calc_c2s_matrix(l_max, co, orbtramat)

    call re_alloc(hfx_contract, 1, maxn_orbnl, 1, nspecies, &
&     "nao2gto_transfer")
    call re_alloc(nao2gto_zeta, 1, maxn_contract, 1, maxn_orbnl, &
&     1, nspecies, "nao2gto_transfer")
    call re_alloc(nao2gto_coefficient, 1, maxn_contract, 1, maxn_orbnl, &
&     1, nspecies, "nao2gto_transfer")

    hfx_contract(:,:) = 0
    nao2gto_zeta(:,:,:) = 0.0_dp
    nao2gto_coefficient(:,:,:) = 0.0_dp

    call nao2gto_read_options(hfx_opts)

    ! -------------------------------------------------------------------------
    !> Step 2: Read Hartree-Fock exchange parameters from FDF input file
    !!
    !! \note See ldau_specs.f to understand how to read blocks
    ! -------------------------------------------------------------------------
    do_fit = .true.
    if (fdf_block('NAO2GTO', bfdf)) then

      do_fit = .false.

      ifdfblk = 1
      do while(fdf_bline(bfdf,pline))

        !> Step 2.a(nofit): Read species
        if ( .not. fdf_bmatch(pline,'ni') ) then
          write(msg,'(A," (line ",I4,")")') &
&           'Wrong format in NAO2GTO', ifdfblk
          call die(trim(msg))
        endif
        isp = label2species(fdf_bnames(pline,1))
        if (isp .eq. 0) then
          write(*,'(a,1x,a)') &
&           'WRONG species symbol in NAO2GTO:', &
&         trim(fdf_bnames(pline,1))
          call die()
        endif
        spp => species(isp)
        ifdfblk = ifdfblk + 1

        !> Step 2.b(nofit): Prepare variables that will tell whether
        !! all coefficients for all orbitals have been read
        call re_alloc(chk_coeffs, 1, spp%n_orbnl, "nao2gto_transfer")
        chk_coeffs(:) = .false.

        !> Step 2.c(nofit): Read data for each (n, l, zeta) triplet
        !! \note Coefficients can be provided in any order.
        do iorbrd=1,spp%n_orbnl

          !> Step 2.c.1(nofit): Read information about which orbital the
          !! coefficients correspond to
          if (.not. fdf_bline(bfdf, pline)) &
&           call die('Not enough information on the Gaussian expansion')
          if (fdf_bmatch(pline,'iiii')) then
            ren = fdf_bintegers(pline,1)
            rel = fdf_bintegers(pline,2)
            rezeta = fdf_bintegers(pline,3)
            recoeffs = fdf_bintegers(pline,4)
          else
            write(msg,'(A," (line ",I4,")")') &
&             'Wrong format in NAO2GTO', ifdfblk
            call die(trim(msg))
          endif

          !> Step 2.c.2(nofit): Translate n, l, zeta into an orbital index
          inlz = orb_nlz_to_index(spp, ren, rel, rezeta)
          if ( inlz == -1 ) then
            write(msg,'(A,3(1X,I2),A," (line ",I4,")")') &
&             'Could not find orbital(', ren, rel, rezeta, ') in SIESTA', &
&             ifdfblk
            call die(trim(msg))
          end if
          hfx_contract(inlz,isp) = recoeffs
          ifdfblk = ifdfblk + 1

          !> Step 2.c.3(nofit): Check for duplicates
          if ( .not. chk_coeffs(inlz) ) then
            chk_coeffs(inlz) = .true.
          else
            write(msg,'(3(A),3(1X,I2),A," (line ",I4,")")') &
&             'Duplicate coefficients for ', trim(spp%symbol), ' orbital(', &
&             ren, rel, rezeta,') in NAO2GTO block', ifdfblk
            call die(trim(msg))
          end if

          !> Step 2.c.3(nofit): Read Gaussian coefficients
          do jnlz=1,hfx_contract(inlz,isp)
            if (.not. fdf_bline(bfdf, pline)) &
&               call die('Not enough information on the Gaussian expansion')
            if (fdf_bmatch(pline,'vv')) then
              nao2gto_zeta(jnlz,inlz,isp) = fdf_bvalues(pline,1)
              nao2gto_coefficient(jnlz,inlz,isp) = fdf_bvalues(pline,2)
            else
              write(msg,'(A," (line ",I4,")")') &
&                 'Wrong format in NAO2GTO', ifdfblk
              call die(trim(msg))
            endif
            ifdfblk = ifdfblk + 1
          enddo

        enddo ! iorbrd=1,spp%n_orbnl

        !> Step 2.d(nofit): Check that all coefficients for the
        !! current species have actually been read
        do iorbrd=1,spp%n_orbnl
          if ( .not. chk_coeffs(iorbrd) ) then
            write(msg,'(A,3(1X,I2),A," (line ",I4,")")') &
&             'Missing coefficients for orbital(', ren, rel, rezeta, &
&             ') in NAO2GTO block', ifdfblk
            call die(trim(msg))
          end if
        end do
        call de_alloc(chk_coeffs, 'chk_coeffs', 'nao2gto_transfer')

      enddo ! while fdf_bline(bfdf,pline)

    else

      !> Step 2.a(fit): Prepare data stores
      call re_alloc(rad_pts, 1, hfx_opts%npts_fit, 'rad_pts', &
        'nao2gto_transfer')
      call re_alloc(orb_pts, 1, hfx_opts%npts_fit, 'orb_pts', &
        'nao2gto_transfer')

      do isp=1,nspecies

        spp => species(isp)
        spp%label = species_label(isp)

        do inlz=1,spp%n_orbnl

          !> Step 2.b(fit): Translate rad_func representation of the orbital
          !! \bug Hard-coded the array size!
          do irad=1,hfx_opts%npts_fit
            rad_pts(irad) = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&                   (real(hfx_opts%npts_fit, dp) - 1)
            call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)
          enddo
          orb_pts(hfx_opts%npts_fit) = 0.0_dp

          !> Step 2.c(fit): Perform Gaussian fitting and filter results
          call orb_data%init(rad_pts, orb_pts, gfstat=errno, gfmsg=errmsg)
          if ( errno /= GAUFRE_ERR_OK ) then
            write(*, fmt='("GAUFRE ERROR: ",A)') trim(errmsg)
          end if
          call orb_data%init_physics(trim(spp%label), spp%orbnl_n(inlz), &
&           spp%orbnl_l(inlz), spp%orbnl_z(inlz), spp%orbnl_ispol(inlz), &
&           spp%orbnl_pop(inlz), gfstat=errno, gfmsg=errmsg)
          if ( errno /= GAUFRE_ERR_OK ) then
            write(*, fmt='("GAUFRE ERROR: ",A)') trim(errmsg)
          end if
          call fit_pilot%init(method=1, &
&           gfstat=errno, gfmsg=errmsg)
          if ( errno /= GAUFRE_ERR_OK ) then
            write(*, fmt='("GAUFRE ERROR: ",A)') trim(errmsg)
          end if
          call fit_pilot%exec(orb_data, gfstat=errno, gfmsg=errmsg)
          ipgf = orb_data%get_ngfs()
          if ( errno == GAUFRE_ERR_OK ) then
            hfx_contract(inlz,isp) = ipgf
            do jnlz=1,ipgf
              nao2gto_zeta(jnlz,inlz,isp) = &
&               orb_data%trial(2*jnlz-1)
              nao2gto_coefficient(jnlz,inlz,isp) = &
&               orb_data%trial(2*jnlz)
            end do
          else
            write(*, fmt='("GAUFRE ERROR: ",A)') trim(errmsg)
          end if
          call fit_pilot%free()

        end do   ! inlz

      end do   ! isp

    endif   ! NAO2GTO block

    !> Step 2.d(fit): Clean-up the mess
    call de_alloc(rad_pts, 'rad_pts', 'nao2gto_transfer')
    call de_alloc(orb_pts, 'orb_pts', 'nao2gto_transfer')

    ! -------------------------------------------------------------------------
    !> Step 3: Process input data for each species
    ! -------------------------------------------------------------------------
    do isp=1,nspecies
      spp => species(isp)
      inlz = 0

      !> Step 3.a: Store GTOs for normal orbitals (m == -l)
      !!
      !! \bug Did not check the use of negative Z numbers
      !!      in original implementation
      do io=1,spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m .ne. -l ) cycle   ! not a normal orbital
        inlz = inlz + 1
        gtos(isp)%orbnl_contract(inlz) = hfx_contract(inlz,isp)
        do ipgf=1,gtos(isp)%orbnl_contract(inlz)
          gtos(isp)%orbnl_zeta(ipgf,inlz) = nao2gto_zeta(ipgf,inlz,isp)
          gtos(isp)%orbnl_coefficient(ipgf,inlz) = nao2gto_coefficient(ipgf,inlz,isp)
        enddo
        gtos(isp)%orbnl_adjoined_zeta(inlz) = &
&         minval(gtos(isp)%orbnl_zeta(1:gtos(isp)%orbnl_contract(inlz),inlz))
      enddo

      !> Step 3.b: Add cutoff radius for primitive GTOs and contracted
      !! GTOs, actually PAO. PAO itself has cutoff radius in siesta, we
      !! use the shell cutoff just for comparison
      !! (added by Xinming Qin, Oct. 2018).
      gtos(isp)%pgf_radius(1:maxn_contract,1:maxn_orbnl)=0.0d0
      gtos(isp)%shell_radius(1:maxn_orbnl)=0.0d0
      inlz = 0
      do io=1,spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m /= -l ) cycle   ! Not a normal orbital
        inlz = inlz+1
        do ipgf=1,gtos(isp)%orbnl_contract(inlz)
          gcca = nao2gto_coefficient(ipgf,inlz,isp)
          zeta = nao2gto_zeta(ipgf,inlz,isp)
          gtos(isp)%pgf_radius(ipgf,inlz)= exp_radius(l, zeta, 1.0E-5_dp, gcca)
        enddo
        gtos(isp)%shell_radius(inlz) = &
&         maxval(gtos(isp)%pgf_radius(1:gtos(isp)%orbnl_contract(inlz),inlz))
      enddo
      gtos(isp)%kind_radius = maxval(gtos(isp)%shell_radius(1:maxn_orbnl))

      !> Step 3.c: Compute the number of cphi coefficients for each orbital
      !! and store their total number in i_cphi
      num_cphi(:) = 0
      i_cphi = 1
      do io=1,spp%norbs
        if ( spp%orb_m(io) .lt. 2 ) then
          num_cphi(i_cphi) = io
          i_cphi = i_cphi + 1
        else if ( spp%orb_m(io) .eq. 2 ) then
          num_cphi(i_cphi) = io
          num_cphi(i_cphi+1) = io
          i_cphi = i_cphi + 2
        else if ( spp%orb_m(io) .eq. 3 ) then
          num_cphi(i_cphi) = io
          num_cphi(i_cphi+1) = io
          num_cphi(i_cphi+2) = io
          num_cphi(i_cphi+3) = io
          i_cphi = i_cphi + 4
        endif
      enddo
      i_cphi = i_cphi - 1

      !> Step 3.d: Initialize indices for cartesian coefficients
      !!
      !! \note Here, we only consider the lmax <= 3 case (s, p, d, and f)
      if ( spp%lmax_basis .lt. 2 ) then
        gtos(isp)%norbs_cphi = spp%norbs   !< Coefficients for s and p orbitals
      else
        gtos(isp)%norbs_cphi = i_cphi !< Coefficients for d and f orbitals
      endif

      do i_cphi=1,gtos(isp)%norbs_cphi
        io = num_cphi(i_cphi)
        gtos(isp)%orb_n_cphi(i_cphi) = spp%orb_n(io)
        gtos(isp)%orb_l_cphi(i_cphi) = spp%orb_l(io)
        gtos(isp)%orb_index_cphi(i_cphi) = spp%orb_index(io)
!       write(*,'(A,A,5(1X,A,"=",I4))') &
!         "[DEBUG][orb_index] ", trim(spp%label), &
!         "norbs_cphi", gtos(isp)%norbs_cphi, &
!         "io", io, &
!         "n", spp%orb_n(io), &
!         "l", spp%orb_l(io), &
!         "orb_index", spp%orb_index(io)
      enddo

      io = 0
      do inlz=1,spp%n_orbnl
        l = spp%orbnl_l(inlz)

        if ( inlz .eq. 1 ) then
          gtos(isp)%orbnl_index_cphi(inlz) = 1
          gtos(isp)%orbnl_index_sphi(inlz) = 1
        else
          gtos(isp)%orbnl_index_cphi(inlz) = gtos(isp)%orbnl_index_cphi(inlz-1) &
&                                    + nco(spp%orbnl_l(inlz-1))
          gtos(isp)%orbnl_index_sphi(inlz) = gtos(isp)%orbnl_index_sphi(inlz-1) &
&                                    + nso(spp%orbnl_l(inlz-1))
        endif

        do ico=ncosum(l-1)+1,ncosum(l)
          io = io + 1
          gtos(isp)%orb_cphi_lx(io) = indco(1,ico)
          gtos(isp)%orb_cphi_ly(io) = indco(2,ico)
          gtos(isp)%orb_cphi_lz(io) = indco(3,ico)
        enddo
      enddo

      !> Step 3.e: Compute norms of cartesian GTOs
      do io=1,gtos(isp)%norbs_cphi
        l = gtos(isp)%orb_l_cphi(io)
        expzet = 0.5_dp*REAL(2*l + 3,dp)
        fnorm = 0.0_dp
        inlz = gtos(isp)%orb_index_cphi(io)
        do ipgf=1,gtos(isp)%orbnl_contract(inlz)
          gcca = gtos(isp)%orbnl_coefficient(ipgf,inlz)
          zeta = gtos(isp)%orbnl_zeta(ipgf,inlz)
          do jpgf=1,gtos(isp)%orbnl_contract(inlz)
            gccb = gtos(isp)%orbnl_coefficient(jpgf,inlz)
            zetb = gtos(isp)%orbnl_zeta(jpgf,inlz)
            fnorm = fnorm + gcca*gccb/((zeta + zetb)**expzet)
          end do
        end do

        fnorm = (0.5_dp**l)*(pi**1.5_dp)*fnorm
        lx = gtos(isp)%orb_cphi_lx(io)
        ly = gtos(isp)%orb_cphi_ly(io)
        lz = gtos(isp)%orb_cphi_lz(io)
        prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
        gtos(isp)%norm_cphi(io) = 1.0_dp/SQRT(prefac*fnorm)
      enddo

      !> Step 3.f: Compute cphi and sphi
      gtos(isp)%cphi(1:ncon_max,1:gtos(isp)%norbs_cphi) = 0.0_dp
      gtos(isp)%sphi(1:ncon_max,1:spp%norbs) = 0.0_dp

      do io=1,gtos(isp)%norbs_cphi
        inlz = gtos(isp)%orb_index_cphi(io)
        lx = gtos(isp)%orb_cphi_lx(io)
        ly = gtos(isp)%orb_cphi_ly(io)
        lz = gtos(isp)%orb_cphi_lz(io)
        do jnlz=1,gtos(isp)%orbnl_contract(inlz)
          gtos(isp)%cphi((jnlz-1)*nco(gtos(isp)%orb_l_cphi(io))+co(lx,ly,lz),io) =  &
&           gtos(isp)%orbnl_coefficient(jnlz,inlz) * gtos(isp)%norm_cphi(io)
        enddo
      enddo

      call cphi2sphi(ncon_max, gtos(isp)%norbs_cphi, spp%norbs, l_max, &
&                    spp%n_orbnl, spp%orbnl_l, nco, nso, &
&                     gtos(isp)%orbnl_index_cphi, gtos(isp)%orbnl_index_sphi, &
&                     gtos(isp)%cphi, gtos(isp)%sphi, orbtramat)

      inlz = 0
      do io=1,spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m .ne. -l ) cycle     ! not a normal orbital
        inlz = inlz+1
        gtos(isp)%orbnl_contraction_coeff(inlz) =  &
&         maxval([(sum(abs(gtos(isp)%sphi(1:nco(l)*gtos(isp)%orbnl_contract(inlz),jnlz))), &
&           jnlz=io,io+nso(l)-1)])
      enddo

    enddo   ! is=1,nspecies

    ! -------------------------------------------------------------------------
    !> Step 4: Report about fitting parameters
    ! -------------------------------------------------------------------------
    call nao2gto_print_info(gtos, hfx_opts)

    if ( hfx_opts%dump_fit_data ) then

      call io_assign(fit_fd)
      open(unit=fit_fd, file="nao2gto_fit.yml", form="formatted", &
&       position="rewind", action="write", status="unknown")
      write(unit=fit_fd, fmt='(A/"---")') "%YAML 1.1"

      write(unit=fit_fd, fmt='(/A,":")') "hfx_options"
      write(unit=fit_fd, fmt='(2X,A,": ",I2)') "potential_type", &
&       hfx_opts%potential_type
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "omega", &
&       hfx_opts%omega
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_farfield", &
&       hfx_opts%eps_farfield
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_pairlist", &
&       hfx_opts%eps_pairlist
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_schwarz", &
&       hfx_opts%eps_schwarz
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_stored", &
&       hfx_opts%eps_stored
      if ( hfx_opts%DM_trunc ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "DM_trunc", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "DM_trunc", "False"
      end if
      if ( hfx_opts%dump_fit_data ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "dump_fit_data", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "dump_fit_data", "False"
      end if

      if ( hfx_opts%farfield ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "farfield", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "farfield", "False"
      end if

      do isp=1,nspecies
        spp => species(isp)
        spp%label = species_label(isp)

        do inlz=1,spp%n_orbnl
          write(unit=fit_fd, fmt='(/"---"//A,":")') "orbital"
          write(unit=fit_fd, fmt='(2X,A,": ",A)') "species", &
&           trim(adjustl(spp%label))
          write(unit=fit_fd, fmt='(2X,A,": ",I4)') "qn_n", spp%orbnl_n(inlz)
          write(unit=fit_fd, fmt='(2X,A,": ",I4)') "qn_l", spp%orbnl_l(inlz)
          write(unit=fit_fd, fmt='(2X,A,": ",I4)') "zeta", spp%orbnl_z(inlz)
          if ( spp%orbnl_ispol(inlz) ) then
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "polarized", "True"
          else
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "polarized", "False"
          end if
          write(unit=fit_fd, fmt='(2X,A,": ",E15.5)') &
&           "population", spp%orbnl_pop(inlz)

          if ( do_fit ) then
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "do_fit", "True"
          else
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "do_fit", "False"
          end if

          write(unit=fit_fd, fmt='(2X,A,":")') "raw_data"
          do irad=1,hfx_opts%npts_fit
            fit_rad = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&               real(hfx_opts%npts_fit, dp)
            call rad_get(spp%orbnl(inlz), fit_rad, fit_orb, dummy)
            write(unit=fit_fd, fmt='(4X,"- [",E20.8,", ",E20.8,"]")') &
&               fit_rad, fit_orb
          enddo

          write(unit=fit_fd, fmt='(2X,A,":")') "trial"
          do jnlz=1,hfx_contract(inlz,isp)
            write(unit=fit_fd, fmt='(4X,"- [",E24.8,", ",E24.8,"]")') &
&             gtos(isp)%orbnl_zeta(jnlz,inlz), &
              gtos(isp)%orbnl_coefficient(jnlz,inlz)
          end do
        end do
      end do

      write(unit=fit_fd, fmt='(/A)') "..."
      close(unit=fit_fd)

    end if   ! hfx_opts%dump_fit_data

    ! -------------------------------------------------------------------------
    !> Step 5: Check the completeness of the fitting
    ! -------------------------------------------------------------------------
    do_stop = .false.
    do isp=1,nspecies
      spp => species(isp)
      spp%label = species_label(isp)
      do inlz=1,spp%n_orbnl
        if ( hfx_contract(inlz,isp) == 0 ) then
          write(msg, fmt='(A," (",A,", n=",I2,", l=",I2,", zeta=",I2,")")') &
&           "Gaussian fitting of the orbitals failed for", &
&           trim(spp%label), spp%orbnl_n(inlz), spp%orbnl_l(inlz), spp%orbnl_z(inlz)
          write(*, fmt='(A)') msg
          do_stop = .true.
        end if
      end do
    end do

    ! -------------------------------------------------------------------------
    !> Step 6: Properly free memory
    ! -------------------------------------------------------------------------
    do l=0,l_max
      write(msg,'("orbtramat(",I1,")%c2s")') l
      call de_alloc(orbtramat(l)%c2s, name=trim(msg), &
&       routine='nao2gto_transfer')
      nullify(orbtramat(l)%c2s)
    enddo
    deallocate(orbtramat)
    nullify(orbtramat)

    call de_alloc(hfx_contract, name="hfx_contract", &
&     routine="nao2gto_transfer")
    call de_alloc(nao2gto_zeta, name="nao2gto_zeta", &
&     routine="nao2gto_transfer")
    call de_alloc(nao2gto_coefficient, name="nao2gto_coefficient", &
&     routine="nao2gto_transfer")

    ! -------------------------------------------------------------------------
    !> Step 7: Stop program if the fitting is incomplete
    ! -------------------------------------------------------------------------
    if ( do_stop ) then
      write(msg, fmt='(A)') &
&       "The fitting procedure was not successfully completed"
      call die(trim(msg))
    end if

    call nao2gto_libint_dump(hfx_libint)
    call nao2gto_libderiv_dump(hfx_libderiv)

    write(*,'(A,/)') "************************ End: HYBRID XC INITIALIZATION ************************"

  end subroutine nao2gto_transfer

  ! ***************************************************************************
  ! *** Private routines                                                    ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Displays the values of the specified Hartree-Fock exchange data
  !!        structure
  !!
  !! \param[in] gtos: data structure storing Gaussian fitting information
  !! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
  ! ***************************************************************************
  subroutine nao2gto_print_info(gtos, hfx_opts)

    use atm_types, only: nspecies, species, species_info
    use chemical,  only: species_label
    use parallel, only: IOnode
#ifdef BSC_CELLXC
    use bsc_xcmod, only: nXCfunc, XCfunc, XCauth
#else
    use siestaXC, only: getXC
#endif /* BSC_CELLXC */
    use nao2gto_common
    use nao2gto_types, only: gto_info_type, hfx_options_type

    implicit none

    ! Arguments
    type(gto_info_type), intent(in) :: gtos(nspecies)
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    integer :: inlz, isp, jnlz, nf
    type(species_info), pointer :: spp => null()
#ifndef BSC_CELLXC
    character(len=20) :: XCfunc(10), XCauth(10)
    integer :: nXCfunc
#endif /* BSC_CELLXC */

    ! -------------------------------------------------------------------------

    if ( .not. IOnode ) return


    write(*,'(/,A,/)') "nao2gto_print_info: NAO2GTO fitting information -------------------------------"

    ! Mark the beginning of the FDF data for automatic extraction
    write(*,'("# %%%",1X,A,/)') "HYBRID XC FDF BEGIN"

    ! Display XC parameters
#ifndef BSC_CELLXC
    call getXC(nXCfunc, XCfunc, XCauth)
#endif /* BSC_CELLXC */
    do nf=1,nXCfunc
      select case(hfx_opts%potential_type)
        case(do_hfx_potential_coulomb)
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "full Coulomb"
        case(do_hfx_potential_short)
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "short-range"
        case(do_hfx_potential_truncated)
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "truncated Coulomb"
        case default
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "*UNKNOWN*"
      end select
      write(*,'("XC.functional",3X,A)') trim(XCfunc(nf))
      write(*,'("XC.authors",6X,A)') trim(XCauth(nf))
    enddo

    ! Display the NAO2GTO block
    write(*,'(/,"#",1X,A,/,"#")') &
      "The NAO2GTO block is structured as following:"
    write(*,'("#",5X,A)') "species1 norbs_species1"
    write(*,'("#",5X,A)') "n1 l1 zeta1 ngaussians_orbital1"
    write(*,'("#",5X,A)') "exponent_gaussian1 coefficient_gaussian1"
    write(*,'("#",5X,A)') "exponent_gaussian2 coefficient_gaussian2"
    write(*,'("#",5X,A)') "exponent_gaussian3 coefficient_gaussian3"
    write(*,'("#",5X,A)') "..."
    write(*,'("#",5X,A)') "n2 l2 zeta2 ngaussians_orbital2"
    write(*,'("#",5X,A)') "exponent_gaussian1 coefficient_gaussian1"
    write(*,'("#",5X,A)') "..."
    write(*,'("#",5X,A)') "species2 norbs_species2"
    write(*,'("#",5X,A,/,"#")') "..."
    write(*,'(A)') "%block NAO2GTO"
    do isp=1,nspecies
      spp => species(isp)
      spp%label = species_label(isp)

      write(*,'(A,1X,I3)') trim(spp%label), spp%n_orbnl
      do inlz=1,spp%n_orbnl
        write(*,'(I1,3(1X,I2))') spp%orbnl_n(inlz), spp%orbnl_l(inlz), &
          spp%orbnl_z(inlz), gtos(isp)%orbnl_contract(inlz)
        do jnlz=1,gtos(isp)%orbnl_contract(inlz)
          write(*,'(E18.8,1X,E18.8)') gtos(isp)%orbnl_zeta(jnlz,inlz), &
            gtos(isp)%orbnl_coefficient(jnlz,inlz)
        enddo
      enddo
    enddo
    write(*,'(A)') "%endblock NAO2GTO"

    ! Display Hartree-Fock options
    write(*,'(/,"#",1X,A)') "Hartree-Fock exchange options"
    write(*,'(A,9X,L12)') "HFX.TruncateDM", hfx_opts%DM_trunc
    write(*,'(A,8X,L12)') "HFX.DumpFitData", hfx_opts%dump_fit_data
    write(*,'(A,11X,L12)') "HFX.FarField", hfx_opts%farfield
    write(*,'(A,6X,I12)') "HFX.FitDataPoints", hfx_opts%npts_fit
    write(*,'(A,2X,E12.3)') "HFX.FarFieldTolerance", hfx_opts%eps_farfield
    write(*,'(A,2X,E12.3)') "HFX.PairListTolerance", hfx_opts%eps_pairlist
    write(*,'(A,3X,E12.3)') "HFX.SchwarzTolerance", hfx_opts%eps_schwarz
    write(*,'(A,1X,E12.3)') "HFX.StoreERIsTolerance", hfx_opts%eps_stored
    write(*,'(A,14X,E12.3)') "HFX.Omega", hfx_opts%omega

    ! Mark the end of the FDF data for automatic extraction
    write(*,'(/,"# %%%",1X,A)') "HYBRID XC FDF END"

    write(*, '(/,a,/)') "nao2gto_print_info: END -------------------------------------------------------"

  end subroutine  nao2gto_print_info

! *****************************************************************************
!> \brief Reads Hartree-Fock exchange parameters from the FDF input
!!
!! This routine reads the values of the Hartree-Fock exchange parameters
!! from the FDF input of SIESTA and strores them into the specified data
!! structure. It also sets the default values for parameters not present
!! in the input file.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[out] hfx_opts: data structure storing Hartree-Fock exchange
!!                            parameters
! *****************************************************************************
  subroutine nao2gto_read_options(hfx_opts)

    use fdf, only: fdf_get
#ifdef BSC_CELLXC
    use bsc_xcmod, only: nXCfunc, XCfunc, XCauth
#else
    use siestaXC, only: getXC
#endif /* BSC_CELLXC */
    use atmparams, only: NTBMAX
    use nao2gto_common
    use nao2gto_types, only: hfx_options_type

    implicit none

    ! Arguments
    type(hfx_options_type), intent(out) :: hfx_opts

    ! Local variables
    integer :: nf
#ifndef BSC_CELLXC
    integer :: nXCfunc
    character(len=20) :: XCfunc(10), XCauth(10)
#endif /* BSC_CELLXC */

    ! -------------------------------------------------------------------------

#ifndef BSC_CELLXC
    call getXC(nXCfunc, XCfunc, XCauth)
#endif /* BSC_CELLXC */

    hfx_opts%DM_trunc = fdf_get("HFX.TruncateDM", .true.)
    hfx_opts%dump_fit_data = fdf_get("HFX.DumpFitData", .true.)
    hfx_opts%farfield = fdf_get("HFX.FarField", .true.)
    hfx_opts%npts_fit = fdf_get("HFX.FitDataPoints", NTBMAX)

    do nf = 1,nXCfunc
      if ( ((XCauth(nf).eq.'pbe0') .or. (XCauth(nf).eq.'PBE0')) &
&          .and. (XCfunc(nf).eq.'GGA') ) then
        hfx_opts%potential_type = do_hfx_potential_coulomb
      else
        if ( ((XCauth(nf).eq.'hse06') .or. (XCauth(nf).eq.'HSE06')) &
&            .and. (XCfunc(nf).eq.'GGA') ) then
          hfx_opts%potential_type = do_hfx_potential_short
        else
          hfx_opts%potential_type = do_hfx_potential_truncated
        endif
      endif
    enddo
    hfx_opts%omega = fdf_get("HFX.Omega", 0.11d0)

    hfx_opts%eps_farfield = fdf_get('HFX.FarFieldTolerance', 1.0d-6)
    hfx_opts%eps_pairlist = fdf_get('HFX.PairListTolerance', 1.0d-6)
    hfx_opts%eps_schwarz = fdf_get('HFX.SchwarzTolerance', 1.0d-6)
    hfx_opts%eps_stored = fdf_get('HFX.StoreERIsTolerance', 1.0d-6)

  end subroutine nao2gto_read_options

  ! ***************************************************************************
  !> \brief Finds the orbital index of a (n, l, zeta) triplet
  ! ***************************************************************************
  function orb_nlz_to_index(spp, orb_n, orb_l, orb_z) result(orb_index)

    use atm_types, only: species_info

    implicit none

    ! Arguments
    type(species_info), intent(in) :: spp
    integer, intent(in) :: orb_n
    integer, intent(in) :: orb_l
    integer, intent(in) :: orb_z

    ! Local variables
    integer :: orb_index
    integer :: iorb

    orb_index = -1

    do iorb=1,spp%n_orbnl
      if ( spp%orbnl_n(iorb) == orb_n ) then
        if ( spp%orbnl_l(iorb) == orb_l ) then
          if ( spp%orbnl_z(iorb) == orb_z ) then
            orb_index = iorb
            exit
          end if
        end if
      end if
    end do

  end function orb_nlz_to_index

end module nao2gto_io
