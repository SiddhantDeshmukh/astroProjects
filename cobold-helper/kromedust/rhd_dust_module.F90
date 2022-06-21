!***************************************************************************************************
!
!   #####   #    #  #####         #####   #    #   ####  #####
!   #    #  #    #  #    #        #    #  #    #  #        #
!   #    #  ######  #    #        #    #  #    #   ####    #
!   #####   #    #  #    #        #    #  #    #       #   #
!   #   #   #    #  #    #        #    #  #    #  #    #   #
!   #    #  #    #  #####  ###### #####    ####    ####    # _module
!
!   Radiation-Hydrodynamics-Code: dust module
!
!***************************************************************************************************
!   This module provides the master source-term routine that calls the individual routines
!   to handle the source terms to treat dust formation, molecular networks, or hydrogen ionization.
!***************************************************************************************************
!   Fortran90
!   2020-01-24
!***************************************************************************************************
! With contributions from:
!   Bernd Freytag:         Uppsala, Lyon
!   Matthias Steffen:      Potsdam
!   Susanne Hoefner:       Uppsala
!   Sven Wedemeyer-Boehm:  Freiburg
!   Hans-Guenter Ludwig:   Kiel, Lund, Meudon, Heidelberg
!***************************************************************************************************
! MODULES:
!   rhd_dust_module:       Routines for handling of source terms due to dust or molecules
!
!***************************************************************************************************

!# /* --- Default values for some preprocessor macros ------------------------------------------- */
!#
!# /* --- Allocation of local 3D arrays                                                       --- */
!# /*     0: default: declaration in header of subroutine                                         */
!# /*     1: pointer arrays                                                                       */
!# /*     2: allocatable arrays                                                                   */
# ifndef macro_arrays_LocalVarType_Case
#   define macro_arrays_LocalVarType_Case 0
# endif
!#
!# /* ------------------------------------------------------------------------------------------- */


!------***************------------------------------------------------------------------------------
module rhd_dust_module
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_module ('rhd_dust_module')
!
! PURPOSE:
!   Provide routines for dust handling
!
! CATEGORY:
!   Hydrodynamics, Dust
!
! CALLING SEQUENCE:
!   use rhd_dust_module
!
! VARIABLES:
!
! ROUTINES: (contained)
!   rhd_dust_SwitchInfo:   Print information about used compiler switches
!   rhd_dust_simpleSourceStep: Update model by applying simple (test) dust source term step
!   rhd_dust_COSourceStep: Update model by applying CO source term step
!   rhd_dust_K3MONSourceStep: Update model by applying dust source term step (one moment only)
!   rhd_dust_K3MON2SourceStep: Update model by applying dust source term step (one moment only)
!   rhd_dust_M4C2SourceStep: Update model by applying dust source term step (4 moments)
!   rhd_dust_Mg2SiO4AGBSourceStep: Update model by applying forsterite source term step
!   rhd_dust_Mg2SiO4pAl2O3AGBSourceStep: Update model by applying forsterite+corundum source terms
!
! MODIFICATION HISTORY:
!   2002-10-29 (B.F.) First version
!   2002-11-20 (B.F.) More source step subroutines
!   2002-12-17 (B.F. & S.H.) First version of M4C2
!   2005-01-25 (S.W.-B.) case 'chemreacnet'
!   2006-10-11 (H.G.L.)  dustform_4 (carbon) replaced by dustform_5 (fosterite)
!   2008-08-12 (M.S., B.F.) eosinter_T gets additional arguments (m1,n1,...)
!   2012-02-13 (B.F.) eosinter_T -> eosinter_TOMP
!   2012-03-08 (B.F.) Add routine rhd_dust_Mg2SiO4AGBSourceStep
!   2017-05-22 (B.F. Uppsala) Add routine rhd_dust_Mg2SiO4pAl2O3AGBSourceStep
!   2017-09-17 (B.F.) eosinter_TOMP gets aditional arguments (nghost1,...)
!   2019-11-11 (B.F.) Add routine rhd_dust_SwitchInfo
!   2019-11-24 (B.F.) Move routines to rhd_source_module and rename them
!---------------------------------------------------------------------------------------------------
implicit none
!
contains


!----------*******************----------------------------------------------------------------------
subroutine rhd_dust_SwitchInfo(ncp)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_SwitchInfo ('rhd_dust_Switch_Information')
!
! PURPOSE:
!   Print information about used compiler switches.
!
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_SwitchInfo(ncp)
!
! INPUT:
!   ncp:          ('number_channel_print') integer with output channel number (usually 6)
!
! OUTPUT:
!   None
!
! LOCAL VARIABLES:
!
! ROUTINES:
!   None
!
! MODULES:
!   None
!
! SIDE EFFECTS:
!   None
!
! RESTRICTIONS:
!   If another compiler switch is used anywhere in the module it has to be explicitly
!   included in this routine.
!
! PROCEDURE:
!   Just print the value of all compiler switches used for the module.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2011-11-11 (B.F.) Adapted from rhd_shortrad_SwitchInfo
!   2019-11-24 (B.F. Uppsala) Remove some compiler switches
!---------------------------------------------------------------------------------------------------
implicit none
!
! --- I/O parameters ---
integer,                   intent(in)   :: ncp
!
! --- Local variables ---
character(len=*), parameter             :: defstr='     defined', undefstr='   undefined'
!---------------------------------------------------------------------------------------------------
!
write(ncp,'(A)')    'Compiler switches: rhd_dust_module ........................'
!
write(ncp,'(A,I6)') '  macro_'//'arrays_LocalVarType_Case:     ',macro_arrays_LocalVarType_Case
!
!#
# ifdef rhd_dust_t01
!#
write(ncp,'(A,A)' ) '  rhd_'//'dust_t01:                       ',defstr
!#
# else
!#
write(ncp,'(A,A)' ) '  rhd_'//'dust_t01:                       ',undefstr
!#
# endif
!#
!
end subroutine rhd_dust_SwitchInfo


!----------*************************----------------------------------------------------------------
subroutine rhd_dust_simpleSourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_simpleSourceStep ('rhd_dust_simple_Source_Step')
!
! PURPOSE:
!   Update model by applying simple (test) dust source term step.
!   
! CATEGORY:
!   Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_simpleSourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                      e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('n_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('n_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('n_3') integer, upper index, 3rd dimension
!
! ROUTINES:
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-10-29 (B.F.) First version
!   2002-11-20 (B.F.) Put into own routine
!   2009-10-15 (B.F.) Parallelize call to eosinter
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3, iquc
real                                    :: factor
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
integer                                 :: ierr0
character(len=80)                       :: outstr0
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
!
! --- Solve EOS, find out Helium abundance (number fraction) ---
!
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i1, i2, iquc, factor)
!
select case(par%dustScheme)
  case('dust_simple_01')
    ! --- Simple test 'dust' formation (one component: nquc=1) ---
    !
    iquc=1
    !
    !$OMP DO
    do i3=m3,n3
      do i2=m2,n2
        do i1=m1,n1
          factor=exp((2800.0-T(i1,i2,i3))/1000.0*dtime/20.0)
          model%quc(i1,i2,i3,iquc)=max(min( ( model%quc(i1,i2,i3,iquc) + &
                                              0.1*model%rho(i1,i2,i3) * max(factor-1.0,0.0) &
                                            ) * factor, &
                                            model%rho(i1,i2,i3) &
                                          ), &
                                       0.0 &
                                      )
        end do ! i1
      end do ! i2
    end do ! i3
    !$OMP END DO
    !
end select ! par%dustScheme
!
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
end subroutine rhd_dust_simpleSourceStep


!----------*********************--------------------------------------------------------------------
subroutine rhd_dust_COSourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_COSourceStep ('rhd_dust_CO_Source_Step')
!
! PURPOSE:
!   Update model by applying CO source term step.
!   
! CATEGORY:
!   Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_COSourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                      e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('n_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('n_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('n_3') integer, upper index, 3rd dimension
!
! ROUTINES:
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-10-29 (B.F.) First version
!   2002-11-04 (B.F., M.S.) First CO one-component model
!   2002-11-20 (B.F.) Put into own routine
!   2009-10-15 (B.F.) Parallelize call to eosinter
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3, iquc
real                                    :: T_red, T_red_06, T_red_22, nH15_red, k1, k2
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
integer                                 :: ierr0
character(len=80)                       :: outstr0
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
!
! --- Solve EOS, find out Helium abundance (number fraction) ---
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
! --- Typical loop ---
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i1, i2, iquc, T_red, T_red_06, T_red_22, nH15_red, k1, k2) 
!
select case(par%dustScheme)
  case('co_component01_01')
    ! --- Simple CO formation (one component: nquc=1) ---
    !
    iquc=1
    !
    !$OMP DO
    do i3=m3,n3
      do i2=m2,n2
        do i1=m1,n1
          T_red   =2.0E-04   * min(max(T(i1,i2,i3), 10.0), 1.0E+05)
          T_red_06=T_red**0.6
          T_red_22=T_red**22.2
          nH15_red=4.215e+08 * model%rho(i1,i2,i3)
          k1      =1.0E+07   * nH15_red**2 * T_red_06
          k2      =2.5E-05   * nH15_red    * T_red_06 * (1 + 40.0 * T_red_22)
          !
          model%quc(i1,i2,i3,iquc)=k1/k2 + (model%quc(i1,i2,i3,iquc) - k1/k2) * exp(-k2*dtime)
        end do ! i1
      end do ! i2
    end do ! i3
    !$OMP END DO
    !
end select ! par%dustScheme
!
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
end subroutine rhd_dust_COSourceStep


!----------************************-----------------------------------------------------------------
subroutine rhd_dust_K3MONSourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_K3MONSourceStep ('rhd_dust_K3_monomer_Source_Step')
!
! PURPOSE:
!   Update model by applying dust source term step (one moment only).
!   
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_K3MONSourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                      e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('n_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('n_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('n_3') integer, upper index, 3rd dimension
!
! ROUTINES:
!   DUSTFORM_3:            Compute rates for simplified forsterite dust formation model
!   DUSTFORM_6:            Compute condensation or evaporation rate for forsterite
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   tabinter_getinfo:      Give requested information about EOS interpolation table
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   dust_k3mon_module:     Simplified forsterite dust formation model (1 moment type)
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Control parameters:
!     par%C_dust01: CzuO
!     par%C_dust02: xnDmax
!     par%C_dust03: Tredfactor
!    (par%C_dust04: Settling velocity factor)
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-11-20 (B.F.) First version
!   2002-11-22 (B.F.) Because of instabilities due to short time scales
!                     approximate linear equation introduced and solved exactly
!   2003-08-27 (B.F.) Change of meaning of control parameters: par%C_dust01: CzuO
!                     par%C_dust01 -> par%C_dust02, 02 -> 03, 03 -> 04
!   2007-12-03 (B.F. Lyon) Call dustform_6 instead of dustform_5
!   2009-10-15 (B.F.) Parallelize call to eosinter
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
use dust_k3mon_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3
integer                                 :: npt
real                                    :: abuY, CzuO, xnDmax, Tredfactor, r1DUST
real, dimension(model%m1:model%n1)      :: taui, DUSTK0, xnCtot, xnCgas, rateK3, equK3
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
integer                                 :: ierr0
character(len=80)                       :: outstr0
integer, parameter                      :: i_DUSTK3=1, i_nCgas=2
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
npt=n1-m1+1
!
! --- Solve EOS, find out Helium abundance (number fraction) ---
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
call tabinter_getinfo(abuY=abuY)
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
!!$CzuO      =1.8
!!$xnDmax    =par%C_dust01
!!$Tredfactor=par%C_dust02
CzuO      =par%C_dust01
xnDmax    =par%C_dust02
Tredfactor=par%C_dust03
!
! --- Typical loop ---
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i2, taui, DUSTK0, xnCtot, xnCgas, r1DUST, rateK3, equK3) 
!
select case(par%dustScheme)
  case('dust_k3mon_01')
    ! --- Simple dust formation (one moment (nquc=1): K3) ---
    !
    !$OMP DO
    do i3=m3,n3
      do i2=m2,n2
        ! --- Compute dust formation rate ---
        call DUSTFORM_3(npt, &
                        model%rho(m1:n1,i2,i3), Tredfactor*T(m1:n1,i2,i3), &
                        model%quc(m1:n1,i2,i3,i_DUSTK3), &
                        CzuO, xnDmax,  &
                        taui, DUSTK0, xnCtot, xnCgas, r1DUST, &
                        rateK3, equK3)
!!$if (i3 == 16) then
!!$  write(6,'(A,11(1PE11.3))') 'rateK3: ', rateK3(20:30)
!!$  write(6,'(A,11(1PE11.3))') 'equK3:  ', equK3(20:30)
!!$endif
        !
        ! --- Prevent existence of dust grains with "impossible" sizes ---
        model%quc(m1:n1,i2,i3,i_DUSTK3)=min(max(model%quc(m1:n1,i2,i3,i_DUSTK3), &
                                                DUSTK0), &
                                            xnCtot)
        !
        ! --- Add source term ---
!!$        model%quc(m1:n1,i2,i3,i_DUSTK3)=model%quc(m1:n1,i2,i3,i_DUSTK3) + dtime * &
!!$                                          (taui * model%quc(m1:n1,i2,i3,i_DUSTK3)**(2.0/3.0)  &
!!$                                                * DUSTK0**(1.0/3.0))
        model%quc(m1:n1,i2,i3,i_DUSTK3)=equK3 + &
                                      (model%quc(m1:n1,i2,i3,i_DUSTK3)-equK3) * exp(-rateK3*dtime)
        !
        ! --- Prevent existence of dust grains with "impossible" sizes ---
        model%quc(m1:n1,i2,i3,i_DUSTK3)=min(max(model%quc(m1:n1,i2,i3,i_DUSTK3), &
                                                DUSTK0), &
                                            xnCtot)
        !
      end do ! i2
    end do ! i3
    !$OMP END DO
    !
  case('dust_k3mon_02')
    ! --- Simple dust formation (two components (nquc=2): K3, xnCtot) ---
    !
    !$OMP DO
    do i3=m3,n3
      do i2=m2,n2
        ! --- Compute dust formation rate ---
        call DUSTFORM_6(npt, &
                        model%quc(m1:n1,i2,i3,i_nCgas), Tredfactor*T(m1:n1,i2,i3), &
                        model%quc(m1:n1,i2,i3,i_DUSTK3), &
                        CzuO, xnDmax,  &
                        taui, DUSTK0, xnCtot, r1DUST, &
                        rateK3, equK3)
!!$ if (i3 == 120) then
!!$    write(6,'(A,11(1PE11.3))') 'xncgas: ', model%quc(20:30,1,i3,i_nCgas)
!!$    write(6,'(A,11(1PE11.3))') 'k3:     ', model%quc(20:30,1,i3,i_DUSTK3)
!!$    write(6,'(A,11(1PE11.3))') 'rateK3: ', rateK3(20:30)
!!$    write(6,'(A,11(1PE11.3))') 'equK3:  ', equK3(20:30)
!!$    write(6,'(A,11(1PE11.3))') 'dustk0: ', dustk0(20:30)
!!$    write(6,'(A,11(1PE11.3))') 'taui:   ', taui(20:30)
!!$    write(6,'(A,11(1PE11.3))') 'r1DUST: ', r1dust  
!!$    write(6,'(A,11(1PE11.3))') 'xnCtot: ', xnctot(20:30)
!!$    write(6,'(A,11(1PE11.3))') 'v3:     ', model%v3(20:30,1,i3)
!!$  endif
        !
        ! --- Prevent existence of dust grains with "impossible" sizes ---
        model%quc(m1:n1,i2,i3,i_DUSTK3)=min(max(model%quc(m1:n1,i2,i3,i_DUSTK3), &
                                                DUSTK0), &
                                            xnCtot)
        !
        ! --- Add source term ---
!!$        model%quc(m1:n1,i2,i3,i_DUSTK3)=model%quc(m1:n1,i2,i3,i_DUSTK3) + dtime * &
!!$                                          (taui * model%quc(m1:n1,i2,i3,i_DUSTK3)**(2.0/3.0)  &
!!$                                                * DUSTK0**(1.0/3.0))
        model%quc(m1:n1,i2,i3,i_DUSTK3)=equK3 + &
                                      (model%quc(m1:n1,i2,i3,i_DUSTK3)-equK3) * exp(-rateK3*dtime)
        !
        ! --- Prevent existence of dust grains with "impossible" sizes ---
        model%quc(m1:n1,i2,i3,i_DUSTK3)=min(max(model%quc(m1:n1,i2,i3,i_DUSTK3), &
                                                DUSTK0), &
                                            xnCtot)
        model%quc(m1:n1,i2,i3,i_nCgas)=xnCtot-model%quc(m1:n1,i2,i3,i_DUSTK3)
!!$  if (i3 == 120) then
!!$    write(6,'(A,11(1PE11.3))') '-> xncgas: ', model%quc(20:30,1,i3,i_nCgas)
!!$    write(6,'(A,11(1PE11.3))') '-> k3:     ', model%quc(20:30,1,i3,i_DUSTK3)
!!$  endif
        !
      end do ! i2
    end do ! i3
    !$OMP END DO
    !
end select ! (par%dustScheme)
!
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
end subroutine rhd_dust_K3MONSourceStep


!----------*************************----------------------------------------------------------------
subroutine rhd_dust_K3MON2SourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_K3MON2SourceStep ('rhd_dust_K3_monomer_Source_Step')
!
! PURPOSE:
!   Update model by applying dust source term step (one moment only).
!   
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_K3MON2SourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                      e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('n_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('n_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('n_3') integer, upper index, 3rd dimension
!
! ROUTINES:
!   DUSTFORM_6:            Compute condensation or evaporation rate for forsterite
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   const_module:          Global physical and mathematical constants and units
!   dust_k3mon_module:     Simplified forsterite dust formation model (1 moment type)
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Control parameters:
!     par%ar_dustGrainRadius(2): Rmax; maximum grain radius
!     par%C_dust01:              RHOD; density of grain material
!     par%C_dust02:              AGW;  atomic weight of dust monomer
!     par%C_dust03:              eps;  number fraction of rarest component, e.g. Mg2
!     par%C_dust04:              Lower dust limit fraction (for settling velocity)
!     par%C_dust05:              Minimum number of monomers during condensation
!     par%C_dust06:              STICK1; sticking coefficient
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-11-20 (B.F.) First version
!   2002-11-22 (B.F.) Because of instabilities due to short time scales
!                     approximate linear equation introduced and solved exactly
!   2003-08-27 (B.F.) Change of meaning of control parameters: par%C_dust01: CzuO
!                     par%C_dust01 -> par%C_dust02, 02 -> 03, 03 -> 04
!   2007-12-03 (B.F.) Call dustform_6 instead of dustform_5
!   2007-12-11 (B.F.) Derive routine from rhd_dust_K3MON2SourceStep,
!                     use different control parameters, source for constants, swap qucs
!   2007-12-18 (B.F.) Re-arrangement of input parameters
!   2008-04-16 (B.F.) Private i1
!   2008-04-23 (B.F.) quctot
!   2009-10-15 (B.F.) Parallelize call to eosinter
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
use const_module, only: pi, amu
use dust_k3mon_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3
integer                                 :: npt
real                                    :: CzuO, xnDmax, Tredfactor, r1DUST, Rmax, AGW, RHOD
real, dimension(model%m1:model%n1)      :: taui, DUSTK0, xnCtot, xnCgas, rateK3, equK3, quctot
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
integer                                 :: ierr0
character(len=80)                       :: outstr0
integer, parameter                      :: i_mono=1, i_dust=2
!real, parameter                         :: AGW    = 140.71, & ! atomic weight of forsterite
!                                           RHOD   = 3.3E+00   ! density of grain material
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
npt=n1-m1+1
!
! --- Solve EOS, find out Helium abundance (number fraction) ---
!
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
! --- Set to meaningless dummy value ---
CzuO  =1.8
!
! === Take basic parameters from parameter file ===
! --- Maximum dust grain radius ---
Rmax  =par%ar_dustGrainRadius(i_dust)
! --- Density of grain material ---
RHOD  =par%C_dust01
! --- Atomic weight of dust monomer ---
AGW   =par%C_dust02
!
! --- Maximumum number of monomers per grain ---
xnDmax=(4.0*pi/3.0)*Rmax**3*RHOD/(AGW*amu)
!
! --- Typical loop ---
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i1, i2, taui, DUSTK0, xnCtot, xnCgas, r1DUST, rateK3, equK3, quctot)
!
select case(par%dustScheme)
  case('dust_k3mon_03')
    ! --- Simple dust formation; two components (nquc=2): mass density of monomers and dust ---
    !
    !$OMP DO
    do i3=m3,n3
      do i2=m2,n2
        ! --- Compute dust formation rate (routine works with number densities) ---
        call DUSTFORM_6(npt, &
                        model%quc(m1:n1,i2,i3,i_mono)/(AGW*amu), T(m1:n1,i2,i3), &
                        model%quc(m1:n1,i2,i3,i_dust)/(AGW*amu), &
                        CzuO, xnDmax,  &
                        taui, DUSTK0, xnCtot, r1DUST, &
                        rateK3, equK3)
!if ((i2 == m2) .and. (i3 >= 120) .and. (i3 <= 150)) then
!  write(6,*) i3
!  write(6,'(A,11(1PE11.3))') 'v3:      ', model%v3(m1,m2,i3)
!  write(6,'(A,11(1PE11.3))') 'T:       ', T(m1,m2,i3)
!  write(6,'(A,11(1PE11.3))') 'rho:     ', model%rho(m1,m2,i3)
!!  write(6,'(A,11(1PE11.3))') 'qucmono: ', model%quc(m1,m2,i3,i_mono)
!!  write(6,'(A,11(1PE11.3))') 'qucdust: ', model%quc(m1,m2,i3,i_dust)
!  write(6,'(A,11(1PE11.3))') 'qucmono/:', model%quc(m1,m2,i3,i_mono)/(AGW*amu)
!  write(6,'(A,11(1PE11.3))') 'qucdust/:', model%quc(m1,m2,i3,i_dust)/(AGW*amu)
!  write(6,'(A,11(1PE11.3))') 'dustk0:  ', dustk0(m1)
!  write(6,'(A,11(1PE11.3))') 'xnCtot:  ', xnctot(m1)
!  write(6,'(A,11(1PE11.3))') 'rateK3:  ', rateK3(m1)
!  write(6,'(A,11(1PE11.3))') 'equK3:   ', equK3(m1)
!!  write(6,'(A,11(1PE11.3))') 'taui:    ', taui(m1)
!!  write(6,'(A,11(1PE11.3))') 'r1DUST:  ', r1dust
!!  write(6,'(A,11(1PE11.3))') 'xnDmax:  ', xnDmax
!endif
        !
        quctot=model%quc(m1:n1,i2,i3,i_mono)+model%quc(m1:n1,i2,i3,i_dust)
        !
        ! --- Prevent existence of dust grains with "impossible" sizes ---
        model%quc(m1:n1,i2,i3,i_dust)=min(max(model%quc(m1:n1,i2,i3,i_dust), &
                                              (AGW*amu)*DUSTK0), &
                                          quctot)
        !
        ! --- Add source term ---
        do i1=m1,n1
          if (abs(rateK3(i1)*dtime) > 1.0E-05) then
            model%quc(i1,i2,i3,i_dust)=((AGW*amu)*equK3(i1) + &
                                        (model%quc(i1,i2,i3,i_dust)-(AGW*amu)*equK3(i1)) * &
                                         exp(-rateK3(i1)*dtime) &
                                       )
          else
            model%quc(i1,i2,i3,i_dust)=model%quc(i1,i2,i3,i_dust)+ &
                                      (model%quc(i1,i2,i3,i_dust)-(AGW*amu)*equK3(i1))* &
                                      (-rateK3(i1)*dtime)
          endif
        end do ! i1
        !
        ! --- Prevent existence of dust grains with "impossible" sizes ---
        model%quc(m1:n1,i2,i3,i_dust)=min(max(model%quc(m1:n1,i2,i3,i_dust), &
                                              (AGW*amu)*DUSTK0), &
                                          quctot)
        model%quc(m1:n1,i2,i3,i_mono)=quctot-model%quc(m1:n1,i2,i3,i_dust)
        !
      end do ! i2
    end do ! i3
    !$OMP END DO
    !
end select ! par%dustScheme
!
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
end subroutine rhd_dust_K3MON2SourceStep


!----------***********************------------------------------------------------------------------
subroutine rhd_dust_M4C2SourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_M4C2SourceStep ('rhd_dust_moment_4_carbon_2_Source_Step')
!
! PURPOSE:
!   Update model by applying dust source term step (4 moments).
!   
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_M4C2SourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                      e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('n_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('n_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('n_3') integer, upper index, 3rd dimension
!
! ROUTINES:
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   MOMENT:                Calculate the 4 moments for carbon dust at the new time level
!   tabinter_getinfo:      Give requested information about EOS interpolation table
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   const_module:          Global physical and mathematical constants and units
!   dust_momentc2_module:  Routines for 4-moment carbon dust scheme
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!   Dust temperature is assumed to be gas temperature.
!   There is no time step control: changes in gust densities per time step are assumed to
!   be small.
!
! PROCEDURE:
!   Control parameters:
!     par%C_dust01: CzuO   (C to O ratio (typically 1.4 to 1.8))
!     par%C_dust02: OHS    (abundance of O (solar: 1.0D+01**(-3.18D+00) = 6.606934E-04))
!     par%C_dust03: epsint (cutoff for integration of degree of condensation, typically 1.0D-05) 
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-11-20 (B.F.) First version
!   2002-11-22 (B.F.) Because of instabilities due to short time scales
!                     approximate linear equation introduced and solved exactly
!   2003-08-25 (B.F. & S.H.) Integration by call of MOMENT
!   2003-08-26 (B.F. & S.H.) OpenMP calls of MOMENT
!   2003-08-27 (B.F.) Read par%C_dust01,2,3
!   2007-12-18 (B.F.) Set pmass=amu
!   2008-06-24 (B.F.) Different allocation methods for local 3D array
!   2009-10-15 (B.F.) Parallelize call to eosinter
!   2017-09-15 (B.F.) Account for possible ghost cells in input model in MPI case
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
use const_module, only: amu
use dust_momentc2_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
include 'DINDEX.INC'
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer, parameter                      :: srk=kind(0.0E+00), drk=kind(0.0D+00)
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3
integer                                 :: npt
integer, dimension(3)                   :: i_c
real                                    :: abuY, avg_c, max_c
real(kind=drk)                          :: HeH, R0dust, pmass, OHS, CzuO, consNC, SGmin, epsint, &
                                           dtime_drk
real(kind=drk), &
      dimension(model%m1:model%n1)      :: Temp, rhoG
real(kind=drk), &
      dimension(model%m1:model%n1,NEQ)  :: Yold, Ynew
real(kind=drk), &
      dimension(model%m1:model%n1,11)   :: DV_scratch
real(kind=drk), &
      dimension(model%m1:model%n1,36)   :: SCRD
!#
# if (macro_arrays_LocalVarType_Case == 0)
!#
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
!#
# elif (macro_arrays_LocalVarType_Case == 1)
!#
real, dimension(:,:,:), pointer         :: T
!#
# else
!#
real, dimension(:,:,:), allocatable     :: T
!#
# endif
!#
integer, parameter                      :: i_DUSTK0=1, i_DUSTK1=2, i_DUSTK2=3, i_DUSTK3=4
integer                                 :: ierr0
character(len=80)                       :: outstr0
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
npt=n1-m1+1
!
!#
# if (macro_arrays_LocalVarType_Case != 0)
!#
! --- Allocate temporary 3D array ---
allocate(T(m1:n1,m2:n2,m3:n3))
!#
# endif
!#
!
! --- Solve EOS: compute temperature; find out Helium abundance (number fraction) ---
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
call tabinter_getinfo(abuY=abuY)
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
! --- Set some parameters ---
!                   --- Abundance of He (default:10.D00**(-1.00D00)) ---
HeH   =abuY
!                   --- C to O ratio (typically 1.4 to 1.8) ---
CzuO  =par%C_dust01
!                   --- Abundance of O (solar: 1.0D+01**(-3.18D+00) = 6.606934E-04) ---
OHS   =par%C_dust02
!                   --- Cutoff for integration of degree of condensation (typically 1.0D-05) ---
epsint=par%C_dust03
!                   --- Minimum supersaturation ratio for nucleation ---
SGmin =3.0D+00
!                   --- Proton mass ---
!pmass =1.6726D-24
pmass =amu
!                   --- abundance ratio abuC/rho ---
consNC=(CzuO - 1.0D00) * OHS / (pmass * (1.0D00 + HeH*4.0D00))
!                   --- time step
dtime_drk=real(dtime, kind=drk)
!
!
! --- Dust formation: update the 4 dust moments ---
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i2, Temp, rhoG, Yold, Ynew, DV_scratch, SCRD) 
!
!$OMP DO
do i3=m3,n3
  do i2=m2,n2
    ! ---Transform to double precision ---
    Temp       =real(T(        m1:n1,i2,i3)         , kind=drk)
    rhoG       =real(model%rho(m1:n1,i2,i3)         , kind=drk)
    Yold(:,IK0)=real(model%quc(m1:n1,i2,i3,i_DUSTK0), kind=drk)/rhoG
    Yold(:,IK1)=real(model%quc(m1:n1,i2,i3,i_DUSTK1), kind=drk)/rhoG
    Yold(:,IK2)=real(model%quc(m1:n1,i2,i3,i_DUSTK2), kind=drk)/rhoG
    Yold(:,IK3)=real(model%quc(m1:n1,i2,i3,i_DUSTK3), kind=drk)/rhoG
    !
    ! --- Limit temperature to reasonable values ---
    Temp=min(Temp, 5000.0_drk)
    !
    ! --- Integrate moment rate equations ---
    call MOMENT(npt, Yold, dtime_drk, rhoG, Temp, consNC, &
                HeH, SGmin, epsint, &
                R0dust, Ynew, DV_scratch, SCRD)
    !
    ! --- Put temporary data back into model arrays ---
    model%quc(m1:n1,i2,i3,i_DUSTK0)=real(Ynew(:,IK0)*rhoG, kind=srk)
    model%quc(m1:n1,i2,i3,i_DUSTK1)=real(Ynew(:,IK1)*rhoG, kind=srk)
    model%quc(m1:n1,i2,i3,i_DUSTK2)=real(Ynew(:,IK2)*rhoG, kind=srk)
    model%quc(m1:n1,i2,i3,i_DUSTK3)=real(Ynew(:,IK3)*rhoG, kind=srk)
    !
  end do ! i2
end do ! i3
!$OMP END DO
!
!$OMP END PARALLEL
!
! --- Find maximum degree of condensation ---
! --- ToDo: Use OpenMP and MPI            ---
i_c  =maxloc(model%quc(m1:n1,m2:n2,m3:n3,i_DUSTK3)/model%rho) +(/m1, m2, m3/)-1
avg_c=sum(   model%quc(m1:n1,m2:n2,m3:n3,i_DUSTK3)/model%rho) / &
      real(size(model%rho)) / real(consNC,kind=srk)
max_c=model%quc(i_c(1),i_c(2),i_c(3),i_DUSTK3)/(consNC*model%rho(i_c(1),i_c(2),i_c(3)))
if (gl%nc_m >= 0) then
  write(gl%nc_m,'(A,F8.6,A,F8.6,A,3(1X,I4),A,1PE10.4,1X,1PE10.4)') &
    'DUST: dc_avg ', avg_c, &
    '  max ', max_c, &
    ' at', i_c, &
    ' w. T,rho ', T(i_c(1),i_c(2),i_c(3)), model%rho(i_c(1),i_c(2),i_c(3))
endif
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
!#
# if (macro_arrays_LocalVarType_Case != 0)
!#
! --- Decallocate temporary 3D arrays ---
deallocate(T)
!#
# endif
!#
!
end subroutine rhd_dust_M4C2SourceStep


!----------*****************************------------------------------------------------------------
subroutine rhd_dust_Mg2SiO4AGBSourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_Mg2SiO4AGBSourceStep ('rhd_dust_Mg2_Si_O4_AGB_Source_Step')
!
! PURPOSE:
!   Update model by applying dust source term step.
!   
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_Mg2SiO4AGBSourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                     e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('n_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('n_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('n_3') integer, upper index, 3rd dimension
!   i1:           ('index_1') integer, loop index for first spatial dimension (x1 direction)
!   i2:           ('index_2') integer, loop index for second spatial dimension (x2 direction)
!   i3:           ('index_3') integer, loop index for third spatial dimension (x3 direction)
!   npt:          ('number_points') integer with number of grid points in vector
!   R0dust:       ('Radius_0_dust') double-precision real
!   epsint:       ('epsilon_integration') double-precision real
!   dtime_drk:    ('delta_time_double_real_kind') double-precision version of time step
!   CNDNH:        ('CNDNH') double-precision real vector with dust grain abundances
!   Temp:         ('Temperature') double-precision real vector with gas temperature
!   rhoG:         ('rho_Gas') double-precision real vector with gas density
!   Yold:         ('Y_old') double-precicion real 2D array with
!                 monomer_number_density/gas_mass_density at old time level
!   Ynew:         ('Y_new') double-precicion real 2D array with
!                 monomer_number_density/gas_mass_density at new time level
!   DV_scratch:   ('Dust_Vector_scratch') double-precicion real 2D scratch array
!   T:            ('Temperature') real 3D array with temperature
!   ierr0:        ('integer_error_0') integer with error status
!   outstr0:      ('output_string_0') character with error report
!
! ROUTINES:
!   dust_Mg2SiO4AGB_1Dchange: Calculate the new relative amount of forsterite monomers
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   const_module:          Global physical and mathematical constants and units
!   dust_Mg2SiO4AGB_module: Calculate dust (forsterite) concentration at the new time level
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!   Counters in timing routines may be modified.
!
! RESTRICTIONS:
!   Dust temperature is assumed to be gas temperature.
!   There is no time step control: changes in gust densities per time step are assumed to
!   be small.
!
! PROCEDURE:
!   Control parameters:
!     par%C_dust03: epsint (cutoff for integration of degree of condensation, typically 1.0E-05) 
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-11-20 (B.F.) First version
!   2002-11-22 (B.F.) Because of instabilities due to short time scales
!                     approximate linear equation introduced and solved exactly
!   2003-08-25 (B.F. & S.H.) Integration by call of MOMENT
!   2003-08-26 (B.F. & S.H.) OpenMP calls of MOMENT
!   2003-08-27 (B.F.) Read par%C_dust01,2,3
!   2007-12-18 (B.F.) Set pmass=amu
!   2008-06-24 (B.F.) Different allocation methods for local 3D array
!   2009-10-15 (B.F.) Parallelize call to eosinter
!   2012-03-08 (B.F.) M4C2 -> Mg2SiO4AGB
!   2017-09-15 (B.F.) Account for possible ghost cells in input model in MPI case
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
use const_module, only: amu
use dust_Mg2SiO4AGB_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer, parameter                      :: srk=kind(0.0E+00), drk=kind(0.0D+00)
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3, iquc
integer                                 :: npt
real(kind=drk)                          :: R0dust, epsint, dtime_drk
real(kind=drk), &
      dimension(1:model%nquc)           :: CNDNH
real(kind=drk), &
      dimension(model%m1:model%n1)      :: Temp, rhoG
real(kind=drk), &
      dimension(model%m1:model%n1, &
                1:model%nquc)           :: Yold, Ynew
real(kind=drk), &
      dimension(model%m1:model%n1,7)    :: DV_scratch
!#
# if (macro_arrays_LocalVarType_Case == 0)
!#
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
!#
# elif (macro_arrays_LocalVarType_Case == 1)
!#
real, dimension(:,:,:), pointer         :: T
!#
# else
!#
real, dimension(:,:,:), allocatable     :: T
!#
# endif
!#
integer                                 :: ierr0
character(len=80)                       :: outstr0
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
npt=n1-m1+1
!
!#
# if (macro_arrays_LocalVarType_Case != 0)
!#
! --- Allocate temporary 3D array ---
allocate(T(m1:n1,m2:n2,m3:n3))
!#
# endif
!#
!
! --- Solve EOS: compute temperature; find out Helium abundance (number fraction) ---
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
! --- Set some parameters ---
!                   --- Cutoff for integration of degree of condensation (typically 1.0D-05) ---
epsint=par%C_dust03
!                   --- Vector of dust grain abundances ---
CNDNH =par%ar_dustGrainAbu(1:model%nquc)
!
!                   --- time step
dtime_drk=real(dtime, kind=drk)
!
!
! --- Dust formation: update the nquc dust types ---
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i2, iquc, Temp, rhoG, Yold, Ynew, DV_scratch) 
!
!$OMP DO
do i3=m3,n3
  do i2=m2,n2
    ! ---Transform to double precision ---
    rhoG       =real(model%rho(m1:n1,i2,i3)        , kind=drk)
    Temp       =real(T(        m1:n1,i2,i3)        , kind=drk)
    !
    ! --- Limit temperature to reasonable values ---
    Temp=min(Temp, 5000.0_drk)
    !
    ! --- Compute temporary data ---
    do iquc=1,model%nquc
      Yold(:,iquc)=real(model%quc(m1:n1,i2,i3,iquc), kind=drk)/rhoG
    end do ! iquc
    !
    ! --- Integrate rate equation ---
    call dust_Mg2SiO4AGB_1Dchange(npt, model%nquc, Yold, dtime_drk, rhoG, Temp, &
                                  CNDNH, epsint, R0dust, Ynew, DV_scratch)
    !
    ! --- Put temporary data back into model arrays ---
    do iquc=1,model%nquc
      model%quc(m1:n1,i2,i3,iquc)=real(Ynew(:,iquc)*rhoG, kind=srk)
    end do ! iquc
    !
  end do ! i2
end do ! i3
!$OMP END DO
!
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
!#
# if (macro_arrays_LocalVarType_Case != 0)
!#
! --- Decallocate temporary 3D arrays ---
deallocate(T)
!#
# endif
!#
!
end subroutine rhd_dust_Mg2SiO4AGBSourceStep


!----------***********************************------------------------------------------------------
subroutine rhd_dust_Mg2SiO4pAl2O3AGBSourceStep(model, action, dtime, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dust_Mg2SiO4pAl2O3AGBSourceStep ('rhd_dust_Mg2_Si_O4_plus_Al2_O3_AGB_Source_Step')
!
! PURPOSE:
!   Update model by applying forsterite+corundum source term step.
!   
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_dust_Mg2SiO4pAl2O3AGBSourceStep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action:       ('action') derived type with collection of control parameters
!   dtime:        ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:      ('model') derived type, contains all model data (e.g. model%rho)
!     m1:         ('m_number_1') integer, 1. dimension, start index
!     n1:         ('number_1')   integer, 1. dimension, end index
!     m2:         ('m_number_2') integer, 2. dimension, start index
!     n2:         ('number_2')   integer, 2. dimension, end index
!     m3:         ('m_number_3') integer, 3. dimension, start index
!     n3:         ('number_3')   integer, 3. dimension, end index
!     nghost1:    ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:    ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:    ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:      ('index_time') integer, time step number
!     nquc:       ('number_quantities_center') integer, number of additional fields
!     time:       ('time') real, time
!     xc1:        ('x_center_1') real, 1. coordinate of cell centers
!     xc2:        ('x_center_2') real, 2. coordinate of cell centers
!     xc3:        ('x_center_3') real, 3. coordinate of cell centers
!     xb1:        ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:        ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:        ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:        ('rho') real, dimension(:,:,:), density
!     ei:         ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:         ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:         ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:         ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:        ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities,
!                                     e.g., quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!   m1:           ('m_1') integer, lower index, 1st dimension
!   n1:           ('number_1') integer, upper index, 1st dimension
!   m2:           ('m_2') integer, lower index, 2nd dimension
!   n2:           ('number_2') integer, upper index, 2nd dimension
!   m3:           ('m_3') integer, lower index, 3rd dimension
!   n3:           ('number_3') integer, upper index, 3rd dimension
!   i1:           ('index_1') integer, loop index for first spatial dimension (x1 direction)
!   i2:           ('index_2') integer, loop index for second spatial dimension (x2 direction)
!   i3:           ('index_3') integer, loop index for third spatial dimension (x3 direction)
!   npt:          ('number_points') integer with number of grid points in vector
!   R0dust:       ('Radius_0_dust') double-precision real
!   epsint:       ('epsilon_integration') double-precision real
!   dtime_drk:    ('delta_time_double_real_kind') double-precision version of time step
!   CNDNH:        ('CNDNH') double-precision real vector with dust grain abundances
!   Temp:         ('Temperature') double-precision real vector with gas temperature
!   rhoG:         ('rho_Gas') double-precision real vector with gas density
!   Yold:         ('Y_old') double-precicion real 2D array with
!                 monomer_number_density/gas_mass_density at old time level
!   Ynew:         ('Y_new') double-precicion real 2D array with
!                 monomer_number_density/gas_mass_density at new time level
!   DV_scratch:   ('Dust_Vector_scratch') double-precicion real 2D scratch array
!   T:            ('Temperature') real 3D array with temperature
!   ierr0:        ('integer_error_0') integer with error status
!   outstr0:      ('output_string_0') character with error report
!
! ROUTINES:
!   dust_Al2O3AGB_1Dchange: Calculate the new relative amount of corundum monomers
!   dust_Mg2SiO4AGB_1Dchange: Calculate the new relative amount of forsterite monomers
!   eosinter_TOMP:         Solve EOS by bicubic interpolation: rho, e -> T
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   const_module:          Global physical and mathematical constants and units
!   dust_Al2O3AGB_module:  Calculate dust (corundum) concentration at the new time level
!   dust_Mg2SiO4AGB_module: Calculate dust (forsterite) concentration at the new time level
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   tabinter_module:       Routines for handling of interpolation coefficients
!   timing_module:         Timing routines and data for performance measurements and benchmarks
!
! SIDE EFFECTS:
!   Counters in timing routines may be modified.
!
! RESTRICTIONS:
!   Dust temperature is assumed to be gas temperature.
!   There is no time step control: changes in gust densities per time step are assumed to
!   be small.
!
! PROCEDURE:
!   Control parameters:
!     par%C_dust03: epsint (cutoff for integration of degree of condensation, typically 1.0E-05) 
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-11-20 (B.F.) First version
!   2002-11-22 (B.F.) Because of instabilities due to short time scales
!                     approximate linear equation introduced and solved exactly
!   2003-08-25 (B.F. & S.H.) Integration by call of MOMENT
!   2003-08-26 (B.F. & S.H.) OpenMP calls of MOMENT
!   2003-08-27 (B.F.) Read par%C_dust01,2,3
!   2007-12-18 (B.F.) Set pmass=amu
!   2008-06-24 (B.F.) Different allocation methods for local 3D array
!   2009-10-15 (B.F.) Parallelize call to eosinter
!   2012-03-08 (B.F.) M4C2 -> Mg2SiO4AGB
!   2017-05-22 (B.F. Uppsala) Adapt from dust_Mg2SiO4AGB_SourceStep
!   2017-09-15 (B.F.) Account for possible ghost cells in input model in MPI case
!   2017-12-13 (B.F.) Bugfix: Change parameter in calls to dust_XXX_1Dchange: model%nquc -> 1
!---------------------------------------------------------------------------------------------------
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
use const_module, only: amu
use dust_Al2O3AGB_module
use dust_Mg2SiO4AGB_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
type(box_type),            intent(inout):: model
type(action_type),         intent(in)   :: action
real,                      intent(in)   :: dtime
integer,         optional, intent(out)  :: ierr
character(len=*),optional, intent(out)  :: outstr
!
! --- Local variables ---
integer, parameter                      :: srk=kind(0.0E+00), drk=kind(0.0D+00)
integer                                 :: m1,n1, m2,n2, m3,n3, i1, i2, i3, iquc
integer                                 :: npt
real(kind=drk)                          :: R0dust, epsint, dtime_drk
real(kind=drk), &
      dimension(1:model%nquc)           :: CNDNH
real(kind=drk), &
      dimension(model%m1:model%n1)      :: Temp, rhoG
real(kind=drk), &
      dimension(model%m1:model%n1, &
                1:model%nquc)           :: Yold, Ynew
real(kind=drk), &
      dimension(model%m1:model%n1,7)    :: DV_scratch
!#
# if (macro_arrays_LocalVarType_Case == 0)
!#
real, dimension(model%m1:model%n1, &
                model%m2:model%n2, &
                model%m3:model%n3)      :: T
!#
# elif (macro_arrays_LocalVarType_Case == 1)
!#
real, dimension(:,:,:), pointer         :: T
!#
# else
!#
real, dimension(:,:,:), allocatable     :: T
!#
# endif
!#
integer                                 :: ierr0
character(len=80)                       :: outstr0
!---------------------------------------------------------------------------------------------------
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- Make abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
npt=n1-m1+1
!
!#
# if (macro_arrays_LocalVarType_Case != 0)
!#
! --- Allocate temporary 3D array ---
allocate(T(m1:n1,m2:n2,m3:n3))
!#
# endif
!#
!
! --- Solve EOS: compute temperature; find out Helium abundance (number fraction) ---
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST: eosinter')
!#
# endif
!#
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, &
                   m1,n1, m2,n2, m3,n3, model%nghost1, model%nghost2, model%nghost3)
!$OMP END PARALLEL
!#
# ifdef rhd_dust_t01
!#
call timing_stop( 'DUST: eosinter')
call timing_start('DUST: loop')
!#
# endif
!#
!
! --- Set some parameters ---
!                   --- Cutoff for integration of degree of condensation (typically 1.0D-05) ---
epsint   =real(par%C_dust03, kind=drk)
!                   --- Vector of dust-grain abundances ---
CNDNH    =real(par%ar_dustGrainAbu(1:model%nquc), kind=drk)
!
!                   --- time step
dtime_drk=real(dtime, kind=drk)
!
!
! --- Dust formation: update the nquc dust types ---
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP   PRIVATE(i2, iquc, Temp, rhoG, Yold, Ynew, DV_scratch) 
!
!$OMP DO
do i3=m3,n3
  do i2=m2,n2
    ! ---Transform to double precision ---
    rhoG       =real(model%rho(m1:n1,i2,i3)        , kind=drk)
    Temp       =real(T(        m1:n1,i2,i3)        , kind=drk)
    !
    ! --- Limit temperature to reasonable values ---
    Temp=min(Temp, 5000.0_drk)
    !
    ! --- Compute temporary data ---
    do iquc=1,model%nquc
      Yold(:,iquc)=real(model%quc(m1:n1,i2,i3,iquc), kind=drk)/rhoG
    end do ! iquc
    !
    ! --- Integrate rate equation for forsterite (Mg2SiO4) ---
    iquc=1
    call dust_Mg2SiO4AGB_1Dchange(npt, 1, Yold(:,iquc:iquc), dtime_drk, rhoG, Temp, &
                                  CNDNH(iquc:iquc), epsint, R0dust, Ynew(:,iquc:iquc), DV_scratch)
    !
    if (model%nquc >= 2) then
    ! --- Integrate rate equation for corundum (Al2O3) ---
      iquc=2
      call dust_Al2O3AGB_1Dchange(npt, 1, Yold(:,iquc:iquc), dtime_drk, rhoG, Temp, &
                                  CNDNH(iquc:iquc), epsint, R0dust, Ynew(:,iquc:iquc), DV_scratch)
    !
    endif
    !
    ! --- Put temporary data back into model arrays ---
    do iquc=1,model%nquc
      model%quc(m1:n1,i2,i3,iquc)=real(Ynew(:,iquc)*rhoG, kind=srk)
    end do ! iquc
    !
  end do ! i2
end do ! i3
!$OMP END DO
!
!$OMP END PARALLEL
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST: loop')
!#
# endif
!#
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
!#
# if (macro_arrays_LocalVarType_Case != 0)
!#
! --- Decallocate temporary 3D arrays ---
deallocate(T)
!#
# endif
!#
!
end subroutine rhd_dust_Mg2SiO4pAl2O3AGBSourceStep

end module rhd_dust_module
