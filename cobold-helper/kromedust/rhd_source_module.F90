!***************************************************************************************************
!
!   #####  #    # #####          ####   ####  #    # #####   ####  ###### 
!   #    # #    # #    #        #      #    # #    # #    # #    # #      
!   #    # ###### #    #         ####  #    # #    # #    # #      #####  
!   #####  #    # #    #             # #    # #    # #####  #      #      
!   #   #  #    # #    #        #    # #    # #    # #   #  #    # #      
!   #    # #    # #####  ######  ####   ####   ####  #    #  ####  ###### _module
!
!   Radiation-Hydrodynamics-Code: source-term master module
!
!***************************************************************************************************
!   This module provides the master source-term routine that calls the individual routines
!   to handle the source terms to treat dust formation, molecular networks, or hydrogen ionization.
!***************************************************************************************************
!   Fortran90
!   2020-01-25
!***************************************************************************************************
! With contributions from:
!   Bernd Freytag:         Uppsala, Lyon
!   Matthias Steffen:      Potsdam
!   Susanne Hoefner:       Uppsala
!   Sven Wedemeyer-Boehm:  Freiburg
!   Hans-Guenter Ludwig:   Kiel, Lund, Meudon, Heidelberg
!***************************************************************************************************
! MODULES:
!   rhd_source_module:     Provide master routine for source terms due to dust or molecules
!
!***************************************************************************************************


!------*****************----------------------------------------------------------------------------
module rhd_source_module
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_source_module ('rhd_source_module')
!
! PURPOSE:
!   Provide master routine for source terms due to dust or molecules.
!
! CATEGORY:
!   Hydrodynamics, Dust
!
! CALLING SEQUENCE:
!   use rhd_source_module
!
! VARIABLES:
!
! ROUTINES: (contained)
!   rhd_source_SwitchInfo: Print information about used compiler switches
!   rhd_source_PropInit:   Initialize property structure and set e.g. number of dust components
!   rhd_source_Step:       Update model by applying source term step
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
!   2019-11-24 (B.F.) Extract routines from rhd_dust_module and rename them
!---------------------------------------------------------------------------------------------------
implicit none
!
contains

!----------*******************----------------------------------------------------------------------
subroutine rhd_source_SwitchInfo(ncp)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_source_SwitchInfo ('rhd_source_Switch_Information')
!
! PURPOSE:
!   Print information about used compiler switches.
!
! CATEGORY:
!   Dust, Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_source_SwitchInfo(ncp)
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
!   2011-11-11 (B.F. Uppsala) Adapted from rhd_shortrad_SwitchInfo
!   2011-11-24 (B.F.) Adapted from rhd_dust_SwitchInfo
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
write(ncp,'(A)')    'Compiler switches: rhd_source_module ......................'
!
!#
# ifdef rhd_chem01
!#
write(ncp,'(A,A)' ) '  rhd_'//'chem01:                         ',defstr
!#
# else
!#
write(ncp,'(A,A)' ) '  rhd_'//'chem01:                         ',undefstr
!#
# endif
!#
!
!#
# ifdef rhd_krome01
!#
write(ncp,'(A,A)' ) '  rhd_'//'krome01:                         ',defstr
!#
# else
!#
write(ncp,'(A,A)' ) '  rhd_'//'krome01:                         ',undefstr
!#
# endif
!
!#
# ifdef rhd_hion01
!#
write(ncp,'(A,A)' ) '  rhd_'//'hion01:                         ',defstr
!#
# else
!#
write(ncp,'(A,A)' ) '  rhd_'//'hion01:                         ',undefstr
!#
# endif
!#
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
end subroutine rhd_source_SwitchInfo


!----------*******************----------------------------------------------------------------------
subroutine rhd_source_PropInit(quc_prop)
!---------------------------------------------------------------------------------------------------
! NAME:
!    rhd_source_PropInit ('rhd_source_Properties_Init')
!
! PURPOSE:
!   Initialize property structure and set e.g. number of dust components.
!
! CATEGORY:
!   box handling
!
! CALLING SEQUENCE:
!   call rhd_source_PropInit(quc_prop)
!
! INPUT:
!   quc_prop:     ('quantities_centered_properties')
!
! OUTPUT:
!   None
!
! ROUTINES:
!   chem_rn_init:          Read chemistry data, initialize quantities
!   DKINIT:                Set parameters and compute Spline coefficients for C-rich chemistry
!
! MODULES:
!   chem_rn_module:        Routines for handling chemical reaction networks
!   krome_rn_module:       External chemical reaction network solver
!   dust_momentc2_module:  Routines for 4-moment carbon dust scheme
!   rhd_dust_module:       Routines for handling of source terms due to dust or molecules
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   rhd_prop_module:       Properties of quantities in 'box' structure
!
! SIDE EFFECTS:
!   Components of global variable 'prop' are modified.
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2002-10-29 (B.F. Uppsala) Written
!   2002-11-08 (B.F.) Dust property settings moved into rhd_dust_module
!   2003-08-26 (B.F. & S.H.) Call DKINIT
!   2003-08-26 (B.F. & S.H.) Call DKINIT
!   2005-01-25 (S.W.-B.) Case 'chemreacnet'
!   2007-01-27 (B.F. Lyon) Set quc_bound to 2 instead of 1
!   2007-06-04 (B.F.) Case 'dust_bins_01'
!   2007-12-11 (B.F.) Case('dust_k3mon_02'): swap quc_advect: 1<-->2,
!                     case('dust_k3mon_03') added
!                     Numbering for quc_advect changed
!   2012-03-07 (B.F.) Branch for dust scheme 'dust_Mg2SiO4AGB'
!   2012-07-19 (B.F.) Write to channel gl%nc_m only if >=0
!   2014-12-05 (B.F.) Case 'nosource' can handle more than one dust type
!   2017-05-22 (B.F.) Set array quc_positive to control check for positivity of quc arrays
!   2019-11-24 (B.F.) Adapted from rhd_dust_PropInit; stop with error message f. unknown dust scheme
!   2021-06-16 (S.A.D.) Case 'krome'
!---------------------------------------------------------------------------------------------------
use rhd_gl_module
use rhd_prop_module
use rhd_dust_module
use dust_momentc2_module
!#
# ifdef rhd_chem01
!#
use chem_rn_module
!#
# endif 
!#
# ifdef rhd_krome01
!#
use krome_rn_module
!#
# endif
!#
!
implicit none
!
! --- I/O parameters ---
character(len=*),          intent(in)   :: quc_prop
!
! --- Local variables ---
integer                                 :: iquc
character(len=80)                       :: qucstr
!---------------------------------------------------------------------------------------------------
!
prop%quc_prop=quc_prop
!
prop%nquc=0
nullify(prop%quc_ident, prop%quc_name  , prop%quc_unit, &
        prop%quc_bound, prop%quc_advect, prop%quc_positive)
! \label{f90:dustscheme__use3}
select case(quc_prop)
  case('', 'none', 'None')
    ! --- No additional densities (and no source terms) ---
    !
  case('nosource', 'dust_simple_01', 'co_component01_01', 'dust_k3mon_01')
    ! --- No source terms at all ---
    !
    if (par%N_dustGrainRadius > 0) then
      prop%nquc=par%N_dustGrainRadius
    else
      prop%nquc=1
    endif
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:prop%nquc)=2
    !
    prop%quc_advect(1:prop%nquc)=1
    !
    prop%quc_positive(1:prop%nquc)=0
    !
  case('dust_k3mon_02')
    ! --- Simplified 'dust' formation model: one moment + CO/Mg2SiO4 in gas phase ---
    !
    prop%nquc=2
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:2)=2
    !
    prop%quc_advect(1:1)=2
    prop%quc_advect(2:2)=1
    !
    prop%quc_positive(1:prop%nquc)=1
    !
  case('dust_k3mon_03')
    ! --- Simplified 'dust' formation model: Mg2SiO4 in gas phase + one moment ---
    !
    prop%nquc=2
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:2)=2
    !
    prop%quc_advect(1:1)=1
    prop%quc_advect(2:2)=3
    !
    prop%quc_positive(1:prop%nquc)=1
    !
  case('dust_moment04_c2')
    ! --- C-rich dust chemistry, 4 moments ---
    !
    prop%nquc=4
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:prop%nquc)=2
    !
    prop%quc_advect(1:prop%nquc)=1
    !
    prop%quc_positive(1:prop%nquc)=1
    !
    ! --- Initialize spline coefficients (C-rich dust chemistry) ---
    call DKINIT()
    !
  case('dust_moment04_02')
    ! --- 4 moment dust formation ---
    !
    prop%nquc=5
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:4)=2
    prop%quc_bound(5:5)=0
    !
    prop%quc_advect(1:4)=1
    prop%quc_advect(5:5)=0
    !
    prop%quc_positive(1:prop%nquc)=1
    !
  case('dust_bins_01')
    ! --- Dust bin module: one bin for each dust grain radius (+ one for monomers) ---
    !
    prop%nquc=par%N_dustGrainRadius
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:prop%nquc)=2
    !
    prop%quc_advect(1:1)=1
    prop%quc_advect(2:prop%nquc)=4
    !
    prop%quc_positive(1:prop%nquc)=1
    !
  case('dust_Mg2SiO4AGB')
    ! --- Multiple sizes of forsterite (Mg2SiO4), non-interacting ---
    !
    prop%nquc=par%N_dustGrainAbu
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:prop%nquc)=2
    !
    prop%quc_advect(1:prop%nquc)=1
    !
    prop%quc_positive(1:prop%nquc)=0
    !
  case('dust_Mg2SiO4pAl2O3AGB')
    ! --- Forsterite and corundum (Mg2SiO4 and Al2O3), non-interacting ---
    !
    prop%nquc=par%N_dustGrainAbu
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    prop%quc_bound(1:prop%nquc)=2
    !
    prop%quc_advect(1:prop%nquc)=1
    !
    prop%quc_positive(1:prop%nquc)=0
    !
!#
# ifdef rhd_chem01
!#
    !
  case('chemreacnet')
    ! --- CO network ---
    !
    ! --- Warning: Use of test parameter!!! ---
    prop%nquc=par%N_test1
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    ! --- Control boundary condition: no effect in current version ---
    ! --- B.F. (2007-01-27): That was only true because of an evil side effect:
    ! --- The prop structure was modified within the rhd_bound_3Dinoutflow routine.
    ! --- Removed: 2007-01-27, now the array has probably to be set to 2 instead of 0 ---
    prop%quc_bound(1:prop%nquc)=2
    !
    ! --- Control advection of particle densities ---
    prop%quc_advect(1:prop%nquc)=int(par%C_dust01)
    !
    prop%quc_positive(1:prop%nquc)=1
    !
    ! --- Initialise ---
    !
    call chem_rn_init()
    !
!#
# endif
!#
# ifdef rhd_krome01
!#
    !
  case('krome')
    ! --- KROME external solver ---
    !
    ! --- Warning: Use of test parameter!!! ---
    ! Does KROME still need all this stuff? It's pretty standalone, has its own rates
    ! and solver, QUC arrays defined in start model and left up to user to ensure
    ! they are defined in the correct manner to match KROME. What's left to be done
    ! here?
    prop%nquc=par%N_test1
    !
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    ! --- Control boundary condition: no effect in current version ---
    ! --- B.F. (2007-01-27): That was only true because of an evil side effect:
    ! --- The prop structure was modified within the rhd_bound_3Dinoutflow routine.
    ! --- Removed: 2007-01-27, now the array has probably to be set to 2 instead of 0 ---
    prop%quc_bound(1:prop%nquc)=2
    !
    ! --- Control advection of particle densities ---
    prop%quc_advect(1:prop%nquc)=int(par%C_dust01)
    !
    prop%quc_positive(1:prop%nquc)=1
    !
    ! --- Initialise ---
    !
    ! Perhaps call krome_init() here?
    call chem_rn_init()
    !
!#
# endif
!#      /* rhd_chem01 */
    !
!#
# ifdef rhd_hion01
!#
    !
  case('hion')
    ! --- Hydrogen ionization, time independent ---
    !
    if (gl%nc_m >= 0) then
      write(gl%nc_m,*) 'Time independent HION chosen. No advection needed.'
    endif
    !
  case('hiontd')
    ! --- Hydrogen ionization, time dependent ---
    !
    if (gl%nc_m >= 0) then
      write(gl%nc_m,*) 'Time dependent HION chosen.'
    endif
    !
    prop%nquc=6 ! HARD WIRED CHANGE ASAP
    allocate(prop%quc_ident(prop%nquc), prop%quc_name(  prop%nquc), prop%quc_unit(    prop%nquc), &
             prop%quc_bound(prop%nquc), prop%quc_advect(prop%nquc), prop%quc_positive(prop%nquc)  )
    !
    ! --- Control boundary condition: no effect in current version --- 
    prop%quc_bound(1:prop%nquc)=0 
    ! --- Control advection of particle densities ---
    prop%quc_advect(1:prop%nquc)=1
    ! --- Control positivity test ---
    prop%quc_positive(1:prop%nquc)=1
    !
!#
# endif
!#      /* rhd_hion01 */
    !
  case default
    ! --- Stop execution with error message ---
    !
    if (gl%nc_m >= 0) then
      write(gl%nc_m,*) 'Unknown dust scheme: ', adjustl(quc_prop)
      write(gl%nc_m,*) '  prop%nquc', prop%nquc
    endif
    !
    stop 'rhd_source_PropInit: ERROR: unknown dust scheme'
    !
end select ! quc_prop
!
end subroutine rhd_source_PropInit


!----------***************--------------------------------------------------------------------------
subroutine rhd_source_Step(model, action, dtime, dtime_limit, outstr,ierr)
!---------------------------------------------------------------------------------------------------
! NAME:
!   rhd_source_Step ('rhd_source_Step')
!
! PURPOSE:
!   Update model by applying source-term step.
!   
! CATEGORY:
!   Hydrodynamics
!
! CALLING SEQUENCE:
!   call rhd_source_Step(model, action, dtime, dtime_limit, outstr=outstr,ierr=ierr)
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
!   dtime_limit:  ('delta_time_limit') real, proposed maximum time step
!   outstr:       ('output_string') character, optional
!   ierr:         ('index_error') integer, optional error number:
!                       0: normal run
!                    1-99: "small" error (e.g. time step limit violation)
!                   >=100: "severe" error:
!
! VARIABLES:
!
! ROUTINES:
!   chem_rn_SourceStep:    Update model by applying CO source term step
!   dust_bins_SourceStep:  Update model by applying dust source term step (dust bin scheme)
!   hion_sourcestep:       Update model by applying Hydrogen-ionization source terms
!   rhd_dust_COSourceStep: Update model by applying CO source term step
!   rhd_dust_K3MON2SourceStep: Update model by applying dust source term step (one moment only)
!   rhd_dust_K3MONSourceStep: Update model by applying dust source term step (one moment only)
!   rhd_dust_M4C2SourceStep: Update model by applying dust source term step (4 moments)
!   rhd_dust_Mg2SiO4AGBSourceStep: Update model by applying forsterite source term step
!   rhd_dust_Mg2SiO4pAl2O3AGBSourceStep: Update model by applying forsterite+corundum source terms
!   rhd_dust_simpleSourceStep: Update model by applying simple (test) dust source term step
!   timing_start:          Start timer for specified block
!   timing_stop:           Stop timer for specified block
!
! MODULES:
!   chem_rn_module:        Routines for handling chemical reaction networks
!   dust_bins_module:      Dust-bin scheme
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   hion_main_module:      Chromospheric-radiative-transfer routines HION
!   rhd_action_module:     Routines to handle the control structure 'action'
!   rhd_box_module:        Box handling routines for RHD
!   rhd_dust_module:       Routines for handling of source terms due to dust or molecules
!   rhd_gl_module:         Global type definitions and parameters for RHD
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
!   2002-11-20 (B.F.) Separate subroutines for individual dust/molecule models
!   2005-01-25 (S.W.-B.) Case 'chemreacnet'
!   2007-06-04 (B.F.) Case 'dust_bins_01'
!   2012-03-08 (B.F. Lyon) Case for forsterite in AGB stars
!   2017-05-22 (B.F. Uppsala) Case for forsterite+corundum in AGB stars
!   2019-11-24 (B.F.) Adapted from rhd_dust_SourceStep
!---------------------------------------------------------------------------------------------------
use rhd_gl_module
use rhd_action_module
use rhd_box_module
use rhd_dust_module
use dust_bins_module
!#
# ifdef rhd_dust_t01
!#
use timing_module
!#
# endif
!#
# ifdef rhd_chem01
!#
use chem_rn_module
!#
# endif
!#
# ifdef rhd_krome01
!#
use krome_rn_module
!#
# endif
!#
# ifdef rhd_hion01
!#
use hion_main_module
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
real,                      intent(out)  :: dtime_limit
integer,          optional,intent(out)  :: ierr
character(len=*), optional,intent(out)  :: outstr
!
! --- Local variables ---
integer                                 :: ierr0
character(len=80)                       :: outstr0
!---------------------------------------------------------------------------------------------------
!
!#
# ifdef rhd_dust_t01
!#
call timing_start('DUST:')
!#
# endif
!#
!
if (gl%nc_m >= 0) then
  write(gl%nc_m,'(A)') 'Dust source step: '//trim(par%dustScheme)
endif
!
! --- No error occurred ---
ierr0=0
outstr0=''
!
! --- No limit on time step yet ---
dtime_limit=0.1*huge(dtime_limit)
!
! --- Case distinction for different dust models ---\label{f90:dustscheme__use4}
select case(par%dustScheme)
  case('nosource')
    ! --- No source terms at all ---
    !
  case('dust_simple_01')
    ! --- Simple test 'dust' formation ---
    !
    call rhd_dust_simpleSourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('co_component01_01')
    ! --- Simple CO formation ---
    !
    call rhd_dust_COSourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('dust_k3mon_01', 'dust_k3mon_02')
    ! --- Simplified dust formation model ---
    !
    call rhd_dust_K3MONSourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('dust_k3mon_03')
    ! --- Simplified dust formation model ---
    !
    call rhd_dust_K3MON2SourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('dust_moment04_c2')
    ! --- 4 moment dust formation ---
    !
    call rhd_dust_M4C2SourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('dust_Mg2SiO4AGB')
    ! --- Multiple sizes of forsterite (Mg2SiO4), non-interacting ---
    !
    call rhd_dust_Mg2SiO4AGBSourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('dust_Mg2SiO4pAl2O3AGB')
    ! --- Forsterite and corundum (Mg2SiO4 and Al2O3), non-interacting ---
    !
    call rhd_dust_Mg2SiO4pAl2O3AGBSourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
  case('dust_bins_01')
    ! --- Dust bin module: one bin for each dust grain radius (+ one for monomers) ---
    !
    call dust_bins_SourceStep(model, action, dtime, outstr=outstr0,ierr=ierr0)
    !
!#
# ifdef rhd_chem01
!#
  case('chemreacnet')
    ! --- CO network ---
    !
    call chem_rn_SourceStep(model, action, dtime, outstr=outstr0, ierr=ierr0)
!#
# endif
!#
# ifdef rhd_krome_01
!#
  case('krome')
    ! --- KROME solver ---
    !
    call krome_rn_sourcestep(model, action, dtime, outsr=outsr0, ierr=ierr0)
!#
# endif
!#
# ifdef rhd_hion01
!#
  case('hion')
    ! --- Hydrogen ionization, time independent ---
    !
    call hion_sourcestep(model,action, dtime, outstr=outstr0,ierr=ierr0)    
    !
  case('hiontd')
    ! --- Hydrogen ionization, time dependent ---
    !
    call hion_sourcestep(model,action, dtime, outstr=outstr0,ierr=ierr0)    
    !
!#
# endif
!#
!
end select ! (par%dustScheme)
!
if (present(outstr)) outstr=outstr0
if (present(ierr)) ierr=ierr0
!
!#
# ifdef rhd_dust_t01
!#
call timing_stop('DUST:')
!#
#endif

end subroutine rhd_source_Step

end module rhd_source_module
