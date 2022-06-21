!**************************************************************************************************
!
!    #####  #    #  #####  #     #          #####   #    #      
!   #       #    #  #      ##   ##          #    #  ##   #    
!   #       ######  ####   # # # #          #####   # #  #    
!   #       #    #  #      #  #  #          #  #    #  # #    
!   #       #    #  #      #     #          #   #   #   ##    
!    #####  #    #  #####  #     #  ######  #    #  #    # _module
!
!   Radiation-Hydrodynamics-Code: chemistry module for reaction networks
!
!**************************************************************************************************
!   Fortran90
!
!   SW  Sven Wedemeyer, Freiburg/Oslo
!   BF  Bernd Freytag, Uppsala
!   MS  Matthias Steffen, Potsdam
!   MSz Mikolaj Szydlarski, Oslo
!
!   2005-11-21 (SW) First full version
!   2006-08-31 (SW) error handler and minor changes 
!   2008-08-12 (MS, BF) eosinter_T gets additional arguments (m1,n1,...)
!   2019-11-24 (MS, BF) eosinter_TOMP gets additional arguments (nghost1,...)
!   2020-01-25 (BF) Minor updates, first steps towards OpenMP 
!   2021-04-30 (MSz, SW) Usage of updated DVODE package, OpenMP improvements & bug fixes
!**************************************************************************************************
! MODULES:
!   chem_rn_module: routines for handling chemical reaction networks
!
!**************************************************************************************************


!------**************------------------------------------------------------------------------------
module krome_rn_module
!--------------------------------------------------------------------------------------------------
! NAME:
!   chem_rn_module ('chem_rn_module')
!
! PURPOSE:
!   Provide routines for handling chemical reaction networks
!
! CATEGORY:
!   Hydrodynamics, Dust/Molecules
!
! CALLING SEQUENCE:
!
! VARIABLES:
!   (see definition section below)
!
! RESTRICTIONS:
!   ATTENTION! Unit number for reaction table is always 30. To be changed!
!
! ROUTINES: (contained)
!   subroutine chem_rn_readchem()
!   subroutine chem_rn_analysespecies()
!   subroutine chem_rn_getrates(rmetal, rT, outstr, ierr)
!   subroutine chem_rn_assign_density()
!   subroutine chem_rn_init():
!              Read chemistry data, initialise quantities
!   subroutine chem_rn_sourcestep(model, action, dtime, outstr,ierr):
!              Update model by applying source term step
!   subroutine chem_rn_jac (neq, T, y, ml, mu, pd, nrowpd, rpar, ipar)
!   subroutine chem_rn_func (neq, T, y, ydot, rpar, ipar)
!   subroutine chem_rn_get_speciesindex(creqspec, iisp)
!   subroutine chem_rn_get_speciesindices(creqspec, nisp, iisp, nmult)
!   subroutine chem_rn_init_co_indices()
!   subroutine chem_rn_error()
!
! MODIFICATION HISTORY:
!   09.02.04 (SW) First version
!   17.02.05 (SW) revised version
!   10.06.05 (SW) changed call of routine DVODE
!   08.08.05 (SW) pointer instead of allocatable attribute in derived types
!   17.08.06 (SW) error handler
!   2020-01-25 (B.F.) Change layout of some comments; remove blanks; add () after subroutines
!--------------------------------------------------------------------------------------------------
!
implicit none
!
save 
!
!--------------------------------------------------------------------------------
! The parameters below refer to the numbers of reactants (educts) and products
! of chemical reactions according to UMIST standard. The input file is expected 
! to match this format. 
!--------------------------------------------------------------------------------
integer, parameter ::  nreactants =  3
integer, parameter ::  nproducts  =  4
!--------------------------------------------------------------------------------
! --- parameters for single and double precision reals ---
!
integer, parameter                      :: srk=kind(0.0E+00), drk=kind(0.0D+00)
!integer, parameter :: drk=selected_real_kind(10,200)
!
!--------------------------------------------------------------------------------
! --- definition of special species/processes --- 
!
integer, parameter                      :: ncdt_sp=11
character(len=7), dimension(11)         :: ccdt_spname
character(len=80), dimension(11)        :: ccdt_spdesc
integer, dimension(11)                  :: icdt_id
!
DATA ccdt_spname( 1), icdt_id( 1), ccdt_spdesc( 1) /'e-     ', -70, 'electron recombination reaction'/      
DATA ccdt_spname( 2), icdt_id( 2), ccdt_spdesc( 2) /'CRP    ',  -1, 'cosmic-ray ionisation'/
DATA ccdt_spname( 3), icdt_id( 3), ccdt_spdesc( 3) /'CRPHOT ',  -2, 'cosmic-ray-induced photoreaction'/
DATA ccdt_spname( 4), icdt_id( 4), ccdt_spdesc( 4) /'PHOTON ',  -3, 'photoreaction 1'/
DATA ccdt_spname( 5), icdt_id( 5), ccdt_spdesc( 5) /'PHOT   ',  -4, 'photoreaction 2'/
DATA ccdt_spname( 6), icdt_id( 6), ccdt_spdesc( 6) /'CRPA   ',  -5, 'special for some reactions'/
DATA ccdt_spname( 7), icdt_id( 7), ccdt_spdesc( 7) /'FORM   ', -10, 'H2 formation on dust, freezing out on dust'/
DATA ccdt_spname( 8), icdt_id( 8), ccdt_spdesc( 8) /'EVAP   ', -20, 'evaporation from dust'/
DATA ccdt_spname( 9), icdt_id( 9), ccdt_spdesc( 9) /'RADASS ', -25, 'radiative association'/
DATA ccdt_spname(10), icdt_id(10), ccdt_spdesc(10) /'M      ', -35, 'catalytic reaction with metal'/
DATA ccdt_spname(11), icdt_id(11), ccdt_spdesc(11) /'H2*    ', -80, 'vibrationally excited molecular hydrogen'/
!
! --- list of know chemical symbols ---
integer                                 :: nchsym = 103
character*(2), dimension(103)           :: cchsym =     & 
  (/ 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
     'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
     'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
     'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
     'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
     'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Ce', &
     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
     'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
     'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
     'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
     'Md','No','Lr' /)
!
!--------------------------------------------------------------------------------
! --- type definitions ---
!
! --- type for reaction network ---
type reacnet_type
  ! contains data taken from the input file and from that derived further 
  ! information 
  ! ---------------------------------------------------------------------
  ! nreac          : number of reactions
  ! nelements      : number of atomic constituents (chemical elements)
  ! nspecies       : number of chemical species
  !
  ! --- data from input data ---
  ! id             : reaction ID
  ! alpha          : reaction coefficient alpha
  ! beta           : reaction coefficient beta
  ! gamma          : reaction coefficient gamma
  ! desc           : description (e.g., reference) of the reactions in input file
  !
  ! --- derived information on reactions ---
  ! reactype       : reaction type
  ! ispidrat       : indices of reactant species for calculation of total rate
  !
  ! --- derived information on chemical species ---
  ! natomconst     : number of atomic constiuents for each chemical species 
  !                  (0: this element is not contained)
  ! natomspsp      : total number of atomic constituents in each molecule/species
  !                  (needed for total particle densities, e.g. natomspsp=2 for CO)
  ! eleccont       : charge of molecule/species (needed for electron equation)
  ! speciesname    : molecular species (1..nspecies), names of involved chemical 
  !                  species in reaction table file
  ! elementname    : (character) atomic contituents (chemical elements) 1..nelements
  !                  (names of involved chemical elements)
  ! iatomeqsp      : index of species that are identical with a single element 
  ! ---------------------------------------------------------------------
  !
  integer                                      :: nreac, nspecies, nelements 
  integer, dimension(:),   pointer             :: id, reactype, & 
                                                  eleccont, natomspsp, &
                                                  iatomeqsp 
  integer, dimension(:,:), pointer             :: natomconst, ispidrat
  !
  real(kind=drk),    dimension(:), pointer     :: alpha, beta, gamma
  !
  character(len=2),  dimension(:), pointer     :: elementname
  character(len=7),  dimension(:), pointer     :: speciesname
  character(len=16), dimension(:), pointer     :: desc
  !
end type reacnet_type
!
!--------------------------------------------------------------------------------
! --- flags ---
logical                                      :: assdens_flag, metal_flag
!
! --- counter ---
integer                                      :: ncall
!!!$OMP THREADPRIVATE(ncall)
!
! --- ID of species at positions in reaction formulas ---
integer, dimension(:,:), allocatable         :: ircn_ispreac
!
integer, dimension(:,:), allocatable         :: ispecact
!
!
! --- names of chemical species provided in the UIO name entries of the QUC density arrays ---
character(len=7),  dimension(:), allocatable :: cquc_spname
!
!
! --- number of chemical species provided via QUC arrays in the model ---
integer                                      :: nquc_sp
!
! --- indices of chemical species in QUC with regard to species in reaction network ---
integer, dimension(:), allocatable           :: ircn_iquc
integer                                      :: ircn_iqucmetal
!
! --- indices of C and CO refering to provided QUC arrays 
!     (only needed for CO cooling in the radiative transfer) ---
integer                                      :: isp_n, isp_co
integer, dimension(:), allocatable           :: isp_c, nsp_cmult
!
! --- reaction coefficients alpha, beta, gamma ---
type(reacnet_type) :: rn
!
! --- resulting reaction rates ---
real(kind=drk), dimension(:), allocatable   :: rat
!$OMP THREADPRIVATE(rat)
!
!--------------------------------------------------------------------------------------------------
contains

!----------****************------------------------------------------------------------------------
subroutine chem_rn_readchem()
!--------------------------------------------------------------------------------------------------
! NAME:
!   chem_rn_readchem
!
! PURPOSE:
!   Read chemical input data.
!
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! VARIABLES:
!
! MODIFICATION HISTORY:
!   31.01.05 Written by S.Wedemeyer-Boehm (original code by  I. Kamp)
!--------------------------------------------------------------------------------------------------
!                                                                        
use rhd_gl_module
use rhd_box_module
!                                                                       
implicit none 
!
! --- name of input file ---
character(len=80)                                 :: filename
!
! --- index variables ---
integer                                           :: i, j, k
integer                                           :: ios
!                                                                        
! --- auxiliary reals and integers ---
integer                                           :: itmp_id
real(kind=drk)                                    :: rtmp_alpha, rtmp_beta, rtmp_gamma
real(kind=drk), dimension(:), allocatable         :: rtmp_arr
integer, dimension(:), allocatable                :: itmp_arr
integer, dimension(:,:), allocatable              :: itmp2_arr
!
! --- auxiliary strings ---
character(len=1)                                  :: ctmp1
character(len=7), dimension(nreactants+nproducts) :: ctmp_spname
character(len=7), dimension(:), allocatable       :: ctmp7_arr
character(len=16)                                 :: ctmp_text 
character(len=16), dimension(:), allocatable      :: ctmp_arr 
!
! --- local arrays ---
integer,           dimension(:), allocatable, target :: rn_id
real(kind=drk),    dimension(:), allocatable, target :: rn_alpha, rn_beta, rn_gamma
character(len=7),  dimension(:), allocatable, target :: rn_speciesname
character(len=16), dimension(:), allocatable, target :: rn_desc
! 
! --- logical ---
logical                                           :: get 
logical, dimension(nreactants+nproducts)          :: got
logical                                           :: exist_flag
!                                                                       
!--------------------------------------------------------------------------------------------------
# ifdef rhd_krome_debug01
!
do i=1, ncdt_sp
  write(*,'(I2,1X,A7,1X,I6,1X,A60)') i,ccdt_spname(i),icdt_id(i), ccdt_spdesc(i)
end do
!
# endif
!
!--------------------------------------------------------------------------------------------------
! --- format of input data ---
! --- for UMIST database                                                
!1000  FORMAT(1X,I4,1X,4(1A7,1X),A3,1X,A3,1PE8.2,1X,0PF5.2,1X,F8.1,A9)  
! --- for UMIST 99 ratefile                                             
 2000 FORMAT(I4,5(A8,1X),2(1X,A4),1X,1PE8.2,3X,0PF5.2,2X,0PF8.1,A16)                                                
! --- for comparison with FB data my chemistry                          
!1000  FORMAT(1X,I4,1X,4(1A7,1X),A3,1X,A3,1PE8.2,1X,0PF5.2,1X,F12.4,A9) 
!--------------------------------------------------------------------------------------------------
!
! --- number of found species --> will replace nc (read_species) later ---
rn%nspecies=0
!
rn%nreac=0
!
!--------------------------------------------------------------------------------------------------
! --- create file name ---
!
filename=adjustl(trim(par%chem_reacpath))
!
if (filename(1:1).ne.' ') then 
  ctmp1=filename(len_trim(filename):len(filename))
  if (ctmp1(1:1).ne.'/') filename=adjustl(trim(filename))//'/'
end if  
filename=adjustl(trim(filename))//adjustl(trim(par%chem_reacfile))
!
!--------------------------------------------------------------------------------------------------
! --- read reaction table ---
!
write(*,'(A)') '------------------------------------------------------------------------------------'
write(*,'(A,A)') '> Reading reaction table: ',trim(filename)
!   
inquire(file=filename, exist=exist_flag)
if (exist_flag) then 
  !====================================================================================================
  ! ATTENTION! Unit number is always 30 ... Change in later version!
  !====================================================================================================
  ! --- open input file ---
  open(unit=30, file=filename, form='FORMATTED', iostat=ios) 
  !
  if (ios.eq.0) then 
    !
    i=1
    !
    ! --- Read input file line by line, use format of UMIST chemistry ---
10  read(30, 2000, err=11, end=12) itmp_id, (ctmp_spname(j), j=1,nreactants+nproducts), & 
         rtmp_alpha, rtmp_beta, rtmp_gamma, ctmp_text
    !
    get=.true. 
    rn%nreac=rn%nreac+1
    !
    ! --- provide array for reactant/product types ---
    if (rn%nreac == 1) then 
      allocate(ircn_ispreac(rn%nreac,nreactants+nproducts))
      ircn_ispreac(rn%nreac,1:nreactants+nproducts)=-99
    else 
      ! --- make array larger ---
      !
      allocate(itmp2_arr(rn%nreac-1,nreactants+nproducts))  
      itmp2_arr(1:rn%nreac-1,1:nreactants+nproducts)=ircn_ispreac(1:rn%nreac-1,1:nreactants+nproducts)
      deallocate(ircn_ispreac)
      allocate(ircn_ispreac(rn%nreac,nreactants+nproducts))
      ircn_ispreac(1:rn%nreac-1,1:nreactants+nproducts)=itmp2_arr(1:rn%nreac-1,1:nreactants+nproducts)
      ircn_ispreac(rn%nreac,1:nreactants+nproducts)=-99
      deallocate(itmp2_arr)
    end if
    !
    ! --- get type of reactant/product (ircn_ispreac) ---
    !     (loop over reactants and products)
    !
    do k=1,nreactants+nproducts 
      got(k) = .FALSE. 
      !                                                                 
      if (ctmp_spname(k) == '      ') then
        ircn_ispreac(rn%nreac,k)=-99
        got(k) = .TRUE. 
      else  
        ! --- UMIST chemistry ---
        !
        do i=1,ncdt_sp
          if (ctmp_spname(k) == ccdt_spname(i)) then 
            ircn_ispreac(rn%nreac,k)=icdt_id(i)
            got(k) = .TRUE. 
          end if
        end do
        !
        if (.not.got(k)) then 
          ! --- Search for chemical species ---
          if (rn%nspecies.gt.0) then 
            do j=1,rn%nspecies 
              if (ctmp_spname(k) == rn_speciesname(j)) then 
                ! --- species is already known ---
                ircn_ispreac(rn%nreac,k)=j 
                got(k)=.TRUE. 
              end if
            end do
          end if
          !
          if (.not.got(k)) then 
            ! --- found new species (rn_speciesname) ---
            !write(6,'(A,A)') '> Found new species: ', ctmp_spname(k)
            rn%nspecies=rn%nspecies+1
            !
            if (rn%nspecies == 1) then 
              allocate(rn_speciesname(1))
              rn_speciesname(1)=ctmp_spname(k)
            else
              allocate(ctmp7_arr(rn%nspecies-1)) 
              ctmp7_arr(1:rn%nspecies-1)=rn_speciesname(1:rn%nspecies-1)
              deallocate(rn_speciesname)
              !
              allocate(rn_speciesname(rn%nspecies))
              rn_speciesname(1:rn%nspecies-1)=ctmp7_arr(1:rn%nspecies-1)
              rn_speciesname(rn%nspecies)=ctmp_spname(k)
              deallocate(ctmp7_arr)
            end if
            !
            ircn_ispreac(rn%nreac,k)=rn%nspecies
            got(k)=.TRUE. 
          end if
          !
        end if
        !
      end if
      get=get.and.got(k) 
      !write(*,*) k,'>',ctmp_spname(k),'<','>',ctmp7,'<',ilen,' => ircn_ispreac=',ircn_ispreac(rn%nreac,k)
    end do ! k
    !                                                                        
    ! --- get reaction coefficients, IDs, and text ---
    if (get) then 
      !
      if (rn%nreac.le.1) then 
        ! --- first array elements ---
        !
        ! --- reaction coefficients ---
        !
        allocate(rn_alpha(1))
        rn_alpha(1)=rtmp_alpha
        !     
        allocate(rn_beta(1))
        rn_beta(1)=rtmp_beta
        !     
        allocate(rn_gamma(1))
        rn_gamma(1)=rtmp_gamma
        !     
        ! --- reaction IDs ---
        ! 
        allocate(rn_id(1))
        rn_id(1)=itmp_id
        !     
        ! --- reaction text ---
        ! 
        allocate(rn_desc(1))
        rn_desc(1)=ctmp_text
        !     
      else
        ! --- store data in temporary array, destroy old arrays, store data  
        !     in new larger (+1) array ---
        !
        ! --- reaction coefficients ---
        !
        allocate(rtmp_arr(rn%nreac-1))
        ! 
        rtmp_arr(1:rn%nreac-1)=rn_alpha(1:rn%nreac-1)
        deallocate(rn_alpha)
        allocate(rn_alpha(rn%nreac))
        rn_alpha(1: rn%nreac-1)=rtmp_arr(1:rn%nreac-1)
        rn_alpha(rn%nreac)=rtmp_alpha
        !
        rtmp_arr(1:rn%nreac-1)=rn_beta(1:rn%nreac-1)
        deallocate(rn_beta)
        allocate(rn_beta(rn%nreac))
        rn_beta(1: rn%nreac-1)=rtmp_arr(1:rn%nreac-1)
        rn_beta(rn%nreac)=rtmp_beta
        !
        rtmp_arr(1:rn%nreac-1)=rn_gamma(1:rn%nreac-1)
        deallocate(rn_gamma)
        allocate(rn_gamma(rn%nreac))
        rn_gamma(1: rn%nreac-1)=rtmp_arr(1:rn%nreac-1)
        rn_gamma(rn%nreac)=rtmp_gamma
        !
        deallocate(rtmp_arr)
        !
        ! --- reaction IDs ---
        ! 
        allocate(itmp_arr(rn%nreac-1))
        !
        itmp_arr(1:rn%nreac-1)=rn_id(1:rn%nreac-1)
        deallocate(rn_id)
        allocate(rn_id(rn%nreac))
        rn_id(1: rn%nreac-1)=itmp_arr(1:rn%nreac-1)
        rn_id(rn%nreac)=itmp_id
        !    
        deallocate(itmp_arr)
        !
        ! --- reaction text ---
        ! 
        allocate(ctmp_arr(rn%nreac-1))
        !
        ctmp_arr(1:rn%nreac-1)=rn_desc(1:rn%nreac-1)
        deallocate(rn_desc)
        allocate(rn_desc(rn%nreac))
        rn_desc(1: rn%nreac-1)=ctmp_arr(1:rn%nreac-1)
        rn_desc(rn%nreac)=ctmp_text
        !     
        deallocate(ctmp_arr)
        !
      end if
      !
    end if
    !
    ! --- Finished line, read next line ---
    i=i+1
    goto 10 
    !
    ! --- Error occured ---
11   write(6,'(A)') 'KROME> ERROR in chemistry input file chem.dat!' 
    !ierr=101
    stop 
    !
    ! --- End of file encountered ---
12  close(30) 
    !
  else 
    write(6,'(A)') 'KROME> ERROR while opening chemistry input file!'
    stop
  end if
else 
  write(6,'(A)') 'KROME> ERROR: Chemistry input file not found!'
  stop
end if
!
! --- allocate pointers in structure RN ---
allocate(rn%speciesname(rn%nspecies))
rn%speciesname(1:rn%nspecies)=rn_speciesname(1:rn%nspecies)
allocate(rn%alpha(rn%nreac))
rn%alpha(1:rn%nreac)=rn_alpha(1:rn%nreac)
allocate(rn%beta(rn%nreac))
rn%beta(1:rn%nreac)=rn_beta(1:rn%nreac)
allocate(rn%gamma(rn%nreac))
rn%gamma(1:rn%nreac)=rn_gamma(1:rn%nreac)
allocate(rn%id(rn%nreac))
rn%id(1:rn%nreac)=rn_id(1:rn%nreac)
allocate(rn%desc(rn%nreac))
rn%desc(1:rn%nreac)=rn_desc(1:rn%nreac)
!
! --- deallocate temporary arrays ---
deallocate(rn_speciesname)
deallocate(rn_alpha)
deallocate(rn_beta)
deallocate(rn_gamma)
deallocate(rn_id)
deallocate(rn_desc)
!
end subroutine chem_rn_readchem


!----------**********************------------------------------------------------------------------
subroutine chem_rn_analysespecies()
!--------------------------------------------------------------------------------------------------
! NAME:
!   chem_rn_module ('chem_rn_analysespecies')
!
! PURPOSE:
!   Analyse chemical input data.
!
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! VARIABLES:
!   rn%elementname       : (character) atomic contiuents (chemical elements) 1..rn%nelements
!   two      : true if the element symbol is composed of two letters (e.g., He)
!   rn%natomspsp       : total number of atomic constituents in each molecule/species
!            (needed for total particle densities, e.g. rn%natomspsp=2 for CO)
!   rn%natomconst       : atomic constituents for each molecule/species
!   rn%eleccont       : charge of molecule/species (needed for electron equation)
!   rn%nelements       : number of  atomic constituents (elements)
!   rn%speciesname       : molecular species (1..rn%nspecies)
!   rn%nspecies       : number of molecular species
!   rn%iatomeqsp      : index of species which is identical with a single element 
!
! MODIFICATION HISTORY:
!   31.01.05 Written by S.Wedemeyer-Boehm (original code by  I. Kamp)
!--------------------------------------------------------------------------------------------------
!                                                                        
implicit none
!                                                                       
integer                                 :: i, j, k, l
integer, dimension(:), allocatable      :: nreac, nprod
integer                                 :: nlen
character(len=80)                       :: ctmp80
character(len=7)                        :: ctmp7
character*(2), dimension(:), allocatable:: ctmp2_arr
character*(2), &
         dimension(nreactants+nproducts):: ctlet
character(len=2)                        :: ctmpsym 
character*(1)                           :: clet(2)
character*(1), &
         dimension(nreactants+nproducts):: cslet
logical                                 :: two, get, known_flag 
logical                                 :: existflag
!                                                                       
! --- local arrays ---
character(len=2),  dimension(:), allocatable, target :: rn_elementname
!--------------------------------------------------------------------------------
! --- output ---
 2001 FORMAT(I4,1X,I4,1X,A65,1X,1PE8.2,3X,0PF5.2,2X,0PF8.1,A16)
!--------------------------------------------------------------------------------
!
! --- text output ---
!
write(6,'(A)')         ' '
write(6,'(A6,I3,A11)') 'Found ', rn%nreac, ' reactions.'
write(6,'(A)')         ' '
write(6,'(A)') ' nr   ID  reactants                          products                       alpha       beta     gamma  ref.'
write(*,'(A)') '------------------------------------------------------------------------------------------------------------'
!
do l=1,rn%nreac
  ! --- construct string for reaction equation ---
  ctmp80(:)=' '
  do k=1,nreactants+nproducts
    if (ircn_ispreac(l,k).gt.0) then  
      ! --- true chemical species ---
      ctmp7=rn%speciesname(ircn_ispreac(l,k))
    else 
      if (ircn_ispreac(l,k) == -99) then 
        ! --- blank field ---
        ctmp7(:)=' '
      else 
        ! --- search for special species/processes ---
        do i=1,ncdt_sp
          if (ircn_ispreac(l,k) == icdt_id(i))  ctmp7(:)=ccdt_spname(i)
        end do
      end if
    end if
    ! 
    ! --- start position ---
    i=(k-1)*(7+3)+1
    ! --- add one to account for the longer right-arrow ---
    if (k.gt.nreactants) i=i+1
    ! 
    ! --- paste string ---
    ctmp80(i:i+6)=ctmp7
    !
    ! --- add plus sign ---
    if ((len_trim(ctmp7).gt.0).and.(k.gt.1).and.(k.ne.nreactants+1)) then 
      ! --- start position ---
      i=(k-1)*(7+3)-2
      if (k.gt.nreactants) i=i+1
      ctmp80(i:i+2)=' + '
    end if
    !
  end do
  ctmp80(nreactants*(7+3)-3+1:nreactants*(7+3)+1)= ' -> '
  !
  write(*,2001) l, rn%id(l), ctmp80, rn%alpha(l), rn%beta(l), rn%gamma(l), rn%desc(l)              
end do
write(*,'(A)') & 
  '------------------------------------------------------------------------------------------------------------'   
!
! --- write out symbols for species ---
!
write(*,'(A6,I3,A9)') 'Found ', rn%nspecies, ' species:'
do i=1,rn%nspecies
  write(6,'(I3,2X,A8)') i, rn%speciesname(i) 
end do
!
!--------------------------------------------------------------------------------
! --- Analyse input ---
!
! --- find atomic constituents --
rn%nelements=0
!
do j=1,rn%nspecies  
  ctmp7=adjustl(rn%speciesname(j))
  nlen=len_trim(ctmp7)  
  !
  if (nlen.gt.0) then 
    do i=1,nlen
      ! 
      ctmpsym='  '
      existflag=.FALSE.
      !
      if ((lge(ctmp7(i:i),'A')).and.(lle(ctmp7(i:i),'Z'))) then 
        ctmpsym(1:1)=ctmp7(i:i)
        !
        ! --- check subsequent character if available ---
        if (i.lt.nlen) then 
          if ((lge(ctmp7(i+1:i+1),'a')).and.(lle(ctmp7(i+1:i+1),'z'))) then 
            ctmpsym(2:2)=ctmp7(i+1:i+1)
          end if      
        end if
      end if
      !
      !write(6,*) '// ',i,'>',ctmpsym,'<'
      !
      if (ctmpsym(1:1).ne.' ') then 
        !        
        ! --- check if atom already exists ---
        if (rn%nelements.gt.0) then 
          do k=1,rn%nelements 
            if (ctmpsym == rn_elementname(k)) then 
                existflag=.TRUE.
                ! --- atom already known ---
                !write(*,*) '                       -> atom already known: ',rn_elementname(k)
            end if 
          end do
        end if
        !
        if (.not.existflag) then 
          !write(*,*) '                       -> new atom: >',ctmpsym,'<'
          rn%nelements=rn%nelements+1
          if (rn%nelements == 1) then 
            allocate(rn_elementname(1))
          else 
            allocate(ctmp2_arr(rn%nelements-1))
            ctmp2_arr(1:rn%nelements-1)=rn_elementname(1:rn%nelements-1)
            deallocate(rn_elementname)
            allocate(rn_elementname(rn%nelements))
            rn_elementname(1:rn%nelements-1)=ctmp2_arr(1:rn%nelements-1)
            deallocate(ctmp2_arr)
          end if
          !
          ! --- add new element entry ---
          rn_elementname(rn%nelements)=ctmpsym
        end if
        !
      end if
      ! 
    end do
    ! 
  end if 
end do
!
! --- allocate pointers in structure RN ---
allocate(rn%elementname(rn%nelements))
rn%elementname(1:rn%nelements)=rn_elementname(1:rn%nelements)
!
! --- deallocate temporary arrays ---
deallocate(rn_elementname)
!
! --- text output ---
write(*,'(A20)') '----------------------------------------'
write(*,'(A11,I3,A10)') 'Composed of', rn%nelements, ' elements:'                           
do i=1,rn%nelements
  write(6,'(I3,2X,A2)') i, rn%elementname(i) 
end do
!
! --- loop over all reactions ---
write(*,'(A)') '----------------------'
write(*,'(A)') 'Non-default reactions:'
write(*,'(A)') ' nr  reactype  description'
!
! --- initialise arrays
allocate(nreac(rn%nreac), nprod(rn%nreac))
nreac(1:rn%nreac)=0 
nprod(1:rn%nreac)=0 
rn%reactype(1:rn%nreac)=0.
!
do l=1,rn%nreac
  ! --- number of reactants and products for each reaction ---
  !nreac(l)=0 
  !nprod(l)=0 
  !rn%reactype(l)=0.
  !
  ! --- find number of reactants and products in the chemical reaction 
  ! --- loop over all reactants/products (index j) for each reaction (index l)
  do j=1,nreactants+nproducts 
    ! --- type of reactant/product
    i=ircn_ispreac(l,j) 
    !
    ! --- count number of reactants and products
    if (i.ge.0) then 
      if (j.le.nreactants) then 
        nreac(l)=nreac(l)+1
      else 
        nprod(l)=nprod(l)+1
      endif
    endif
    !
    ! --- type of reaction 
    if ((i.lt.0).and.(i.ne.-99).and.(j.le.nreactants)) then 
      rn%reactype(l)=i 
    endif 
  end do
  !                                                                       
  ! --- set upper limit for products loop (first NREACTANTS are always reactant)
  nprod(l)=nreactants+nprod(l) 
  !
  ! --- text output ---
  ctmp80(:)=' '
  do j=1, ncdt_sp
    if (icdt_id(j) ==  rn%reactype(l)) then 
      ctmp80(1:len_trim(ccdt_spdesc(j)))=ccdt_spdesc(j)
    end if
  end do
  !
  if (rn%reactype(l).ne.0) write(*,'(I3,2X,I8,2X,A)') l,rn%reactype(l), ctmp80
  !
end do ! l
!                                                                        
!--------------------------------------------------------------------------------------------------
! --- check if chemical symbols are known ---
!
do i=1,rn%nelements 
  known_flag=.FALSE.
  do j=1,nchsym
    !
    if (rn%elementname(i) == cchsym(j)) then 
      known_flag=.TRUE.
    end if
  end do 
  !
  if (.not.known_flag) then 
    write(6,'(A38,A2)') '> WARNING! Chemical symbol not known: ',rn%elementname(i)
  end if
  !
end do
!
!--------------------------------------------------------------------------------
! --- initialise vector which contains the electron contribution of each molecule                                                          
!
allocate(rn%eleccont(rn%nspecies), rn%natomspsp(rn%nspecies))
allocate(rn%iatomeqsp(rn%nelements))
!
rn%eleccont(:)=0 
rn%natomspsp(:)=0 
!
! --- initialise matrix which contains for each element in which molecule
!     it can be found (first index: 1,nspecies, second: 1,nelements)                                               
allocate(rn%natomconst(rn%nspecies,rn%nelements))
rn%natomconst(:,:)=0 
!                                                                       
! --- initialise array rn%iatomeqsp (1:nelements)
rn%iatomeqsp(:)=0
!                                                                       
!--------------------------------------------------------------------------------
! --- find the atomic constituents for each molecule --> RN%NATOMCONST ---
!     loop over all molecules                                           
do j=1,rn%nspecies 
  rn%eleccont(j)=0
  !                                                         
  do k=1,7 
    cslet(k) = rn%speciesname(j)(k:k)
    if (k.lt.7) ctlet(k) = rn%speciesname(j)(k:k+1)
    !
    if (rn%speciesname(j)(k:k) == '+') then 
      rn%eleccont(j)=rn%eleccont(j) + 1 
    else if (rn%speciesname(j)(k:k) == '-') then 
      rn%eleccont(j)=rn%eleccont(j)-1 
    end if
  end do! k
  !
  !write(6,*) 'Electron contribution for ',rn%speciesname(j),': ',rn%eleccont(j)
  !
  !-----------------------------------------------
  ! --- loop over all elements ---
  do i=1,rn%nelements 
    clet(1)=rn%elementname(i)(1:1) 
    clet(2)=rn%elementname(i)(2:2) 
    two = .FALSE. 
    !                                                                      
    ! --- decide wether it is an element with two letters ---
    if (clet(2).ne.' ') two = .TRUE. 
    !                                                            
    if (two) then 
      do k=1,6 
        if (ctlet(k) == rn%elementname(i)) then 
          if (k.lt.6) then 
            ! --- decide wether it is multiple present ---
            if ((ichar(cslet(k+2)).le.ichar('9')).and. & 
                (ichar(cslet(k+2)).ge.ichar('1'))) then                                                  
              rn%natomconst(j,i) = ichar(cslet(k+2)) - 48 
            else 
              rn%natomconst(j,i) = 1 
            end if
          else 
            rn%natomconst(j,i) = 1 
          end if
        end if
      end do ! k
    else 
      do k=1,7 
        if (cslet(k) == clet(1)) then 
          ! --- decide wether cslet(k) is only part of an element with two letters ---
          if (k.lt.7) then 
            get=.TRUE. 
            do l=1,rn%nelements 
              if ((l.ne.i).and.(ctlet(k) == rn%elementname(l))) get=.FALSE.
            end do
          end if
          !                                                    
          if ((k.lt.7).and.(get)) then 
            ! --- decide wether it is multiple present ---
            if ((ichar(cslet(k+1)).le.ichar('9')).and. &
                (ichar(cslet(k+1)).ge.ichar('1'))) then                                                  
              rn%natomconst(j,i) = ichar(cslet(k+1)) - 48 
            else 
              rn%natomconst(j,i) = 1 
            end if
          else if (get) then 
            rn%natomconst(j,i) = 1 
          end if
        end if
      end do ! k
    end if
  end do ! i
end do ! j
!
!================================================================================                 
! --- get total number of atomic constituents in each molecule --> RN%NATOMSPSP 
!     for total particle densities ---  
do j=1,rn%nspecies 
  do i=1,rn%nelements 
    rn%natomspsp(j)=rn%natomspsp(j)+rn%natomconst(j,i) 
  end do
end do
!
!================================================================================                 
! --- find indices of atomic constituents in the field of molecules --> 
!     (for particle conservation) ---
! --- loop over all molecules ---
!
do j=1,rn%nspecies
  ! --- split up the molecule in single letters and groups of two letter ---
  do k=1,7 
    cslet(k)=rn%speciesname(j)(k:k) 
    if (k.lt.7) then 
      ctlet(k) = rn%speciesname(j)(k:k+1) 
    end if 
  end do 
  !                                                                   
  ! --- loop over all elements ---
  do i=1,rn%nelements 
    clet(1) = rn%elementname(i)(1:1) 
    clet(2) = rn%elementname(i)(2:2) 
    two = .FALSE. 
    !                                                                        
    ! --- decide wether it is an element with two letters ---
    if (clet(2).ne.' ') two = .TRUE. 
    !
    if (two) then 
      if (ctlet(1) == rn%elementname(I)) then 
        if (cslet(3) == ' ') then 
          rn%iatomeqsp(i) = j 
        end if 
      end if 
    else 
      if (cslet(1) == clet(1)) then 
        if (cslet(2) == ' ') then 
          rn%iatomeqsp(i) = j 
        end if 
      end if 
    end if 
  end do 
end do 
!                                                                        
end subroutine chem_rn_analysespecies


!----------****************------------------------------------------------------------------------
subroutine chem_rn_getrates(rmetal, rT, outstr, ierr) 
!--------------------------------------------------------------------------------------------------
! NAME:
!    chem_rn_init ('chem_rn_getrates')
!
! PURPOSE:
!   Calculate reaction rates 
!
! CATEGORY:
!
! CALLING SEQUENCE:
!   call chem_rn_getrates(rmetal, rT, outstr, ierr)
!
! INPUT:
!
! OUTPUT:
!   None
!
! ROUTINES:
!   None
!
! MODULES:
!
! SIDE EFFECTS:
!
! TODO:
!   If statement inside loop might decrease performance. Optimisation!?
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2004-02-09 Written by Sven Wedemeyer, KIS Freiburg
!   2021-02-15 (S.W.) Check for negative rates in debug mode 
!--------------------------------------------------------------------------------------------------
!
implicit none
!                                                                       
! --- I/O ---
!real(kind=drk), intent(in)             :: g0, cri, av, w
real,           intent(in)             :: rmetal, rT
integer,        optional, intent(out)  :: ierr
character*(*),  optional, intent(out)  :: outstr
!
! --- constant quantities for photoreaction rates ---
real(kind=drk), parameter                  :: g0  = 1.D-5  ,&
                                              cri = 1.D-17 ,&
                                              av  = 1.D-5  ,&
                                              w   = 0.6D0  ! --- UV ---
! --- local variables ---
integer                                :: l
real(kind=drk)                         :: T300, metal, T
!
!integer :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
!
!--------------------------------------------------------------------------------------------------
!
rat(:)=0.D0
!
metal=dble(rmetal)
T=dble(rT)
!
! --- temperature in units of 300K ---
T300 = T / 300.0D0 
!
! --- loop over all reactions                                           
do l=1, rn%nreac 
  !
  select case(rn%reactype(l))
    case(0)
    ! --- basic rate for chemical reactions --> R                     
    rat(l) = rn%alpha(l) * T300**rn%beta(l) * DEXP(-rn%gamma(l)/T) 
    !
    case(-1)
      ! --- cosmic-ray ionisation rate --> R                            
      rat(l) = cri * rn%alpha(l) 
      !
    case(-2)
      ! --- cosmic-ray-induced photoreaction rate --> R following UMIST 
       rat(l) = cri * rn%alpha(l) * T300**rn%beta(l) * rn%gamma(l) / (1.D0 - w) 
      !
    case(-3)
       ! --- interstellar photoreaction rate 1 --> R following UMIST nomenclature                                
       rat(l) = g0 * rn%alpha(l) * DEXP( - rn%gamma(l) * av) 
       !
    case(-4)
       ! --- interstellar photoreaction rate 2 --> R
       rat(l) = g0 * rn%alpha(l) * DEXP( - av * (rn%beta(l) - rn%gamma(l) * av) ) 
       !
    case(-5,-25)
       ! --- special for some reactions (e.g. radiative association) -->
       rat(l) = rn%alpha(l) * T300**rn%beta(l) * DEXP( - rn%gamma(l) / T) 
       !
    case(-35)
       ! --- catalytic reaction --> R                                    
       rat(l) = metal * rn%alpha(l) * T300**rn%beta(l) * DEXP(-rn%gamma(l)/T) 
       !
    case default
       !
       write(6,'(A,I5,A,I4,A)') 'KROME> ERROR! Unknown reaction type', & 
                                rn%reactype(l),' for reaction ',l,'!'
       ierr=121
       write(outstr,'(A,I5,A,I4,A)') 'KROME> ERROR! Unknown reaction type', & 
                                     rn%reactype(l),' for reaction ',l,'!'
       !stop
       !
  end select
  !
  ! --- check resulting rates                                             
  !write(*,*) "#",omp_get_thread_num(), 'L=',l,'|IRT=',rn%reactype(l),'|RAT=',rat(l)
  !
# ifdef rhd_krome_debug01
  if (isnan(real(rat(l))) .eq. .true.) then
      write(*,*) 'WARNING. Rate value is NaN. Forced to zero. Rate #',l
      rat(l)=0.D0
  end if
# endif

end do                                                                      
!  
end subroutine chem_rn_getrates


!----------**********************------------------------------------------------------------------
subroutine chem_rn_assign_density()
!--------------------------------------------------------------------------------------------------
! NAME:
!    chem_rn_assign_density ('chem_rn_assign_density')
!
! PURPOSE:
!   Check if number densities for the required chemical species are available. 
!   Assign indices of these arrays to indices of chemical species in the reaction network.
!
! CATEGORY:
!
! CALLING SEQUENCE:
!   call chem_rn_assign_density
!
! INPUT:
!
! OUTPUT:
!   None
!
! ROUTINES:
!   None
!
! MODULES:
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
!   25.01.05 Written by Sven Wedemeyer-Boehm, KIS Freiburg
!--------------------------------------------------------------------------------------------------
!
use rhd_gl_module
use rhd_prop_module
use rhd_box_module
!                                                                       
integer                                 :: i, j, k, nlen
logical                                 :: exist_flag
character(len=80)                       :: ctmp
!
!--------------------------------------------------------------------------------------------------
!
metal_flag=.FALSE.
!
!prop%nquc=size(model%quc,4)
!
!write(*,*) 'NQUC = ',prop%nquc
nquc_sp=prop%nquc
!
allocate(cquc_spname(nquc_sp))
allocate(ircn_iquc(nquc_sp))
allocate(isp_c(nquc_sp), nsp_cmult(nquc_sp))
!
! --- initialise index array ---
ircn_iquc(1:nquc_sp)=-99
!
! --- check if number densities for the required chemical species are available ---
write(6,'(A)')    '--- KROME_ASSIGN_DENSITY ---'
write(6,'(A,I3)') 'Number of available density arrays: ',prop%nquc
!
do i=1,prop%nquc
  write(*,'(I3,2X,A6,2X,A50)') i, prop%quc_ident(i), prop%quc_name(i) 
end do
!
do i=1,prop%nquc
  exist_flag=.FALSE.
  ! --- extract name of chemical species from UIO data header (quc_name) --- 
  !     (take characters after "Number density of"
  ctmp=adjustl(prop%quc_name(i))
  nlen=len_trim(ctmp)
  !
  if (nlen.gt.len_trim("Number density of")) then 
    ctmp=adjustl(ctmp(len_trim("Number density of")+1:len(ctmp)))
    cquc_spname(i)=ctmp(1:7)
    !write(6,'(A)') cquc_spname(i)
    !
    nlen=len_trim(ctmp)
    !
    do j=1,rn%nspecies
      if (ctmp(1:nlen) == rn%speciesname(j)) then 
        ! --- found chemical species ---
        !write(*,*) 'MATCH : ',ctmp(1:nlen),rn%speciesname(j),i
        exist_flag=.TRUE.
        ircn_iquc(j)=i
      end if
    end do
    !
    ! --- check if the species is of special type ---
    do k=1,ncdt_sp 
      if (ctmp(1:nlen) == ccdt_spname(j)) then 
        ! --- found chemical species of special type  ---
        exist_flag=.TRUE.
      end if             
    end do 
    ! --- more cases ---
    if ((ctmp(1:nlen) == 'metal').or.(ctmp(1:nlen) == 'M')) then 
      ! --- found array for representative metal  ---
      write(*,'(A)') 'KROME> Found array for representative metal. Ignoring metal abundance parameter (if set).'
      exist_flag=.TRUE.
      metal_flag=.TRUE.
      ircn_iqucmetal=i
    end if
    ! 
    if (.not.exist_flag) then 
        write(6,'(A32,A7,A37)') '> WARNING! The chemical species ',ctmp, & 
                                ' is not used in the reaction network.'
        write(6,'(A21,A8,A23)') '           The array ', prop%quc_ident(i), & 
                                ' will be advected only!' 
    end if
    !
  else 
     write(*,'(A34,A7,A37)') 'KROME> ERROR! UIO name entry for array ',prop%quc_ident(i), & 
                             ' starts not with "Number density of":'
     write(*,'(9X,A60)')     prop%quc_name(i)
     stop
  end if
  !
end do
!
! --- check if all required density arrays are provided ---
!
do j=1,rn%nspecies
  exist_flag=.FALSE.
  do i=1,prop%nquc
    if (rn%speciesname(j) == cquc_spname(i)) then 
      exist_flag=.TRUE.
    end if
  end do
  !
  if (.not.exist_flag) then 
     write(6,'(A44,A7,A12)') 'KROME> ERROR! Density array for chemical species ', & 
                             rn%speciesname(j),' is missing!'
     stop
  end if 
end do
!
if (.not.metal_flag) then 
  write(*,'(A)') 'KROME> No array for representative metal found.'
  write(*,'(A)') '         Taking metal from abundance parameter.'
  write(*,'(A,E12.2)') '         Metal abundance: ',par%chem_abumetal 
end if
!
!return
end subroutine chem_rn_assign_density


!----------************----------------------------------------------------------------------------
subroutine chem_rn_init()
!--------------------------------------------------------------------------------------------------
! NAME:
!    chem_rn_init ('chem_rn_init')
!
! PURPOSE:
!   Read chemistry data, initialise quantities.
!
! CATEGORY:
!
! CALLING SEQUENCE:
!   call chem_rn_init
!
! INPUT:
!
! OUTPUT:
!   None
!
! ROUTINES:
!   None
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! FUTURE IMPROVEMENTS:
!   create matrix ispecact
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   09.02.04 Written by Sven Wedemeyer-Boehm, KIS Freiburg
!--------------------------------------------------------------------------------------------------
!
implicit none
!                                                                       
integer                                 :: i, j
!
!# ifdef rhd_chem_tune01
!common / CHEMREAC / ispecact 
!# endif
!
character(len=7), dimension(7) :: ctmp
!
!--------------------------------------------------------------------------------------------------
! --- init ---
rn%nspecies=0
!
! --- density arrays not assigned yet ---
assdens_flag=.FALSE.
!
! --- Read in the basic chemical data ---
call chem_rn_readchem()
!
! --- now allocate arrays ---
!     (rn%nspecies, rn%nreac are now known)
!
allocate(rn%reactype(rn%nreac))
!
!integer, dimension(nspecies)            :: ircn_iquc
!integer, dimension(nspecies)            :: isp_c, nsp_cmult
!
call chem_rn_analysespecies()
!
!# ifdef rhd_chem_tune01
! --- create matrix ispecact 
!     ispecact: first index: species, second index: reaction
!     multiplicator for construction of derivative
!     +1: product, -1: educt(reactant), 0: not involved 
!
allocate(ispecact(rn%nspecies,rn%nreac))
ispecact(1:rn%nspecies,1:rn%nreac)=0
!
do i=1,rn%nreac 
  ! --- reactants
  do j=1,nreactants
    if (ircn_ispreac(i,j).gt.0) then 
      ispecact(ircn_ispreac(i,j),i)=ispecact(ircn_ispreac(i,j),i)-1 
    endif
  end do
  ! --- products 
  do j=nreactants+1,nreactants+nproducts
    if (ircn_ispreac(i,j).gt.0) then 
      ispecact(ircn_ispreac(i,j),i)=ispecact(ircn_ispreac(i,j),i)+1 
    endif
  end do
end do 
!# endif
!
allocate(rn%ispidrat(rn%nreac,nreactants))
rn%ispidrat(1:rn%nreac,1:nreactants)=rn%nspecies+1
do i=1,rn%nreac
  do j = 1, nreactants
    if (ircn_ispreac(i,j) .ge. 1) then 
      rn%ispidrat(i,j)=ircn_ispreac(i,j)
    end if 
  end do
end do
!
end subroutine chem_rn_init


!----------******************----------------------------------------------------------------------
subroutine krome_rn_sourcestep(model, action, dtime, outstr,ierr)
!--------------------------------------------------------------------------------------------------
! NAME:
!   krome_rn_sourcestep ('krome_rn_source_step')
!
! PURPOSE:
!   Update model by applying CO source term step.
!   
! CATEGORY:
!   Hydrodynamics
!
! CALLING SEQUENCE:
!   call chem_rn_sourcestep(model, action, dtime, outstr=outstr,ierr=ierr)
!
! INPUT:
!   action: ('action') derived type with collection of control parameters
!   dtime:  ('delta_time') real, time step to be used
!
! INPUT/OUTPUT:
!   model:  ('model') derived type, contains all model data (e.g. model%rho)
!     m1:        ('m_number_1') integer, 1. dimension, start index
!     n1:        ('number_1')   integer, 1. dimension, end index
!     m2:        ('m_number_2') integer, 2. dimension, start index
!     n2:        ('number_2')   integer, 2. dimension, end index
!     m3:        ('m_number_3') integer, 3. dimension, start index
!     n3:        ('number_3')   integer, 3. dimension, end index
!     nghost1:   ('number_ghost_1') integer, 1. dimension, number of ghost cells
!     nghost2:   ('number_ghost_2') integer, 2. dimension, number of ghost cells
!     nghost3:   ('number_ghost_3') integer, 3. dimension, number of ghost cells
!     itime:     ('index_time') integer, time step number
!     nquc:      ('number_quantities_center') integer, number of additional fields
!     time:      ('time') real, time
!     xc1:       ('x_center_1') real, 1. coordinate of cell centers
!     xc2:       ('x_center_2') real, 2. coordinate of cell centers
!     xc3:       ('x_center_3') real, 3. coordinate of cell centers
!     xb1:       ('x_boundary_1') real, 1. coordinate of cell boundaries
!     xb2:       ('x_boundary_2') real, 2. coordinate of cell boundaries
!     xb3:       ('x_boundary_3') real, 3. coordinate of cell boundaries
!     rho:       ('rho') real, dimension(:,:,:), density
!     ei:        ('energy_internal') real, dimension(:,:,:), internal energy density per mass
!     v1:        ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!     v2:        ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!     v3:        ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!     quc:       ('quantity_central') real, dimension(:,:,:,:), cell centered set of quantities
!                                     e.g. quc(i1,i2,i3,iquc)
!
! OUTPUT:
!   outstr: ('output_string') character, optional
!   ierr:   ('index_error') integer, optional error number:
!           =0: normal run
!           1-99:  "small" error (e.g. time step limit violation)
!           >=100: "severe" error:
!
! VARIABLES:
!   m1:          ('m_1') integer, lower index, 1st dimension
!   n1:          ('n_1') integer, upper index, 1st dimension
!   m2:          ('m_2') integer, lower index, 2nd dimension
!   n2:          ('n_2') integer, upper index, 2nd dimension
!   m3:          ('m_3') integer, lower index, 3rd dimension
!   n3:          ('n_3') integer, upper index, 3rd dimension
!
!   Quantities which are kept constant for the present application: 
!   (The quantities are only important for photoreactions but there are at the 
!    moment no photoreactions in the database.)
!    g0  : integrated UV radiation field between 912 and 1110 A (check boundaries
!          with UMIST paper);
!          here:  g0   = 1.D-5
!    cri : cosmic rays (typical value, but most probably on the lower edge); 
!          here:  cri  = 1.D-17
!    av  : visual extinction (only necessary for photoreactions and CR reactions)
!          here:  av   = 1.D-5
!    wuv : grain albedo in the far UV, typically 0.6 at 150 nm (Millar et al. 1997)
!          (only necessary for CR induced photoreactions
!          here:  wuv  = 0.6D0
!
!   rpar, ipar: user-defined scalars(/arrays) for solver input, here unused
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! WARNING:
!   Parallelisation is troublesome on some architectures (problems with DVODE).
!   Maybe liw and lrw must be declared shared. Or put all in one large parallel 
!   region. UPDATE: Issue was solved in January 2021. 
!
! EXAMPLE:
!
! TODO:
!   move all if statements outside the loops, use vector masks instead
!
! MODIFICATION HISTORY:
!   09.02.04 (SW) First version
!   03.02.05 (SW) Revised version
!   31.08.06 (SW) minor changes
!   2008-08-12 (M.S., B.F.) eosinter_T gets additional arguments (m1,n1,...)
!   2019-11-24 (M.S., B.F.) eosinter_TOMP gets additional arguments (nghost1,...)
!   2021-04-30 (SW, MSz) Usage of updated DVODE package, OpenMP improvements  
!--------------------------------------------------------------------------------------------------
!
use tabinter_module
use gasinter_module
use rhd_gl_module
use rhd_action_module
use rhd_box_module
!
use krome_main
use krome_user
!
# ifdef rhd_krome_t01
use timing_module
# endif
!
implicit none
!
! --- I/O parameters ---
type(box_type),             intent(inout)  :: model
type(action_type),          intent(in)     :: action
real,                       intent(in)     :: dtime
integer,          optional, intent(out)    :: ierr
character(len=*), optional, intent(out)    :: outstr

!
! --- Local variables ---
integer                                    :: i
integer                                    :: m1, n1, m2, n2, m3, n3, i1, i2, i3
!
! --- error handling ---
integer                                    :: ier_flag, ierr0, ier
character(len=80)                          :: outstr0
!
! --- physical quantities ---
real,           dimension(:,:,:), allocatable :: T, metal
real(kind=drk), dimension(:), allocatable     :: y, rates
! KROME temperature
! Note that because KROME can change the temperature inside itself based on
! heating/cooling, this approach will fail when that is included.
! Instead, we should compile KROME in single precision!
real*8 :: Tgas
!
!--- additional variables for vode ---
integer                                    :: nspec, lrw, liw, &
                                              itol, itask, istate, iopt, mf
integer                                    :: ipar, ng
integer, dimension(:), allocatable         :: iwork
real(kind=drk)                             :: dbltime0, dbltime1
real(kind=drk)                             :: rtol, atol
real(kind=drk)                             :: rpar
real(kind=drk), dimension(:), allocatable  :: rwork
integer, dimension(2)                      :: jroot
integer, dimension(31)                     :: istats
integer, dimension(22)                     :: rstats
! type(VODE_OPTS)                            :: options
!
! --- OpenMP variables --- 
integer :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
!
!--------------------------------------------------------------------------------------------------
!
! --- init error ID/message ---
ierr0=0
outstr0=''
!
! --- Abbreviations ---
m1=model%m1
n1=model%n1
m2=model%m2
n2=model%n2
m3=model%m3
n3=model%n3
!
nspec=rn%nspecies
!
!--------------------------------------------------------------------------------------------------
# ifdef rhd_krome_debug01
do i1=m1,n1
  do i2=m2,n2
     do i3=m3,n3
        do i=1,nspec
           if (isnan(model%quc(i1,i2,i3,i)) .eq. .true.) then
              write(*,*) 'NaN detected: ',i1,i2,i3,i
              model%quc(i1,i2,i3,i)=0.0
           end if
        end do
    end do ! i3
  end do ! i2
end do ! i1
# endif
!--------------------------------------------------------------------------------------------------
! --- assign density arrays if necessary ---
!
if (.not.assdens_flag) then 
  call chem_rn_assign_density()
  assdens_flag=.TRUE.
end if
!
!--------------------------------------------------------------------------------------------------
write(*,'(A)') 'KROME> sourcestep'
! -------------------------------------------------------------------------------
! --- set the modes for the DVode solver (BDF method) ---
!
lrw=22+9*nspec+2*nspec**2
liw=30+nspec
!
rtol   = 1.d-3
atol   = 1.d-30
itol   = 1       ! --- error handling variables rtol,atol only scalars ---
itask  = 1
istate = 1
iopt   = 0
mf     = 22
!
ncall=0
!
allocate(iwork(liw))
allocate(rwork(lrw))
!
! -------------------------------------------------------------------------------
! --- Solve EOS: compute temperature  ---
# ifdef rhd_krome_t01
call timing_start('KROME: eosinter')
# endif
!
allocate(T(m1:n1,m2:n2,m3:n3))
! --- Call OpenMP version of EOS routines (apparently in non-OpenMP region) ---
!$OMP PARALLEL DEFAULT(SHARED)
call eosinter_TOMP(model%rho, model%ei, T, m1,n1, m2,n2, m3,n3, 0,0,0)
!$OMP END PARALLEL
!
# ifdef rhd_krome_t01
call timing_stop( 'KROME: eosinter')
# endif
!
! -------------------------------------------------------------------------------
!# ifdef rhd_krome_debug01
!write(*,*) 'RN%NREAC...... = ',rn%nreac
!write(*,*) 'RN%NSPECIES... = ',rn%nspecies
!write(*,*) 'RN%NELEMENTS.. = ',rn%nelements
!write(*,*) 'RN%SPECIESNAME = ',rn%speciesname
!# endif
! -------------------------------------------------------------------------------
!
! --- metal abundance ---
allocate(metal(m1:n1,m2:n2,m3:n3))
if (metal_flag) then 
  !write(*,*) 'KROME> Taking metal from provided array.'
  metal(:,:,:)=(model%quc(:,:,:, ircn_iqucmetal))
else 
  !write(*,*) 'KROME> Taking metal abundance from parameter file.'
  metal(:,:,:)=par%chem_abumetal*model%rho(:,:,:)/1.67353E-24
end if    
!
! --- time main loop, start ---
!# ifdef rhd_krome_t01
call timing_start('KROME: source step main loop')
!# endif
!
! --- loops over all spatial indices ---
!$OMP PARALLEL DEFAULT(PRIVATE), &
!$OMP          SHARED( m1, n1, m2, n2, m3, n3, nspec, ircn_iquc, & 
!$OMP                 model, dtime, metal, T, lrw, liw, rn) 
!$OMP DO SCHEDULE(STATIC)
do i1=m1,n1
  do i2=m2,n2
     do i3=m3,n3
      ! --- Reset error flag ---
      ier_flag=0
      ! 
      ! --- Get number densities of previous time step ---
      allocate(y(nspec))
      do i=1,nspec
         y(i)=dble(model%quc(i1,i2,i3,ircn_iquc(i)))
      end do
      !
      ! --- Calculate the chemical reaction rates ---
      allocate(rat(rn%nreac))
      !# ifdef rhd_krome_t01 
      !call timing_start('KROME: getrates')
      !# endif                                                                                                                     
      ! KROME has the rates "baked" within it. Do we still need this?
      !call chem_rn_getrates(metal(i1,i2,i3), T(i1,i2,i3), outstr, ierr)
      !# ifdef rhd_krome_t01                                                                                                        
      !call timing_stop('KROME: getrates')
      call timing_start('KROME: solver DLSODES')
      !# endif
      ! 
      ! --- Reset time as initial time will be overwritten after call to DVODE ---
      dbltime0 = 0.0D0
      dbltime1 = dble(dtime)
      !
      ! --- Reset error flag, force solver to continue ---
      ier=1
      istate=1
      !
      ! --- Set accuracy and DVODE control (again) ---
      rtol   = 1.d-4
      atol   = 1.d-30
      itol   = 1       ! --- error handling variables rtol,atol only scalars ---
      itask  = 1
      iopt   = 0
      mf     = 22
      !
      ! KROME doesn't need these options, it has them already when it calls DLSODES
      ! --- Call KROME ---
      ! Have to call with this dummy Tgas for a cell, but how will this work
      ! with OpenMP?
      ! Eventually need to call krome with T(i1,i2,i3) but KROME is compiled
      ! in double precision while T(...) is single
      ! (KROME expects Tgas :: real*8)
      Tgas = T(i1,i2,i3)
      call krome(y(1:nspec), Tgas, dbltime1)
      !# ifdef rhd_krome_t01
      !call timing_stop('KROME: solver dvode')
      call timing_stop('KROME: solver DLSODES')
      !# endif
      !-----------------------------------------------------------------------
      ier=istate
      !
      if (ier .ne. 2) then 
        ier_flag=100
      end if
      !
      ! --- Check for negative abundances ---
      do i=1,nspec
        if (y(i).lt.0.0D0) then 
           write(*,*) & 
                'KROME> SOURCESTEP: negative particle density at ',i1,i2,i3,', species ', & 
                rn%speciesname(ircn_iquc(i)),y(i),' (Forced to zero.)'
           y(i)=0.0D0
        end if
      end do
      !
      ! --- Store resulting number densities ---
      if (ier_flag.lt.100) then 
        ! --- copy results into box ---
        do i=1,nspec
           model%quc(i1,i2,i3,ircn_iquc(i))=real(y(i))
        end do
      else 
        model%quc(i1,i2,i3,1:nspec)=0.0
        call chem_rn_error(ier_flag, ier, i1, i2, i3)
      end if
      deallocate(y)
      deallocate(rat)
      !          
    end do ! i3
  end do ! i2
end do ! i1
!$OMP END DO
!$OMP END PARALLEL

write(*,*) 'KROME LOOP DONE'
!
! --- free memory ---
deallocate(metal, T)
deallocate(iwork, rwork)
!
! --- time main loop, stop ---
# ifdef rhd_krome_t01
call timing_stop('KROME: source step main loop')
# endif
!
!
999 if (present(outstr)) outstr=outstr0
    if (present(ierr)) ierr=ierr0
!
end subroutine krome_rn_sourcestep


!----------***********------------------------------------------------------------------------------
subroutine chem_rn_jac(neq, T, y, ml, mu, pd, nrowpd)
!---------------------------------------------------------------------------------------------------
! NAME:
!   chem_rn_jac
!
! PURPOSE:
!   (Create Jacobian.) 
!   Here (case MF=22) this subroutine only provides a dummy argument. 
!   
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! INPUT:
!   neq   : number of species/equations (here = nspecies)
!   ml, mu: (word) conditional input for the variable IWORK (in routine vode) 
!           IWORK(1) = ML     These are the lower and upper
!           IWORK(2) = MU     half-bandwidths, respectively, of the
!                             banded Jacobian, excluding the main diagonal.
!                             The band is defined by the matrix locations
!                             (i,j) with i-ML .le. j .le. i+MU.  ML and MU
!                             must satisfy  0 .le.  ML,MU  .le. NEQ-1.
!                             These are required if MITER is 4 or 5, and
!                             ignored otherwise.  ML and MU may in fact be
!                             the band parameters for a matrix to which
!                             df/dy is only approximately equal.
!
! OUTPUT:
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
! 
! COMMENT: 
!   Generally, this routine supplies the Jacobian df/dy by loading PD as follows:
!   - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
!     the partial derivative of f(i) with respect to y(j). 
!     (Ignore the ML and MU arguments in this case.)
!   - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
!     df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
!     PD from the top down.
!   - In either case, only nonzero elements need be loaded.
!   - HERE: MF=22
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   25. 1.05 Written by S.Wedemeyer-Boehm, KIS, Freiburg
!   13.01.21 (M.Szydlarski) Updated argument list.  
!--------------------------------------------------------------------------------------------------
!
implicit none
!
integer, intent(in)                                  :: neq, nrowpd, ml, mu
real(kind=drk), intent(in)                           :: T 
real(kind=drk), intent(in), dimension(neq)           :: y
real(kind=drk), intent(inout), dimension(nrowpd,neq) :: pd
!
return
end subroutine chem_rn_jac


!----------************----------------------------------------------------------------------------
subroutine chem_rn_func(neq, T, y, ydot)
!--------------------------------------------------------------------------------------------------
! NAME:
!   chem_rn_func
!
! PURPOSE:
!   calculate temporal derivative vector YDOT by adding rates
!   
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! INPUT:
!   neq        : number of species/equations (here = nspecies)
!   T          : gas temperature 
!   y          : number densities (y has values at index 1:neq)
!   rpar, ipar : user-defined scalars(/arrays) for solver input, here unused
!
! OUTPUT:
!   ydot       : temporal derivative vector
!
! ROUTINES:
!
! MODULES:
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
!   25. 1.05 Written by S.Wedemeyer-Boehm, KIS, Freiburg
!   13.01.21 (M.Szydlarski) Updated argument list.  
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O ---
integer,intent(in)                            :: neq 
real(kind=drk),intent(in)                     :: T 
real(kind=drk),intent(in),   dimension(neq)   :: y
real(kind=drk),              dimension(neq+1) :: ycp
real(kind=drk),intent(out),dimension(neq)   :: ydot
!
! --- local ---
integer                                       :: i, j, jj
real(kind=drk), dimension(rn%nreac)           :: r
!real(kind=drk), dimension()           :: r
!
!--------------------------------------------------------------------------------------------------
! --- count calls of this routine ---
ncall = ncall + 1
!
! --- copy input densities to a new vector with an additional element, 
!     set the additional element to zero ---
ycp(1:neq)=y(1:neq)
ycp(neq+1)=1.D0
!
! --- calculate rate vector ---
r=rat
do i=1,nreactants 
  r=r*ycp(rn%ispidrat(1:rn%nreac,i))
end do
!
! --- initialize derivatives ---
ydot(1:neq)=0.D0
!
! --- calculate derivative vector ---
do i=1,neq
  ydot(i)=sum(ispecact(i,1:rn%nreac)*r)
end do
!
return 
end subroutine chem_rn_func


!----------************************----------------------------------------------------------------
subroutine chem_rn_get_speciesindex(creqspec, iisp)
!--------------------------------------------------------------------------------------------------
! NAME:
!   chem_rn_get_speciesindex 
!
! PURPOSE:
!   get index of QUC array that contains the requested species
!
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! INPUT:
!   creqspec: name of requested chemical species
!
! OUTPUT:
!   iisp    : index of requested species in QUC array
!
! ROUTINES:
!
! MODULES:
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
!   25. 1.05 Written by S.Wedemeyer-Boehm, KIS, Freiburg
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O ---
integer                                 :: iisp
character(len=7)                        :: creqspec
!
! --- Local variables ---
integer                                 :: i
integer                                 :: nlen_rs, nlen_tmp
character(len=7)                        :: ctmp
!
!--------------------------------------------------------------------------------------------------
!
creqspec=adjustl(creqspec)
nlen_rs=len_trim(creqspec)
!
! --- assign number density arrays if necessary ---
!
if (.not.assdens_flag) then 
  call chem_rn_assign_density()
  assdens_flag=.TRUE.
end if
!
iisp=-1
!
do i=1,nquc_sp
  ctmp=adjustl(cquc_spname(i))
  nlen_tmp=len_trim(ctmp)
  !
  if ((ctmp(1:nlen_rs) == creqspec(1:nlen_rs)).and.(nlen_rs == nlen_tmp)) then 
    ! --- found requested species ---
    iisp=i
  endif
end do
!
end subroutine  chem_rn_get_speciesindex


!----------**************************--------------------------------------------------------------
subroutine chem_rn_get_speciesindices(creqspec, nisp, iisp, nmult)
!--------------------------------------------------------------------------------------------------
! NAME:
!   
!
! PURPOSE:
!   get indices (QUC arrays) of all chemical species that include the requested species
!   
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
! 
! INPUT:
!   creqspec: name of requested chemical species
!
! OUTPUT:
!   nisp    : number of species containing the requested species
!   iisp    : indices of QUC arrays for all matching species  
!   nmult   : number of occurence of requested species in each matching species
!             (e.g. nmult=[1,2,1] for requested element H and species set [H,H2,CH])
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!   ATTENTION! This is preliminary version was only tested for chemical species containing 
!              carbon (C)!!!
!              Important: Dimension nquc_sp must be known!
!
! PROCEDURE:
!
! EXAMPLE:
!   Consider a chemical reaction network [H,C,O,CO,CH,OH]
!   call chem_rn_get_speciesindices('C', isp_n, isp_c, nsp_cmult)
!   This call now returns the indices isp_c (a vector with isp_n elements) 
!   of the species C, CO, and CH. nsp_cmult is always 1 as the requested
!   element C is only contained once in all matching species.
!
! MODIFICATION HISTORY:
!   25. 1.05 Written by S.Wedemeyer-Boehm, KIS, Freiburg
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O ---
character(len=7)                        :: creqspec
integer                                 :: nisp
integer, dimension(:), intent(inout)    :: iisp, nmult
!
! --- Local variables ---
integer                                 :: i,j
!
!--------------------------------------------------------------------------------------------------
!
creqspec=adjustl(creqspec)
!
nisp=0
iisp(:)=-1
nmult(:)=0
!
! --- assign number density arrays if necessary ---
!
if (.not.assdens_flag) then 
  call chem_rn_assign_density()
  assdens_flag=.TRUE.
end if
!
! --- find index of wanted chemical element ---
!
do i=1,rn%nelements
  if (creqspec(1:2) == rn%elementname(i)) then 
    do j=1,rn%nspecies
      if (rn%natomconst(j,i).gt.0) then
        nisp=nisp+1
        iisp(nisp)=ircn_iquc(j)
        nmult(nisp)=rn%natomconst(j,i)
      end if
    end do
  end if   
end do
!
return
end subroutine  chem_rn_get_speciesindices


!----------***********************-----------------------------------------------------------------
subroutine chem_rn_init_co_indices()
!--------------------------------------------------------------------------------------------------
! NAME:
!   
!
! PURPOSE:
!   prepare indices of carbon containing chemical species C and CO (QUC arrays).
!   
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! INPUT:
!
! OUTPUT:
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE:
!
! EXAMPLE:
!   In rhd_rad_module: 
!   call chem_rn_init_co_indices()
!   Gets indices of all C containing species in order to calculate 
!   n(CO)/n(C) that is needed for deriving the CO opacity for 
!   radiative transfer with non-equilibrium CO opacity band. 
!
! MODIFICATION HISTORY:
!   16. 2.05 Written by S.Wedemeyer-Boehm, KIS, Freiburg
!--------------------------------------------------------------------------------------------------
!
implicit none
!
integer                                 :: iquc
character(len=7)                        :: ctmp
!
!--------------------------------------------------------------------------------------------------
!
write(6,'(A)') ' Using provided chemical species for CO radiative cooling.'
!
! --- get index of CO in QUC array ---
isp_co=-1
ctmp(:)=' '
ctmp(1:2)='CO'
call chem_rn_get_speciesindex(ctmp, isp_co)
if (isp_co.lt.0) then 
  write(6,'(A)') 'KROME> ERROR! Required chemical species CO not found!'
  stop
end if
!
! --- get index of all carbon containing chemical species in QUC array ---
isp_c(:)=-1
isp_n=0
ctmp(:)=' '
ctmp(1:1)='C'
call chem_rn_get_speciesindices(ctmp, isp_n, isp_c, nsp_cmult)
if (isp_n.le.0)  then 
   write(6,'(A)') 'KROME> ERROR! Required chemical species not found!'
   stop
else 
  write(6,'(A,I3,A)') '> Found ',isp_n,' carbon containing chemical species:'
  do iquc=1,isp_n 
    write(6,'(5X,A7," (N_C = ",I1,")")') cquc_spname(isp_c(iquc)), nsp_cmult(iquc)
  end do
end if
!
end subroutine chem_rn_init_co_indices


!----------*************---------------------------------------------------------------------------
subroutine chem_rn_error(ier_flag, ier, i1, i2, i3)
!--------------------------------------------------------------------------------------------------
! NAME:
!   
!
! PURPOSE:
!   Handle errors.
!   
! CATEGORY:
!   Chemistry
!
! CALLING SEQUENCE:
!
! INPUT:
!
! OUTPUT:
!
! ROUTINES:
!
! MODULES:
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
!   17. 8.06 Written by S.Wedemeyer-Boehm, KIS, Freiburg
!--------------------------------------------------------------------------------------------------
use rhd_gl_module
!
implicit none
!
integer, intent(in) :: ier_flag, ier, i1, i2, i3
!
!--------------------------------------------------------------------------------------------------
!
write(gl%nc_p,'(A,I3,A,"[",3(I3,",",1X),"]")') 'KROME> SOURCESTEP: ERROR! === IER = ',ier,' at ',i1,i2,i3
!
select case(ier) 
  case (1)
    write(gl%nc_p,'(A)') '         Nothing was done, as initial and final time were equal or dtime=0.'
  case (-1) 
    write(gl%nc_p,'(A)') '         Excess work done on this call. (Perhaps wrong MF.)'
  case (-2) 
    write(gl%nc_p,'(A)') '         Excess accuracy requested. (Tolerances too small.)'
  case (-3) 
    write(gl%nc_p,'(A)') '         Illegal input detected.'
  case (-4) 
    write(gl%nc_p,'(A)') '         Repeated error test failures. (Check all input.)'
  case (-5) 
    write(gl%nc_p,'(A)') '         Repeated convergence failures. (Perhaps bad Jacobian supplied or '
    write(gl%nc_p,'(A)') '         wrong choice of MF or tolerances.)'
  case (-6) 
    write(gl%nc_p,'(A)') '         Error weight became zero during problem. (Solution component i '
    write(gl%nc_p,'(A)') '         vanished, and ATOL or ATOL(i) = 0.)'
  case default 
    write(gl%nc_p,'(A)') '         Unknown error.'
end select
!
stop 
!
end subroutine chem_rn_error

end module krome_rn_module
