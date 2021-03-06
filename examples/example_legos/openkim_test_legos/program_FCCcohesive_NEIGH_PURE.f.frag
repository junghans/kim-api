!*******************************************************************************
!**
!**  PROGRAM TEST_NAME_STR
!**
!**  KIM compliant program to find (using the Golden section search algorithm)
!**  the minimum energy of one particle in a periodic FCC crystal (spec="SPECIES_NAME_STR") as a
!**  function of lattice spacing.
!**
!**  Works with the following NBC methods:
!**        NEIGH_PURE_H
!**        NEIGH_PURE_F
!**
!*******************************************************************************

!-------------------------------------------------------------------------------
!
! Main program
!
!-------------------------------------------------------------------------------
program TEST_NAME_STR
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  ! parameters controlling behavior of test
  !
  character(len=KIM_KEY_STRING_LENGTH), parameter :: testname = "TEST_NAME_STR"
  character(len=KIM_KEY_STRING_LENGTH), parameter :: testkimfile = "descriptor.kim"
  character(len=2), parameter :: specname    = 'SPECIES_NAME_STR'
  real(c_double),   parameter :: TOL         = 1.0e-8_cd
  real(c_double),   parameter :: FCCspacing  = FCC_SPACING_STR
  real(c_double),   parameter :: MinSpacing  = 0.800_cd*FCCspacing
  real(c_double),   parameter :: MaxSpacing  = 1.200_cd*FCCspacing
  integer(c_int),   parameter :: DIM         = 3
  integer(c_int),   parameter :: SupportHalf = 1  ! True

  ! significant local variables
  !
  real(c_double)                      :: rcut  ! cutoff radius of the potential
  type(neighObject_type), target      :: neighObject   ! neighbor list object
  real(c_double)                      :: FinalSpacing  ! lattice parameter
  real(c_double)                      :: FinalEnergy   ! energy/particle of crystal
                                                       ! at current spacing
  integer(c_int)                      :: CellsPerRcut  ! number of unit cells
                                                       ! along box (of size 
                                                       ! rcut) side

  ! KIM variables
  !
  character(len=KIM_KEY_STRING_LENGTH) :: modelname  ! KIM-compliant model name
  type(c_ptr)    :: pkim          ! pointer to KIM API object
  integer(c_int) :: N             ! number of particles

  ! other variables
  !
  real(c_double), external  ::  get_model_cutoff_firsttime

!========================= END VARIABLE DEFINITIONS ==========================


  ! Read in KIM Model name to use
  !
  print '("Please enter a valid KIM model name: ")'
  read(*,*) modelname


  ! Get model cutoff radius and compute number of particles needed
  ! (We need 2*cutoff, use 2.125*cutoff for saftey)
  !
  rcut = get_model_cutoff_firsttime(testkimfile, modelname)
  CellsPerRcut = ceiling(rcut/MinSpacing)
  N = 4*((2*CellsPerRcut)**3)


  ! Setup the KIM API object
  !
  call setup_KIM_API_object(pkim, testkimfile, modelname, N, specname, SupportHalf)


  ! allocate storage for neighbor lists, compute them for the first time,
  ! and store necessary pointers in KIM API object
  !
  allocate(neighObject%neighborList(N+1, N))
  call setup_neighborlist_KIM_access(pkim, neighObject)


  ! find equilibrium spacing by minimizing coheseive energy with respect
  ! to the periodic box size
  !
  call NEIGH_PURE_compute_equilibrium_spacing(pkim, &
         DIM,CellsPerRcut,MinSpacing,MaxSpacing, &
         TOL,N,neighObject,.false.,              &
         FinalSpacing,FinalEnergy)

  ! print results to screen
  !
  print '(80(''-''))'
  print '("This is Test          : ",A)', trim(testname)
  print '("Results for KIM Model : ",A)', trim(modelname)
  print *
  print '("Found minimum energy configuration to within",ES25.15)', TOL
  print *
  print '("Energy/part = ",ES25.15,"; Spacing = ",ES25.15)', FinalEnergy, &
        FinalSpacing
  print '(80(''-''))'


  ! Don't forget to free and/or deallocate
  !
  deallocate(neighObject%neighborList)
  call free_KIM_API_object(pkim)

  stop

end program TEST_NAME_STR

!-------------------------------------------------------------------------------
!
! NEIGH_PURE_compute_equilibrium_spacing :
!
!    Use the Golden section search algorithm to find the equilibrium spacing by
!    minimizing the energy of the system with respect to the periodic box size.
!
!-------------------------------------------------------------------------------
subroutine NEIGH_PURE_compute_equilibrium_spacing(pkim, &
             DIM,CellsPerRcut,MinSpacing,MaxSpacing, &
             TOL,N,neighObject,verbose,              &
             RetSpacing,RetEnergy)
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use mod_neighborlist
  implicit none
  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !-- Transferred variables
  type(c_ptr),    intent(in)            :: pkim
  integer(c_int), intent(in)            :: DIM
  integer(c_int), intent(in)            :: CellsPerRcut
  real(c_double), intent(in)            :: MinSpacing
  real(c_double), intent(in)            :: MaxSpacing
  real(c_double), intent(in)            :: TOL
  integer(c_int), intent(in)            :: N
  type(neighObject_type), intent(inout) :: neighObject
  logical,        intent(in)            :: verbose
  real(c_double), intent(out)           :: RetSpacing
  real(c_double), intent(out)           :: RetEnergy

  !-- Local variables
  real(c_double), parameter :: Golden = (1.0_cd + sqrt(5.0_cd))/2.0_cd
  integer(c_int)            :: ier, idum
  real(c_double)            :: Spacings(4)
  real(c_double)            :: Energies(4)
  integer(c_int)            :: MiddlePartId
  real(c_double), pointer   :: energy;         type(c_ptr) :: penergy
  real(c_double), pointer   :: coords(:,:);    type(c_ptr) :: pcoor
  real(c_double), pointer   :: cutoff;         type(c_ptr) :: pcutoff
  real(c_double), parameter :: cutpad = CUTOFF_PADDING_STR ! cutoff radius padding
  logical :: halfflag  ! .true. = half neigh list; .false. = full neigh list
  character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method

  ! Unpack data from KIM object
  !
  call kim_api_getm_data(pkim, ier, &
       "energy",      penergy,  1, &
       "coordinates", pcoor,    1, &
       "cutoff",      pcutoff,  1)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_getm_data", ier)
     stop
  endif
  call c_f_pointer(penergy, energy)
  call c_f_pointer(pcoor, coords, [DIM,N])
  call c_f_pointer(pcutoff, cutoff)

  ! determine which neighbor list type to use
  !
  ier = kim_api_get_nbc_method(pkim, NBC_Method)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_get_nbc_method", ier)
     stop
  endif
  if (index(NBC_Method,"NEIGH_PURE_H").eq.1) then
     halfflag = .true.
  elseif (index(NBC_Method,"NEIGH_PURE_F").eq.1) then
     halfflag = .false.
  else
     ier = KIM_STATUS_FAIL
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "Unknown NBC method", ier)
     return
  endif

  ! Initialize for minimization
  !
  Spacings(1) = MinSpacing
  call create_FCC_configuration(Spacings(1), 2*CellsPerRcut, .true., coords, &
                                MiddlePartId)
  ! compute new neighbor lists (could be done more intelligently, I'm sure)
  call NEIGH_PURE_periodic_neighborlist(halfflag, N, coords, (cutoff+cutpad), &
                                        MiddlePartId, neighObject)
  ier = kim_api_model_compute(pkim)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_model_compute", ier)
     stop
  endif
  Energies(1) = energy
  if (verbose) &
     print '("Energy/part = ",ES25.15,"; Spacing = ",ES25.15)', Energies(1), &
           Spacings(1)

  ! setup and compute for max spacing
  Spacings(3) = MaxSpacing
  call create_FCC_configuration(Spacings(3), 2*CellsPerRcut, .true., coords, &
                                MiddlePartId)
  ! compute new neighbor lists (could be done more intelligently, I'm sure)
  call NEIGH_PURE_periodic_neighborlist(halfflag, N, coords, (cutoff+cutpad), &
                                        MiddlePartId, neighObject)
  ! Call model compute
  ier = kim_api_model_compute(pkim)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_model_compute", ier)
     stop
  endif
  Energies(3) = energy
  if (verbose) &
     print '("Energy/part = ",ES25.15,"; Spacing = ",ES25.15)', Energies(3), &
           Spacings(3)

  ! setup and compute for first intermediate spacing
  Spacings(2) = MinSpacing + (2.0_cd - Golden)*(MaxSpacing - MinSpacing)
  call create_FCC_configuration(Spacings(2), 2*CellsPerRcut, .true., coords, &
                                MiddlePartId)
  ! compute new neighbor lists (could be done more intelligently, I'm sure)
  call NEIGH_PURE_periodic_neighborlist(halfflag, N, coords, (cutoff+cutpad), &
                                        MiddlePartId, neighObject)
  ! Call model compute
  ier = kim_api_model_compute(pkim)
  if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_model_compute", ier)
     stop
  endif
  Energies(2) = energy
  if (verbose) &
     print '("Energy/part = ",ES25.15,"; Spacing = ",ES25.15)', Energies(2), &
           Spacings(2)


  ! iterate until convergence.
  !
  do while (abs(Spacings(3) - Spacings(1)) .gt. TOL)
     ! set new spacing
     Spacings(4) = (Spacings(1) + Spacings(3)) - Spacings(2)
     ! compute new particle coordinates based on new spacing
     call create_FCC_configuration(Spacings(4), 2*CellsPerRcut, .true., coords,&
                                   MiddlePartId)
     ! compute new neighbor lists (could be done more intelligently, I'm sure)
     call NEIGH_PURE_periodic_neighborlist(halfflag,N,coords,(cutoff+cutpad), &
                                           MiddlePartId, neighObject)
     ! Call model compute
     ier = kim_api_model_compute(pkim)
     if (ier.lt.KIM_STATUS_OK) then
        idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                    "kim_api_model_compute", ier)
        stop
     endif
     Energies(4) = energy
     if (verbose) &
        print '("Energy/part = ",ES25.15,"; Spacing = ",ES25.15)', Energies(4),&
              Spacings(4)

     ! determine the new interval
     if (Energies(4) .lt. Energies(2)) then
        ! We want the right-hand interval
        Spacings(1) = Spacings(2); Energies(1) = Energies(2)
        Spacings(2) = Spacings(4); Energies(2) = Energies(4)
     else
        ! We want the left-hand interval
        Spacings(3) = Spacings(1); Energies(3) = Energies(1)
        Spacings(1) = Spacings(4); Energies(1) = Energies(4)
     endif
  enddo

  ! pull out results and return
  !
  RetSpacing = Spacings(2)
  RetEnergy  = Energies(2)

  return

end subroutine NEIGH_PURE_compute_equilibrium_spacing
