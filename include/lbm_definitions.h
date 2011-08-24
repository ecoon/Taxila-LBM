  ! da # of dofs enum
  PetscInt, parameter:: ONEDOF = 1
  PetscInt, parameter:: NPHASEDOF = 2
  PetscInt, parameter:: NPHASEXBDOF = 3
  PetscInt, parameter:: NFLOWDOF = 4
  PetscInt, parameter:: NSPECIEDOF = 5
  PetscInt, parameter:: NSPECIEXBDOF = 6

  ! discretizations
  PetscInt, parameter:: NULL_DISCRETIZATION = 0
  PetscInt, parameter:: D3Q19_DISCRETIZATION = 1
  PetscInt, parameter:: D2Q9_DISCRETIZATION = 2

  ! walls type enum
  PetscInt, parameter:: WALLS_TYPE_PETSC = 1
  PetscInt, parameter:: WALLS_TYPE_FUNCTION = 2

  ! string lengths
  PetscInt, parameter :: MAXHEADERLENGTH = 2048
  PetscInt, parameter :: MAXSTRINGLENGTH = 512
  PetscInt, parameter :: MAXWORDLENGTH = 32
  PetscInt, parameter :: MAXIODIGITS = 3

  ! directions
  PetscInt, parameter :: X_DIRECTION = 1
  PetscInt, parameter :: Y_DIRECTION = 2
  PetscInt, parameter :: Z_DIRECTION = 3

  ! Flow BCs
  PetscInt, parameter :: BC_PERIODIC = 0
  PetscInt, parameter :: BC_PSEUDOPERIODIC = 1
  PetscInt, parameter :: BC_FLUX = 2
  PetscInt, parameter :: BC_DIRICHLET = 3
  PetscInt, parameter :: BC_ZERO_GRADIENT = 4
  PetscInt, parameter :: BC_VELOCITY = 5
  
  ! boundaries
  PetscInt, parameter :: BOUNDARY_XM = 1
  PetscInt, parameter :: BOUNDARY_XP = 2
  PetscInt, parameter :: BOUNDARY_YM = 3
  PetscInt, parameter :: BOUNDARY_YP = 4
  PetscInt, parameter :: BOUNDARY_ZM = 5
  PetscInt, parameter :: BOUNDARY_ZP = 6

  ! cardinal directions: 
  ! CARDINAL_NORMAL cross CARDINAL_CROSS = CARDINAL_RESULTANT
  ! under a right-hand-rule
  PetscInt, parameter :: CARDINAL_NORMAL = 1
  PetscInt, parameter :: CARDINAL_CROSS = 2
  PetscInt, parameter :: CARDINAL_RESULTANT = 3

  ! fluid types
  PetscInt, parameter :: NULL_FLUID_TYPE = 0
  PetscInt, parameter :: PHASE_FLUID_TYPE = 0
  PetscInt, parameter :: SPECIE_FLUID_TYPE = 0

  ! relaxation types
  PetscInt, parameter :: RELAXATION_MODE_SRT = 0
  PetscInt, parameter :: RELAXATION_MODE_MRT = 1

  ! max cruft for bags, these can be arbitrarily
  ! increased with no memory problems as nothing
  ! is actually allocated of this size
  PetscInt, parameter :: NMAX_PHASES = 5
  PetscInt, parameter :: NMAX_DIRECTIONS = 27
  

  ! cruft for use with pflotran
  PetscInt, parameter :: LBM_MODE = 8

  ! EOS types, see Yuan and Schaeffer '06
  PetscInt, parameter :: EOS_NULL = 0
  PetscInt, parameter :: EOS_DENSITY = 1    ! psi = rho
  PetscInt, parameter :: EOS_SC = 2         ! psi = rho0*(1-exp(-rho)), Shan & Chen '93/'94
  PetscInt, parameter :: EOS_PR = 3         ! psi = cubic law from Peng & Robinson '??, 
  PetscInt, parameter :: EOS_THERMO = 4     ! psi = psi0*exp(-rho0/rho), see Shan & Chen '94
