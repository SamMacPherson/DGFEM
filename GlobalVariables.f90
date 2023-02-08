Module GlobalVariables
implicit none

! Use double precision for the real numbers
integer, parameter :: dp = 8

! Gas parameters
real(dp)                                :: gammaGas 
real(dp)                                :: GasConstant, mu0,T0, S, constMu, thermalConductivity
logical                                 :: useSutherland

! Time
real(dp) :: totalTime, dt, residualStop, CFL
integer  :: outputFrequency, RKOrder, fluxScheme, BCScheme,testProblem

! Order of polynomial, cubature, and Gauss
real(dp) :: N
integer  :: Np, Nfp,Ncub, NGauss

! Limiter
integer :: useLimiter
! Nodal locations
real(dp), dimension(:,:), allocatable   :: x, y

! Grid details
integer, dimension(:,:), allocatable    ::  EToE, EToV, BCType
real(dp), dimension(:,:), allocatable   ::  VXMatrix, VYMatrix
real(dp), dimension(:), allocatable     ::  VX, VY
real(dp), dimension(:,:), allocatable   ::  NodexMatrix, NodeyMatrix
real(dp), dimension(:), allocatable     ::  Nodex, Nodey  

! Number of elements
integer:: K

! Matrices to read in cubature data
real(dp), dimension(:,:), allocatable   :: cubDr, cubDs, cubV, cubW, cubrx, cubry, cubsx, cubsy 

! Matrices to read in Gauss data
real(dp), dimension(:,:), allocatable   :: gaussinterp, gaussW, gaussnx, gaussny
integer, dimension(:,:), allocatable    :: gaussmapMMatrix, gaussmapPMatrix 
integer, dimension(:,:), allocatable    :: gaussmapBMatrix, gaussmapIMatrix, gaussmapOMatrix ,gaussmapWMatrix ,gaussmapMWMatrix, gaussmapSMatrix, gaussmapNMatrix, gaussmapDMatrix

! Matrices to read in Misc.
real(dp), dimension(:,:), allocatable   :: V, J, Fscale, massMatrix, divJ
integer, dimension(:,:), allocatable    :: vmapMMatrix, vmapPMatrix, vmapIMatrix, vmapOMatrix ,vmapWMatrix ,vmapMWMatrix, vmapSMatrix, vmapNMatrix, vmapDMatrix
! Define a vector for maps
integer, dimension(:), allocatable      :: vmapM, vmapP, vmapI, vmapO, vmapW, vmapMW, vmapS, vmapN, vmapD
integer, dimension(:), allocatable      :: gaussmapI, gaussmapO, gaussmapW, gaussmapMW, gaussmapS, gaussmapN, gaussmapD, gaussmapB, gaussmapM, gaussmapP

! Boundary condition markers
integer, parameter  ::  Inflow = 1, Outflow = 2, Wall = 3, Slip = 4, Dirichlet = 5, Neumann = 6, Moving = 7

! LSERK Coefficients
real(dp), dimension(5)  :: aL, bL, cL

end module GlobalVariables