Module BCandIC
! Module for boundary initial condition functions

use GlobalVariables
use PrimConTransformation
contains

subroutine ShearFlow2DIC(Q)
    ! Computes initial conditions for Hesthaven and Warbution CNS test case I
    ! Arguments:    Q   -   Array of conserved variables

    ! Author Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(out)             ::  Q
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  ::  rho, rhou, rhov, ener
    
    real(dp)    ::  mu

    mu   = 0.01

    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    rho = 1
    rhou = y**2
    rhov = 0
    ener = (2*mu*x + 10)/(gammaGas - 1) + 0.5*y**4

    Q(:,:,1) = rho
    Q(:,:,2) = rhou
    Q(:,:,3) = rhov
    Q(:,:,4) = ener
end subroutine ShearFlow2DIC

subroutine ShearFlow2DBC(Q,mapB, xin, yin)
    ! Computes boundary conditions for Hesthaven and Warbution CNS test case I
    ! Arguments:    Q       -   Array of conserved variables
    !               mapB    -   Map of Boundary nodes
    !               xin     -   Corresponding array of x-values
    !               yin     -   Corresponding array of y-values

    ! Author Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(out)             ::  Q
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  ::  rho, rhou, rhov, ener, xin, yin
    integer, dimension(:), intent(in)                   ::  mapB

    real(dp)    :: mu
    ! Flat arrays
    real(dp), dimension(:), allocatable ::  rhoflat, rhouflat, rhovflat, enerflat, xflat, yflat

    ! Mask for unpacking variables
    logical, dimension(size(Q(:,1,1)),size(Q(1,:,1))) :: unPackMask 
    unPackMask = .true.

    mu = 0.01
    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    ! Flatten arrays to be accessed by map
    rhoflat  = pack(rho,.true.)
    rhouflat = pack(rhou,.true.)
    rhovflat = pack(rhov,.true.)
    enerflat = pack(ener,.true.)
    xflat = pack(xin,.true.)
    yflat = pack(yin,.true.)

    ! Apply BC's
    rhoflat(mapB)  = 1
    rhouflat(mapB) = yflat(mapB)**2
    rhovflat(mapB) = 0
    enerflat(mapB) = (2*mu*xflat(mapB) + 10)/(gammaGas - 1) + 0.5*(yflat(mapB)**4)

    ! Reform matrices
    rho = unpack(rhoflat, unPackMask, rho)
    rhou = unpack(rhouflat, unPackMask, rhou)
    rhov = unpack(rhovflat, unPackMask, rhov)
    ener = unpack(enerflat, unPackMask, ener)

    ! Reform Q-vector
    Q(:,:,1) = rho
    Q(:,:,2) = rhou
    Q(:,:,3) = rhov
    Q(:,:,4) = ener
end subroutine ShearFlow2DBC

subroutine ZeitounTestCaseIC(Q)
    ! Computes initial condition for Zeitoun 2009 test case for Microshock tube
    ! Arguments:    Q   -   Array of conserved variables

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) :: Q

    ! Conserved variables
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  :: rho, rhou, rhov, ener

    ! Primitive variables
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  :: u, v, p

    ! Diaphragm location
    real(dp)    :: xd

    ! Flow conditions
    real(dp)    :: p4, p1, T

    integer :: i, j

    ! Extract conserved variables
    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    ! Convert to primitive variables
    call ConToPrim(Q,rho,u,v,p)


    xd = 29.6e-3_dp
    ! Initial conditions
   

    ! Flow Conditions
    p1 = 44.2_dp
    p4 = 11.9_dp*p1
    T = 300.0_dp
    ! Loop over elements and set initial conditions depending on whether the node is in the driver or driven section
    do i = 1, size(Q(:,1,1))
        do j = 1, size((Q(1,:,1)))

            if(x(i,j).le.xd) then
                rho(i,j) = p4/GasConstant/T
                p(i,j)   = p4
            else
                rho(i,j) = p1/GasConstant/T
                p(i,j)   = p1
            end if

        end do
    end do

    ! Initial velocity everywhere is zero
    u = 0.0_dp
    v = 0.0_dp

    ! Convert primitive back to conserved variables
    call PrimToCon(Q,rho,u,v,p)

end subroutine ZeitounTestCaseIC

subroutine ZeitounTestCaseNoSlipBC(Q, T, mapW)
    ! Computes no-slip boundary conditions for Zeitoun test case
    ! Arguments:    Q       -   Array of conserved variables
    !               mapW    -   Map of wall nodes   

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) ::  Q
    real(dp), dimension(:,:), intent(out)   :: T
    integer, dimension(:), intent(in)       ::  mapW

    ! Conserved variables
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  :: rho, rhou, rhov, ener

    ! Primitive variables
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  :: u, v, p

    ! unPackMask
    logical, dimension(size(Q(:,1,1)),size(Q(1,:,1)))   :: unPackMask

    ! Flat arrays to be accessed by maps
    real(dp), dimension(:), allocatable ::  uflat, vflat, rhoflat, pflat, xflat, Tflat

    real:: start, finish
    integer:: i

    unPackMask = .true.


    ! Extract conserved variables
    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    ! Convert to primitive variables
    call ConToPrim(Q,rho,u,v,p)

    rhoflat = pack(rho,.true.)
    uflat   = pack(u,.true.)
    vflat   = pack(v,.true.)
    pflat   = pack(p,.true.)
    Tflat   = pack(T,.true.)

    ! No slip conditions, and Tw = 300K
    uflat(mapW) = 0.0_dp
    vflat(mapW) = 0.0_dp
    Tflat(mapW) = 300.0_dp
    pflat(mapW) = rhoflat(mapW)*GasConstant*Tflat(mapW)
    
    ! Reform matrices
    rho = unpack(rhoflat,unPackMask,rho)
    u   = unpack(uflat,unPackMask,u)
    v   = unpack(vflat,unPackMask,v)
    p   = unpack(pflat,unPackMask,p)
    T   = unpack(Tflat,unPackMask,T)
    
    ! Convert back to conserved variables
    call PrimToCon(Q,rho,u,v,p)

    deallocate(uflat,vflat,rhoflat,pflat)
end subroutine ZeitounTestCaseNoSlipBC

subroutine ZeitounTestCaseNoSlipBCLimiter(Q, T, mapW)
    ! Computes no-slip boundary conditions for Zeitoun test case
    ! Arguments:    Q       -   Array of conserved variables
    !               mapW    -   Map of wall nodes   

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) ::  Q
    real(dp), dimension(:,:), intent(in)    ::  T
    integer, dimension(:), intent(in)       ::  mapW

    ! Primitive variables
    real(dp), dimension(size(T(:,1)), size(T(1,:)))         ::  mu, u, v, rho, p
 

    ! unPackMask
    logical, dimension(size(Q(:,1,1)),size(Q(1,:,1)))   :: unPackMask

    ! Flat arrays to be accessed by maps
    real(dp), dimension(:), allocatable ::  rhoflat, uflat, vflat, pflat, Tflat

    ! Diaphragm location
    real(dp)    :: xd

    ! Flow conditions
    real(dp)    :: Tw


    integer:: i

    unPackMask = .true.
    Tw      = 300.0_dp

    ! Extract variables
    call ConToPrim(Q,rho,u,v,p)

    if (useSutherland) then
        call Sutherland(T,mu)
    else
        mu = constMu
    end if




    rhoflat     = pack(rho,.true.)
    uflat       = pack(u,.true.)
    vflat       = pack(v,.true.)
    pflat       = pack(p,.true.)
    Tflat       = pack(T,.true.)
 
    ! Ghost element approach
    uflat(mapW) = - uflat(mapW)
    vflat(mapW) = - vflat(mapW)
    Tflat(mapW) = 2*Tw - Tflat(mapW)
    pflat(mapW) = rhoflat(mapW)*GasConstant*Tflat(mapW)

    
    ! Reform matrices
    rho = unpack(rhoflat,unPackMask,rho)
    u   = unpack(uflat,unPackMask,u)
    v   = unpack(vflat,unPackMask,v)
    p   = unpack(pflat,unPackMask,p)

    call PrimToCon(Q,rho,u,v,p)

    deallocate(uflat,vflat,rhoflat,pflat)
end subroutine ZeitounTestCaseNoSlipBCLimiter

subroutine ZeitounTestCaseSlipBCGaussNodes(bQ, bT, cu, cT, mapW)
    ! Computes no-slip boundary conditions for Zeitoun test case
    ! Arguments:    Q       -   Array of conserved variables
    !               mapW    -   Map of wall nodes   

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) ::  bQ
    real(dp), dimension(:,:),   intent(out) ::  bT
    real(dp), dimension(:,:), intent(in)    ::  cT, cU
    integer, dimension(:), intent(in)       ::  mapW

    ! Conserved variables
    real(dp), dimension(size(bQ(:,1,1)),size(bQ(1,:,1)))  :: rho, rhou, rhov, ener

    ! Primitive variables
    real(dp), dimension(size(bQ(:,1,1)),size(bQ(1,:,1)))  :: u, v, p, mu, gdudy, gdTdy

    ! Derivatives
    real(dp), dimension(Np,K)  :: dudx, dudy, dTdx, dTdy

    ! unPackMask
    logical, dimension(size(bQ(:,1,1)),size(bQ(1,:,1)))   :: unPackMask

    ! Flat arrays to be accessed by maps
    real(dp), dimension(:), allocatable ::  uflat, vflat, rhoflat, pflat, Tflat, dTdyflat, dudyflat, lambdaflat, muflat, Prflat

    ! Diaphragm location
    real(dp)    :: xd

    ! Flow conditions
    real(dp)    :: Tw, alphau, alphaT, cp, pi

    integer:: i

    unPackMask = .true.


    alphau  = 1.142_dp
    alphaT  = 0.5865_dp
    cp      = 520.3_dp
    pi      = 3.1415927_dp
    Tw      = 300.0_dp

    ! Extract conserved variables
    rho  = bQ(:,:,1)
    rhou = bQ(:,:,2)
    rhov = bQ(:,:,3)
    ener = bQ(:,:,4)

    ! Convert to primitive variables
    call ConToPrim(bQ,rho,u,v,p)

    if (useSutherland) then
        call Sutherland(bT,mu)
    else
        mu = constMu
    end if

    ! Cubature and gauss nodes for temperature and u-velocity
    call DGGradientBC(cu,u,u,dudx,dudy)
    call DGGradientBC(cT,bT,bT,dTdx,dTdy)

   
    gdudy = abs(matmul(gaussinterp,dudy))
    gdTdy = abs(matmul(gaussinterp,dTdy))

    
    
    rhoflat     = pack(rho,.true.)
    uflat       = pack(u,.true.)
    vflat       = pack(v,.true.)
    pflat       = pack(p,.true.)
    Tflat       = pack(bT,.true.)
    dTdyflat    = pack(gdTdy,.true.)
    dudyflat    = pack(gdudy,.true.)
    muflat      = pack(mu,.true.)
    
    ! Prandtl number, and variable lambda as seen in Zeitoun 2009
    Prflat = cp*muflat/thermalConductivity
    lambdaflat = sqrt(pi/2.0_dp)*muflat/sqrt(pflat*rhoflat)

    ! slip conditions
    uflat(mapW) = alphau*lambdaflat(mapW)*dudyflat(mapW)
    vflat(mapW) = 0.0_dp
    ! Temperature jump conditions, with ideal gas law
    Tflat(mapW) = Tw + alphaT*gammaGas/(gammaGas-1.0_dp)*lambdaflat(mapW)/Prflat(mapW)*dTdyflat(mapW)
    pflat(mapW) = rhoflat(mapW)*GasConstant*Tflat(mapW)
    ! Reform matrices
    rho = unpack(rhoflat,unPackMask,rho)
    u   = unpack(uflat,unPackMask,u)
    v   = unpack(vflat,unPackMask,v)
    p   = unpack(pflat,unPackMask,p)
    bT  = unpack(Tflat,unPackMask,bT)
    
    ! Convert back to conserved variables
    call PrimToCon(bQ,rho,u,v,p)

    deallocate(dTdyflat,dudyflat,lambdaflat,Prflat)
    deallocate(uflat,vflat,pflat,rhoflat)
end subroutine ZeitounTestCaseSlipBCGaussNodes

subroutine ZeitounTestCaseSlipBCLimiter(Q, T, cu, gu, cT, gT, grho, mapW)
    ! Computes no-slip boundary conditions for Zeitoun test case
    ! Arguments:    Q       -   Array of conserved variables
    !               mapW    -   Map of wall nodes   

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) ::  Q
    real(dp), dimension(:,:), intent(in)    ::  cu, T, gu, cT, gT, grho
    integer, dimension(:), intent(in)       ::  mapW

    ! Primitive variables
    real(dp), dimension(size(T(:,1)), size(T(1,:)))         ::  u, v, rho, p, mu

    ! Derivatives
    real(dp), dimension(Np,K)   :: dudx, dudy, dTdx, dTdy
    real(dp), dimension(1,K)    :: dudyC0, dTdyC0
    real(dp), dimension(3,K)    :: dudyC, dTdyC   

    real(dp), dimension(3*NGauss,K)     :: gdudy, gdTdy, gmu
    real(dp), dimension(NGauss,3*K)     :: facesgdudy, facesgdTdy, facesgrho, facesgT, facesgmu
    logical, dimension(NGauss, 3*K)     :: unPackMaskFaces

    ! unPackMask
    logical, dimension(size(Q(:,1,1)),size(Q(1,:,1)))   :: unPackMask

    ! Flat arrays to be accessed by maps
    real(dp), dimension(:), allocatable ::  rhoflat, uflat, vflat, pflat, Tflat, dudyflat, dTdyflat, lambdaflat, muflat, Prflat
    real(dp), dimension(:), allocatable ::  gdudyflat, gdTdyflat, gmuflat, grhoflat, gTflat

    ! Diaphragm location
    real(dp)    :: xd

    ! Flow conditions
    real(dp)    :: Tw, alphau, alphaT, cp, pi

    real(dp), dimension(Np)                     :: aveflat
    real(dp), dimension(1,Np)                   :: ave
    integer:: i

    unPackMask = .true.
    unPackMaskFaces = .true.

    aveflat = sum(massMatrix,dim=1)/2.0_dp
    ave(1,:) = aveflat

    alphau  = 1.142_dp
    alphaT  = 0.5865_dp
    cp      = 520.3_dp
    pi      = 3.141592653_dp
    Tw      = 300.0_dp

    ! Extract variables
    call ConToPrim(Q,rho,u,v,p)

    if (useSutherland) then
        call Sutherland(T,mu)
    else
        mu = constMu
    end if

    ! Calculate gradients at polynomial nodes
    call DGGradientBC(cu,gu,gu,dudx,dudy)
    call DGGradientBC(cT,gT,gT,dTdx,dTdy)

  
    ! ! Use Gauss nodes to find average value of variables on each face
    ! gdudy =abs(matmul(gaussinterp,dudy))
    ! gdTdy =abs(matmul(gaussinterp,dTdy))

    ! gdudyflat   = pack(gdudy,.true.)
    ! gdTdyflat   = pack(gdTdy,.true.)
    ! gmuflat     = pack(gmu, .true.)
    ! grhoflat    = pack(grho,.true.)
    ! gTflat      = pack(gT,.true.) 

    ! facesgdudy  = unpack(gdudyflat, unPackMaskFaces, facesgdudy)
    ! facesgdTdy  = unpack(gdTdyflat, unPackMaskFaces, facesgdTdy)
    ! facesgmu    = unpack(gmuflat, unPackMaskFaces, facesgmu)
    ! facesgrho   = unpack(grhoflat, unPackMaskFaces,facesgrho)
    ! facesgT     = unpack(gTflat, unPackMaskFaces,facesgT)

    ! ! Find average of Gauss nodes on each face
    ! dudyflatboundary = sum(facesgdudy, dim = 1)/real(NGauss)
    ! dTdyflatboundary = sum(facesgdTdy, dim = 1)/real(NGauss)
    ! muflatboundary   = sum(facesgmu, dim = 1)/real(NGauss)
    ! rhoflatboundary  = sum(facesgrho, dim = 1)/real(NGauss)
    ! Tflatboundary  = sum(facesgT, dim = 1)/real(NGauss)

    
    dudy = abs(dudy)
    dTdy = abs(dTdy)
    
    dudyC0 = matmul(ave,dudy)
    
    dudyC(1,:) = dudyC0(1,EToE(:,1))
    dudyC(2,:) = dudyC0(1,EToE(:,2))
    dudyC(3,:) = dudyC0(1,EToE(:,3))

    dTdyC0 = matmul(ave,dTdy)
    
    dTdyC(1,:) = dTdyC0(1,EToE(:,1))
    dTdyC(2,:) = dTdyC0(1,EToE(:,2))
    dTdyC(3,:) = dTdyC0(1,EToE(:,3))

    ! Primitive variables at element centers
    rhoflat     = pack(rho,.true.)
    uflat       = pack(u,.true.)
    vflat       = pack(v,.true.)
    pflat       = pack(p,.true.)
    Tflat       = pack(T,.true.)
    dudyflat    = pack(dudyC,.true.)
    dTdyflat    = pack(dTdyC,.true.)
    muflat      = pack(mu,.true.)
    
    ! Prandtl number, and variable lambda(mean free path) as seen in Zeitoun 2009
    Prflat = cp*muflat/thermalConductivity
    lambdaflat = sqrt(pi/2.0_dp)*muflat/sqrt(pflat*rhoflat)

    ! slip conditions ghost element approach
    uflat(mapW)   = 2*alphau*lambdaflat(mapW)*dudyflat(mapW) - uflat(mapW)
    vflat(mapW)   = -vflat(mapW)
    Tflat(mapW)   = 2*(Tw + alphaT*gammaGas/(gammaGas-1.0_dp)*lambdaflat(mapW)/Prflat(mapW)*dTdyflat(mapW)) - Tflat(mapW)
    pflat(mapW)   = rhoflat(mapW)*GasConstant*Tflat(mapW)

    
    
    
    
    ! Reform matrices
    rho = unpack(rhoflat,unPackMask,rho)
    u   = unpack(uflat,unPackMask,u)
    v   = unpack(vflat,unPackMask,v)
    p   = unpack(pflat,unPackMask,p)

    call PrimToCon(Q,rho,u,v,p)

    
    deallocate(uflat,vflat,rhoflat,pflat)
end subroutine ZeitounTestCaseSlipBCLimiter

subroutine ZeitounTestCaseConstantSlipBC(Q, T, mapW)
    ! Computes no-slip boundary conditions for Zeitoun test case
    ! Arguments:    Q       -   Array of conserved variables
    !               mapW    -   Map of wall nodes   

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) ::  Q
    real(dp), dimension(:,:), intent(out)   :: T
    integer, dimension(:), intent(in)       ::  mapW

    ! Conserved variables
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  :: rho, rhou, rhov, ener

    ! Primitive variables
    real(dp), dimension(size(Q(:,1,1)),size(Q(1,:,1)))  :: u, v, p

    ! unPackMask
    logical, dimension(size(Q(:,1,1)),size(Q(1,:,1)))   :: unPackMask

    ! Flat arrays to be accessed by maps
    real(dp), dimension(:), allocatable ::  uflat, vflat, rhoflat, pflat, xflat, Tflat

    real:: start, finish
    integer:: i

    unPackMask = .true.


    ! Extract conserved variables
    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    ! Convert to primitive variables
    call ConToPrim(Q,rho,u,v,p)

    rhoflat = pack(rho,.true.)
    uflat   = pack(u,.true.)
    vflat   = pack(v,.true.)
    pflat   = pack(p,.true.)
    Tflat   = pack(T,.true.)

    ! No slip conditions, and Tw = 300K
    uflat(mapW) = 24.98_dp
    vflat(mapW) = 0.0_dp
    Tflat(mapW) = 300.0_dp
    pflat(mapW) = rhoflat(mapW)*GasConstant*Tflat(mapW)
    
    ! Reform matrices
    rho = unpack(rhoflat,unPackMask,rho)
    u   = unpack(uflat,unPackMask,u)
    v   = unpack(vflat,unPackMask,v)
    p   = unpack(pflat,unPackMask,p)
    T   = unpack(Tflat,unPackMask,T)
    
    ! Convert back to conserved variables
    call PrimToCon(Q,rho,u,v,p)

    deallocate(uflat,vflat,rhoflat,pflat)
end subroutine ZeitounTestCaseConstantSlipBC

subroutine ZeitounTestCaseConstantSlipBCLimiter(Q, T, mapW)
    ! Computes no-slip boundary conditions for Zeitoun test case
    ! Arguments:    Q       -   Array of conserved variables
    !               mapW    -   Map of wall nodes   

    ! Author: Sam MacPherson
    implicit none

    ! Define Arguments
    real(dp), dimension(:,:,:), intent(out) ::  Q
    real(dp), dimension(:,:), intent(in)    ::  T
    integer, dimension(:), intent(in)       ::  mapW

    ! Primitive variables
    real(dp), dimension(size(T(:,1)), size(T(1,:)))         ::  mu, u, v, rho, p
 

    ! unPackMask
    logical, dimension(size(Q(:,1,1)),size(Q(1,:,1)))   :: unPackMask

    ! Flat arrays to be accessed by maps
    real(dp), dimension(:), allocatable ::  rhoflat, uflat, vflat, pflat, Tflat

    ! Diaphragm location
    real(dp)    :: xd

    ! Flow conditions
    real(dp)    :: Tw


    integer:: i

    unPackMask = .true.
    Tw      = 300.0_dp

    ! Extract variables
    call ConToPrim(Q,rho,u,v,p)

    if (useSutherland) then
        call Sutherland(T,mu)
    else
        mu = constMu
    end if




    rhoflat     = pack(rho,.true.)
    uflat       = pack(u,.true.)
    vflat       = pack(v,.true.)
    pflat       = pack(p,.true.)
    Tflat       = pack(T,.true.)
 
    ! Ghost element approach
    uflat(mapW) = 2.0_dp*24.98_dp - uflat(mapW)
    vflat(mapW) = - vflat(mapW)
    Tflat(mapW) = 2.0_dp*Tw - Tflat(mapW)
    pflat(mapW) = rhoflat(mapW)*GasConstant*Tflat(mapW)

    
    ! Reform matrices
    rho = unpack(rhoflat,unPackMask,rho)
    u   = unpack(uflat,unPackMask,u)
    v   = unpack(vflat,unPackMask,v)
    p   = unpack(pflat,unPackMask,p)

    call PrimToCon(Q,rho,u,v,p)

    deallocate(uflat,vflat,rhoflat,pflat)
end subroutine ZeitounTestCaseConstantSlipBCLimiter

subroutine DGGradientBC(cU,gU,bcU,dUdx,dUdy)
    ! Computes the DG gradient, used when calculating the viscous stress tensor
    ! Arguments:    cU  -   Field at cubature points
    !               gu  -   Field at Gauss points
    !               bcu -   A copy of gU where boundary conditions have already been applied



    ! Author: Sam MacPherson
    implicit none

    ! define arguments
    real(dp), dimension(:,:), intent(in)    :: cU, gU, bcU
    real(dp), dimension(:,:), intent(out)   :: dUdx, dUdy

    ! Inverse mass matrix
    real(dp) , dimension(:,:), allocatable :: invM

    ! Interior and exterior Gauss nodes and fluxes
    real(dp), dimension(size(gU(:,1)),size(gU(1,:))) :: gUM, gUP, fx, fy
    ! Flat arrays
    real(dp), dimension(:), allocatable :: gUFlat, gUMflat, gUPflat, bcUflat

    ! Mask for unpacking variables
    logical, dimension(size(gU(:,1)),size(gU(1,:))) :: unPackMask
    unPackMask = .true.


    
    ! Calculate cubature part of gradient
    dUdx = matmul(transpose(cubDr),cubW*cubrx*cU) + matmul(transpose(cubDs),cubW*cubsx*cU)
    dUdy = matmul(transpose(cubDr),cubW*cubry*cU) + matmul(transpose(cubDs),cubW*cubsy*cU)
    
    
    
    ! Flatten matrices to be accessed by maps
    gUflat = pack(gU,.true.)
    bcUflat = pack(bcU,.true.)

    ! Access interior and exterior nodes
    gUMflat = gUflat(gaussmapM)
    gUPflat = gUflat(gaussmapP)
    
    ! Apply boundary conditions
    gUPflat(gaussmapB) = bcUflat(gaussmapB)
   
    ! Reform Matrices
    gUM = unpack(gUMflat,unPackMask,gUM)
    gUP = unpack(gUPflat,unPackMask,gUP)
  
    ! Calculate the centred flux terms
    fx = 0.5_dp*(gaussnx*(gUM+gUP))
    fy = 0.5_dp*(gaussny*(gUM+gUP))
    
    ! Add flux terms to gradient
    dUdx = dUdx - matmul(transpose(gaussinterp),gaussW*fx)
    dUdy = dUdy - matmul(transpose(gaussinterp),gaussW*fy)
    
    
    ! Apply inverse mass matrix and Jacobian
    invM = matmul(V,transpose(V))
    

    dUdx = dUdx/J
    dUdy = dUdy/J
    dUdx = matmul(invM,dUdx)
    dUdy = matmul(invM,dUdy)

    ! Correct the sign of the equation
    dUdx = -dUdx
    dUdy = -dUdy
    deallocate(gUFlat, gUMflat, gUPflat, bcUflat,invM)
end subroutine DGGradientBC

end Module