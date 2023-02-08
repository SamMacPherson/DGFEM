Module CNSSubRoutines
use GlobalVariables
use PrimConTransformation
use BCandIC
use Fluxes
use MatrixOperations

Contains

subroutine CNS2DRHS(Q,rhsQ)
    ! Calculates the RHS of the weak formulation of the semi-discrete equation
    ! Arguments:    Q       - Vector of conserved variables
    !               shapeQ  - Shape of Q (given as shape(Q))
    !               rhsQ    - The RHS of the equation

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(in)  :: Q
    real(dp), dimension(:,:,:), intent(out) :: rhsQ

    ! Conserved variables
    real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1))) :: rho, rhou, rhov, ener

    ! Define a temp Q for subroutine use, and Q vector for cubature and Gauss nodes
    real(dp), dimension(:,:,:), allocatable :: cQ, gQ, bQ
   
    ! Gauss definitions
    real(dp), dimension(:,:), allocatable :: grho, grhou, grhov, gener, brho, bu, bv, bp, gu, gv, gp, brhou, brhov, bener, gmu, gT, bT, gx, gy

    ! Cubature definitions
    real(dp), dimension(:,:), allocatable    :: crho, crhou, crhov, cener, cu, cv, cp, cmu, cT

    ! Derivative definitions
    real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1))) :: drhodx, drhody, drhoudx, drhoudy, drhovdx, drhovdy, dTdx, dTdy

    ! Derivative Gauss definitions
    real(dp), dimension(:,:), allocatable   :: gdrhodx, gdrhody, gdrhoudx, gdrhoudy, gdrhovdx, gdrhovdy, gdudx, gdvdx, gdudy, gdvdy, gdTdx, gdTdy

    ! Derivative cubature definitions
    real(dp), dimension(:,:), allocatable   :: cdrhodx, cdrhody, cdrhoudx, cdrhoudy, cdrhovdx, cdrhovdy, cdudx, cdvdx, cdudy, cdvdy, cdTdx, cdTdy

    ! Viscous stress tensor
    real(dp), dimension(:,:), allocatable   :: ct11, ct12, ct22, ct31, ct32, gt11, gt12, gt22, gt31, gt32

    ! Definitions for the call to divergence subroutine
    real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1)),4)    :: div, flux  
    real(dp), dimension(:,:), allocatable       ::  gW, bW, cUsubCon, cUsubVis,  cVsubCon, cVsubVis, gUsubVis, gVsubVis, gUsubCon, gVsubCon, bcUCon, bcVCon,bcUVis, bcVVis
    
    real(dp), dimension(3*NGauss,K,4)   ::  gXconvectiveFluxes, gYconvectiveFluxes, gXviscousFluxes, gYviscousFluxes
    real(dp), dimension(3*NGauss,K,4)   ::  bXconvectiveFluxes, bYconvectiveFluxes, bXviscousFluxes, bYviscousFluxes
    real(dp), dimension(Ncub,K,4)       ::  cXconvectiveFluxes, cYconvectiveFluxes, cXviscousFluxes, cYviscousFluxes
    integer:: i,j,count

    real::  start,finish
    
    ! I know it is slightly messy to have these as allocatabale arrays,
    ! but that is how I first wrote it and if I change it everything
    ! breaks, if it ain't broke...   (Sam MacPherson)
    allocate(crho(Ncub,K),crhou(Ncub,K),crhov(Ncub,K),cener(Ncub,K))
    allocate(cu(Ncub,K), cv(Ncub,K), cp(Ncub,K))
    allocate(cQ(Ncub,K,4))
    allocate(gQ(3*NGauss,K,4))
    allocate(bQ(3*NGauss,K,4))
    allocate(grho(3*Ngauss,K), grhou(3*Ngauss,K), grhov(3*Ngauss,K), gener(3*Ngauss,K), gx(3*Ngauss,K), gy(3*Ngauss,K))
    allocate(gu(3*Ngauss,K), gv(3*Ngauss,K), gp(3*Ngauss,K))
    allocate(brho(3*NGauss,k), bu(3*Ngauss,K), bv(3*Ngauss,K), bp(3*Ngauss,K), bT(3*Ngauss,K))
    allocate(brhou(3*Ngauss,K), brhov(3*Ngauss,K), bener(3*Ngauss,K))
    allocate(cT(Ncub,K), cmu(Ncub,K))
    allocate(gT(3*Ngauss,K), gmu(3*Ngauss,K))

   
    ! Extract conserved variables
    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    ! Generate cubature and Gauss node information
    ! Interpolate cubature nodes
    crho  = matmul(cubV,rho)
    crhou = matmul(cubV,rhou)
    crhov = matmul(cubV,rhov)
    cener = matmul(cubV,ener)


    cQ(:,:,1) = crho
    cQ(:,:,2) = crhou
    cQ(:,:,3) = crhov
    cQ(:,:,4) = cener


    ! Interpolate Gauss nodes
    grho  = matmul(gaussinterp,rho)
    grhou = matmul(gaussinterp,rhou)
    grhov = matmul(gaussinterp,rhov)
    gener = matmul(gaussinterp,ener)

    gQ(:,:,1) = grho
    gQ(:,:,2) = grhou
    gQ(:,:,3) = grhov
    gQ(:,:,4) = gener
    
    gx = matmul(gaussinterp,x)
    gy = matmul(gaussinterp,y)
    
    ! Primitive variables
    call ConToPrim(gQ,grho,gu,gv,gp)
    call ConToPrim(cQ,crho,cu,cv,cp)

    ! Calculate temperatures at Gauss and cubature nodes
    call GetTemperature(cp,crho,cT)
    call GetTemperature(gp,grho,gT)

    
    ! Apply BC's to boundary Gauss vectors
    bQ = gQ
    bT = gT
    

    select case(BCScheme)
        case(0)
            call ZeitounTestCaseNoSlipBC(bQ, bT, gaussmapW)
        case(1)
            call ZeitounTestCaseConstantSlipBC(bQ, bT, gaussmapW)
        case(2)
            call ZeitounTestCaseSlipBCGaussNodes(bQ, bT, cu, cT, gaussmapW)
        case(3)
            call ShearFlow2DBC(bQ,gaussmapB, gx, gy)
        case default
             print*,"Warning: Invalid boundary conditions have been chosen"
    end select

    brho  = bQ(:,:,1)
    brhou = bQ(:,:,2)
    brhov = bQ(:,:,3)
    bener = bQ(:,:,4)

    ! Calculate primitive variables at Gauss nodes with boundary data
    call ConToPrim(bQ,brho,bu,bv,bp)

    

   
    ! Calculate viscosity at Gauss and cubature nodes
    if (useSutherland) then
        call Sutherland(cT,cmu)
        call Sutherland(gT,gmu)
    else
        cmu = constMu
        gmu = constMu
    end if
  

    ! Short summary, now have made conservative and primitve cubature information
    ! Have conservative and primitve variables for Gauss nodes, with and without updated boundary conditions


    ! Calculate gradients
    call DGGradient(crho,grho,brho,drhodx,drhody)
    call DGGradient(crhou,grhou,brhou,drhoudx,drhoudy)
    call DGGradient(crhov,grhov,brhov,drhovdx,drhovdy)
    call DGGradient(cT,gT,bT,dTdx,dTdy)
    

    ! Interpolate derivatives at cubature nodes
    cdrhodx = matmul(cubV,drhodx)
    cdrhody = matmul(cubV,drhody)

    cdrhoudx = matmul(cubV,drhoudx)
    cdrhoudy = matmul(cubV,drhoudy)

    cdrhovdx = matmul(cubV,drhovdx)
    cdrhovdy = matmul(cubV,drhovdy)

    cdTdx = matmul(cubV,dTdx)
    cdTdy = matmul(cubV,dTdy)
    
    ! Use product rule to obtain velocity derivatives
    cdudx = (cdrhoudx - cdrhodx*cu)/crho
    cdudy = (cdrhoudy - cdrhody*cu)/crho

    cdvdx = (cdrhovdx - cdrhodx*cv)/crho
    cdvdy = (cdrhovdy - cdrhody*cv)/crho

  
    ! Compute viscous stress tensor at cubature nodes
    ct11 = cmu*(2.0_dp*cdudx - (2.0_dp/3.0_dp)*(cdudx+cdvdy))
    ct12 = cmu*(cdudy+cdvdx)
    ct22 = cmu*(2.0_dp*cdvdy - (2.0_dp/3.0_dp)*(cdudx+cdvdy))
    ct31 = cu*ct11 + cv*ct12
    ct32 = cu*ct12 + cv*ct22
  
    
    ! Interpolate derivatves at Gauss nodes
    gdrhodx = matmul(gaussinterp,drhodx)
    gdrhody = matmul(gaussinterp,drhody)

    gdrhoudx = matmul(gaussinterp,drhoudx)
    gdrhoudy = matmul(gaussinterp,drhoudy)

    gdrhovdx = matmul(gaussinterp,drhovdx)
    gdrhovdy = matmul(gaussinterp,drhovdy)

    gdTdx = matmul(gaussinterp,dTdx)
    gdTdy = matmul(gaussinterp,dTdy)

    

    ! Use product rule to obtain velocity derivatives
    gdudx = (gdrhoudx - gdrhodx*gu)/grho
    gdudy = (gdrhoudy - gdrhody*gu)/grho

    gdvdx = (gdrhovdx - gdrhodx*gv)/grho
    gdvdy = (gdrhovdy - gdrhody*gv)/grho

    ! Compute viscous stress tensor at cubature nodes
    gt11 = gmu*(2.0_dp*gdudx - (2.0_dp/3.0_dp)*(gdudx+gdvdy))
    gt12 = gmu*(gdudy+gdvdx)
    gt22 = gmu*(2.0_dp*gdvdy - (2.0_dp/3.0_dp)*(gdudx+gdvdy))
    gt31 = gu*gt11 + gv*gt12
    gt32 = gu*gt12 + gv*gt22
    
    
    ! Continuity equation
    ! Convective terms
    cXconvectiveFluxes(:,:,1) = crhou
    cYconvectiveFluxes(:,:,1) = crhov

    gXconvectiveFluxes(:,:,1) = grhou
    gYconvectiveFluxes(:,:,1) = grhov

    bXconvectiveFluxes(:,:,1) = brhou
    bYconvectiveFluxes(:,:,1) = brhov

    ! Viscous terms
    cXviscousFluxes(:,:,1) = crhou*0
    cYviscousFluxes(:,:,1) = crhov*0

    gXviscousFluxes(:,:,1) = grhou*0
    gYviscousFluxes(:,:,1) = grhov*0

    bXviscousFluxes(:,:,1) = brhou*0
    bYviscousFluxes(:,:,1) = brhov*0
    

    ! x-momentum equation
    ! Convective terms
    cXconvectiveFluxes(:,:,2) = crhou*cu + cp
    cYconvectiveFluxes(:,:,2) = crhov*cu

    gXconvectiveFluxes(:,:,2) = grhou*gu + gp
    gYconvectiveFluxes(:,:,2) = grhov*gu

    bXconvectiveFluxes(:,:,2) = brhou*bu + bp
    bYconvectiveFluxes(:,:,2) = brhov*bu

    ! Viscous terms
    cXviscousFluxes(:,:,2) = ct11
    cYviscousFluxes(:,:,2) = ct12

    gXviscousFluxes(:,:,2) = gt11
    gYviscousFluxes(:,:,2) = gt12

    bXviscousFluxes(:,:,2) = gt11
    bYviscousFluxes(:,:,2) = gt12

    ! y-momentum equation
    ! Convective terms
    cXconvectiveFluxes(:,:,3) = crhou*cv
    cYconvectiveFluxes(:,:,3) = crhov*cv + cp 

    gXconvectiveFluxes(:,:,3) = grhou*gv
    gYconvectiveFluxes(:,:,3) = grhov*gv + gp 

    bXconvectiveFluxes(:,:,3) = brhou*bv
    bYconvectiveFluxes(:,:,3) = brhov*bv + bp 

    ! Viscous terms
    cXviscousFluxes(:,:,3) = ct12
    cYviscousFluxes(:,:,3) = ct22

    gXviscousFluxes(:,:,3) = gt12
    gYviscousFluxes(:,:,3) = gt22

    bXviscousFluxes(:,:,3) = gt12
    bYviscousFluxes(:,:,3) = gt22


    

    ! energy equation
    ! Convective terms
    cXconvectiveFluxes(:,:,4) = (cener + cp)*cu
    cYconvectiveFluxes(:,:,4) = (cener + cp)*cv 

    gXconvectiveFluxes(:,:,4) = (gener + gp)*gu
    gYconvectiveFluxes(:,:,4) = (gener + gp)*gv 

    bXconvectiveFluxes(:,:,4) = (bener + bp)*bu
    bYconvectiveFluxes(:,:,4) = (bener + bp)*bv 

    ! Viscous terms
    cXviscousFluxes(:,:,4) = ct31 + thermalConductivity*cdTdx
    cYviscousFluxes(:,:,4) = ct32 + thermalConductivity*cdTdy

    gXviscousFluxes(:,:,4) = gt31 + thermalConductivity*gdTdx
    gYviscousFluxes(:,:,4) = gt32 + thermalConductivity*gdTdy

    bXviscousFluxes(:,:,4) = gt31 + thermalConductivity*gdTdx
    bYviscousFluxes(:,:,4) = gt32 + thermalConductivity*gdTdy

    ! Calculate cubature terms
    call DGCubatureTerm(cXconvectiveFluxes,cXviscousFluxes,cYconvectiveFluxes,cYviscousFluxes,div)
    rhsQ = div

    ! Calculate flux terms
    call DGFluxTerm(gQ,bQ,gXconvectiveFluxes,gXviscousFluxes,gYconvectiveFluxes,gYviscousFluxes, &
                          bXconvectiveFluxes,bXviscousFluxes,bYconvectiveFluxes,bYviscousFluxes,flux)

    rhsQ = div - flux
    deallocate(grho, grhou, grhov, gener, brho, bu, bv, bp, gu, gv, gp, brhou, brhov, bener, gmu, gT, gx, gy)
    deallocate(crho, crhou, crhov, cener, cu, cv, cp, cmu, cT)
    deallocate(gdrhodx, gdrhody, gdrhoudx, gdrhoudy, gdrhovdx, gdrhovdy, gdudx, gdvdx, gdudy, gdvdy)
    deallocate(cdrhodx, cdrhody, cdrhoudx, cdrhoudy, cdrhovdx, cdrhovdy, cdudx, cdvdx, cdudy, cdvdy)
    deallocate(ct11, ct12, ct22, ct31, ct32, gt11, gt12, gt22, gt31, gt32)

end subroutine CNS2DRHS

subroutine DGCubatureTerm(cXconvectiveFluxes,cXviscousFluxes,cYconvectiveFluxes,cYviscousFluxes,cubResult)
    ! Computes the DG Cubature term, used when calculating the weak form of the DG-FEM
    ! This function works for one of the four equations which define the CNS equations
    ! Arguments:    cXconvectiveFluxes  -   x-derivative convective fluxes at cubature nodes
    !               cXviscousFluxes     -   x-derivative viscous fluxes at cubature nodes
    !               cYconvectiveFluxes  -   y-derivative convective fluxes at cubature nodes
    !               cYviscousFluxes     -   y-derivative viscous fluxes at cubature nodes    
    ! Author: Sam MacPherson

    implicit none

    ! define arguments
    real(dp), dimension(:,:,:), intent(in)  ::  cXviscousFluxes, cXconvectiveFluxes, cYconvectiveFluxes, cYviscousFluxes
    real(dp), dimension(:,:,:), intent(out) ::  cubResult
    real(dp), dimension(:,:), allocatable   ::  div

    integer ::  i

    do i = 1,4
        ! Cubature volume terms
        div = matmul(transpose(cubDr), cubW*(cubrx*(cXconvectiveFluxes(:,:,i) - cXviscousFluxes(:,:,i))   & 
            + cubry*(cYconvectiveFluxes(:,:,i) - cYviscousFluxes(:,:,i))) )                               &
            + matmul(transpose(cubDs), cubW*(cubsx*(cXconvectiveFluxes(:,:,i) - cXviscousFluxes(:,:,i) )  & 
            + cubsy*(cYconvectiveFluxes(:,:,i) - cYviscousFluxes(:,:,i) )) )

        ! Apply inverse mass matrix
        div = matmul(transpose(V),div)
        div= matmul(V,div)
        div = div/J

        cubResult(:,:,i) = div
    end do

    deallocate(div)
end subroutine DGCubatureTerm

subroutine DGFluxTerm(gQ,bQ,gXconvectiveFluxes,gXviscousFluxes,gYconvectiveFluxes,gYviscousFluxes,bXconvectiveFluxes,bXviscousFluxes,bYconvectiveFluxes,bYviscousFluxes,fluxResult)
    ! Calculates the flux term for the DG scheme
    ! Arguments:    gQ                  -   Conserved vector at Gauss nodes
    !               bQ                  -   Conserved vector with bc's applied at Gauss nodes
    !               gXconvectiveFluxes  -   x-derivative convective fluxes at Gauss nodes
    !               gXviscousFluxes     -   x-derivative viscous fluxes at Gauss nodes
    !               gYconvectiveFluxes  -   y-derivative convective fluxes at Gauss nodes
    !               gYviscousFluxes     -   y-derivative viscous fluxes at Gauss nodes  
    !               bXconvectiveFluxes  -   x-derivative convective fluxes at Gauss boundary nodes
    !               bXviscousFluxes     -   x-derivative viscous fluxes at Gauss boundary nodes
    !               bYconvectiveFluxes  -   y-derivative convective fluxes at Gauss boundary nodes
    !               bYviscousFluxes     -   y-derivative viscous fluxes at Gauss boundary nodes     

    implicit none
    real(dp), dimension(:,:,:), intent(in)  ::  gQ, bQ
    real(dp), dimension(:,:,:), intent(in)  ::  gXviscousFluxes, gXconvectiveFluxes, gYconvectiveFluxes, gYviscousFluxes
    real(dp), dimension(:,:,:), intent(in)  ::  bXviscousFluxes, bXconvectiveFluxes, bYconvectiveFluxes, bYviscousFluxes
    real(dp), dimension(:,:,:), intent(out) ::  fluxResult
    real(dp), dimension(3*NGauss,K,4)             ::  fluxCon, fluxVis, fluxTemp

    ! Flat arrays
    real(dp), dimension(:), allocatable :: gUConflat, gUMConflat, gUPConflat, bcUConflat
    real(dp), dimension(:), allocatable :: gUVisflat, gUMVisflat, gUPVisflat, bcUVisflat

    real(dp), dimension(:), allocatable :: gVConflat, gVMConflat, gVPConflat, bcVConflat, bcNdotUConflat
    real(dp), dimension(:), allocatable :: gVVisflat, gVMVisflat, gVPVisflat, bcVVisflat, bcNdotUVisflat
    real(dp), dimension(:), allocatable :: gaussnxflat, gaussnyflat

    ! Define normal flux arrays
    real(dp), dimension(:), allocatable     :: gFxMConflat, gFxPConflat
    real(dp), dimension(:,:), allocatable   :: gFxMCon, gFxPCon

    real(dp), dimension(:), allocatable     :: gFxMVisflat, gFxPVisflat
    real(dp), dimension(:,:), allocatable   :: gFxMVis, gFxPVis

    real(dp), dimension(3*NGauss,K,4)       :: normalFluxConM, normalFluxConP, normalFluxVisM, normalFluxVisP

    ! Interior and exterior nodes
    real(dp), dimension(3*NGauss,K)     ::  gWM, gWP
    real(dp), dimension(3*NGauss,K,4)   ::  gQM, gQP
    real(dp), dimension(:), allocatable ::  gWflat, bWflat, gWMflat, gWPflat

    logical, dimension(3*NGauss,K)  :: unPackMask

    integer ::  i

    allocate(gFxMCon(3*NGauss,K), gFxPCon(3*NGauss,K))
    allocate(gFxMVis(3*NGauss,K), gFxPVis(3*NGauss,K))

    unPackMask = .true.

    do i = 1,4
        ! Interior (left) and exterior (right) states at each Gauss node
        gWflat = pack(gQ(:,:,i),.true.)
        bWflat = pack(bQ(:,:,i),.true.)
        gWMflat = gWflat(gaussmapM)
        gWPflat = gWflat(gaussmapP)
    
        ! Apply BC's
        gWPflat(gaussmapB) = bWflat(gaussmapB)
        
        ! Reform matrices
        gWM = unpack(gWMflat,unPackMask,gWM)
        gWP = unpack(gWPflat,unPackMask,gWP)
        gQM(:,:,i) = gWM
        gQP(:,:,i) = gWP
    

    
        ! Flatten convective flux matrices to be accessed by maps
        gUConflat      = pack(gXconvectiveFluxes(:,:,i),.true.)
        gVConflat      = pack(gYconvectiveFluxes(:,:,i),.true.)
        bcUConflat     = pack(bXconvectiveFluxes(:,:,i),.true.)
        bcVConflat     = pack(bYconvectiveFluxes(:,:,i),.true.)

        ! Flatten viscous flux matrices to be accessed by maps
        gUVisflat      = pack(gXviscousFluxes(:,:,i),.true.)
        gVVisflat      = pack(gYviscousFluxes(:,:,i),.true.)
        bcUVisflat     = pack(bXviscousFluxes(:,:,i),.true.)
        bcVVisflat     = pack(bYviscousFluxes(:,:,i),.true.)

        gaussnxflat = pack(gaussnx,.true.)
        gaussnyflat = pack(gaussny,.true.)
        
        ! Access convective fluxes at interior and exterior nodes
        gUMConflat = gUConflat(gaussmapM)
        gUPConflat = gUConflat(gaussmapP)
        gVMConflat = gVConflat(gaussmapM)
        gVPConflat = gVConflat(gaussmapP)

        ! Access viscous fluxes at interior and exterior nodes
        gUMVisflat = gUVisflat(gaussmapM)
        gUPVisflat = gUVisflat(gaussmapP)
        gVMVisflat = gVVisflat(gaussmapM)
        gVPVisflat = gVVisflat(gaussmapP)

        ! Resolve convective fluxes along face normal on interior (left) and exterior (right) states
        gFxMConflat = gaussnxflat*gUMConflat + gaussnyflat*gVMConflat
        gFxPConflat = gaussnxflat*gUPConflat + gaussnyflat*gVPConflat

        ! Resolve viscous fluxes along face normal on interior (left) and exterior (right) states
        gFxMVisflat = gaussnxflat*gUMVisflat + gaussnyflat*gVMVisflat
        gFxPVisflat = gaussnxflat*gUPVisflat + gaussnyflat*gVPVisflat

        ! Convective normal flux for the boundary values
        bcNdotUConflat = gaussnxflat*bcUConflat + gaussnyflat*bcVConflat

        ! Viscous normal flux for the boundary values
        bcNdotUVisflat = gaussnxflat*bcUVisflat + gaussnyflat*bcVVisflat

        ! Apply BC's
        gFxPConflat(gaussmapB) = bcNdotUConflat(gaussmapB)
        gFxPVisflat(gaussmapB) = bcNdotUVisflat(gaussmapB)
    
        ! Reform convective matrices
        gFxMCon = unpack(gFxMConflat,unPackMask,gFxMCon)
        gFxPCon = unpack(gFxPConflat,unPackMask,gFxPCon)

        ! Reform viscous matrices
        gFxMVis = unpack(gFxMVisflat,unPackMask,gFxMVis)
        gFxPVis = unpack(gFxPVisflat,unPackMask,gFxPVis)

        ! Place normal fluxes into normal flux arrays
        normalFluxConM(:,:,i) = gFxMCon
        normalFluxVisM(:,:,i) = gFxMVis

        normalFluxConP(:,:,i) = gFxPCon
        normalFluxVisP(:,:,i) = gFxPVis

    end do

    select case(fluxScheme)
        ! 0 - Lax-Friedrich
        ! 1 - HLL
        ! 2 - HLLC
        case(0)
            ! Both viscous and convective terms can be used with Lax-Friedrich
            call LaxFriedrichFlux( normalFluxConM , normalFluxConP, gQM, gQP,fluxCon)
            call CenteredFlux( -normalFluxVisM , -normalFluxVisP, fluxVis)
            fluxTemp = fluxCon + fluxVis
         
        case(1)
            ! Use centred flux for viscous terms
            call CenteredFlux(-normalFluxVisM,-normalFluxVisP,fluxVis)
            
            ! Use HLL flux for convective terms
            call HLLFlux(normalFluxConM,normalFluxConP,gQM,gQP,fluxCon)

            ! Total flux is convective minis viscous fluxes
            fluxTemp = fluxCon + fluxVis

        case(2)
            ! Use centred flux for viscous terms
            call CenteredFlux(-normalFluxVisM,-normalFluxVisP,fluxVis)
            
            ! Use HLLC flux for convective terms
            call HLLCFlux(normalFluxConM,normalFluxConP,gQM,gQP,fluxCon)

            ! Total flux is convective minis viscous fluxes
            fluxTemp = fluxCon + fluxVis
     
        case default
             print*,"Warning: Invalid flux scheme has been chosen"

    end select
    do i = 1,4
        fluxResult(:,:,i) = matmul(transpose(gaussinterp),gaussW*fluxTemp(:,:,i))

        fluxResult(:,:,i) = matmul(transpose(V),fluxResult(:,:,i))
        fluxResult(:,:,i) = matmul(V,fluxResult(:,:,i))
        fluxResult(:,:,i) = fluxResult(:,:,i)/J
    end do
end subroutine DGFluxTerm

subroutine DGGradient(cU,gU,bcU,dUdx,dUdy)
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
    ! dUdx = matmul(invM,dUdx)/J
    ! dUdy = matmul(invM,dUdx)/J

    dUdx = dUdx/J
    dUdy = dUdy/J
    dUdx = matmul(invM,dUdx)
    dUdy = matmul(invM,dUdy)

    ! call MatrixMultiply(invM,dUdx,dUdx)
    ! call MatrixMultiply(invM,dUdy,dUdy)

    
    ! Correct the sign of the equation
    dUdx = -dUdx
    dUdy = -dUdy
    deallocate(gUFlat, gUMflat, gUPflat, bcUflat,invM)
end subroutine DGGradient

subroutine CNSdt2D(Q,dt)
    ! Calculates the time step for the CNS solver
    ! Arguments:    Q       -   Vector of conserved variables
    !               dt      -   The stable time step 

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(in)  :: Q
    real(dp), intent(out)                   :: dt

    ! Primitive variables and Misc.
    real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1))) :: rho, u, v, p, c, mu, T
    real(dp), dimension(:), allocatable                 :: h
    ! Flattened primitive variables and Misc 
    real(dp), dimension(:), allocatable :: uflat, vflat, pflat, cflat, lam, Fscaleflat, muflat

    ! Extract primitive variables
    call ConToPrim(Q,rho,u,v,p)
    call GetTemperature(p,rho,T)

    if (useSutherland) then
        call Sutherland(T,mu)
    else
        mu = constMu
    end if

    ! Calculate sound speed
    call SoundSpeed(p,rho,c)

    ! As Fscale is calculated at surface nodes, need to extract surface node information
    ! To use the surface map - vmapM - in Fortran the arrays need to first be flattened
    uflat = pack(u,.true.)
    vflat = pack(v,.true.)
    pflat = pack(p,.true.)
    cflat = pack(c,.true.)
    muflat = pack(mu,.true.)
    ! Access the surface nodes
    uflat = uflat(vmapM)
    vflat = vflat(vmapM)
    pflat = pflat(vmapM)
    cflat = cflat(vmapM)
    muflat = muflat(vmapM)

    lam = sqrt(uflat**2.0_dp+vflat**2.0_dp) + cflat
    

    Fscaleflat = pack(Fscale,.true.)
    
    h = 2.0_dp/Fscaleflat
    dt = 0.5_dp*minval( 1.0_dp/( (N+1.0_dp)**2.0_dp*lam/h + (N+1.0_dp)**4.0_dp*muflat/h/h ))

    ! Deallocate the "allocatable" arrays
    deallocate(uflat,vflat,pflat,cflat,lam,Fscaleflat,h,muflat)
end subroutine CNSdt2D

subroutine Limiter2D(Q)
    ! Implements the limiter of Tu (2005)
    ! Arguments:    Q   -   Vector of conserved variables

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(Np,K,4), intent(out)       ::  Q

    ! Conserved variables
    real(dp), dimension(Np,K)                       :: rho, rhou, rhov, ener, p
    real(dp), dimension(Np*K)                       :: rhoflat, rhouflat, rhovflat, enerflat
    
    ! unpack matrix
    logical, dimension(:,:), allocatable            :: unPackMask
    ! Average matrix
    real(dp), dimension(Np)                     :: aveflat
    real(dp), dimension(1,Np)                   :: ave
    real(dp), dimension(Np,Np)                  :: dropAve

    ! Definitions for identity and one matrices
    real(dp), dimension(Np,Np)                  :: eye
    real(dp), dimension(:,:), allocatable       :: ones

    ! Distances to centre of the element
    real(dp), dimension(Np,K)  :: dx, dy

    ! Properties for patch and vertex extraction
    integer, dimension(K)                           :: E1, E2, E3, v1, v2, v3
    real(dp), dimension(3,K)                        :: xv,yv 

    ! Face normals
    real(dp), dimension(3,K)                        :: fnx,fny,fLen                           
    
    ! Element centers
    real(dp), dimension(1,K)                        :: xc0, yc0
    real(dp), dimension(3,K)                        :: xc, yc

    ! Area weights
    real(dp), dimension(1,K)                        :: A0
    real(dp), dimension(3,K)                        :: A

    ! Boundary id's
    integer, dimension(:), allocatable              :: id1, id2, id3

    ! Used for ghost elements
    real(dp), dimension(:), allocatable             :: H1, H2, H3

    ! Loop integers
    integer                                         :: i,n

    ! Average value at center of elements
    real(dp), dimension(1,K)                        :: rhoC, rhouC, rhovC, enerC
    real(dp), dimension(Np,K)                       :: averRho, averRhou, averRhov, averEner

    ! Primitive variables at element centers
    real(dp), dimension(1,K,4)                      :: PC0

    ! Primitive variables neighbours
    real(dp), dimension(3,K,4)                      :: PC
    real(dp), dimension(:), allocatable             :: PCflat
    logical, dimension(3,K)                         :: PCUnpackMask
    real(dp), dimension(3,K)                        :: T    
    ! Boundary condition markers
    integer, dimension(:), allocatable              :: idW, idMW, idI, idO, idS, idD, idN, idB
    integer, dimension(:), allocatable              :: BCTypeFlat

    integer, dimension(6)                           :: idV

    ! Vertex map
    integer, dimension(6*K)                         :: vmapMVertex, vmapPVertex

    ! Vertex array
    real(dp), dimension(6,K,4)                      :: PVVertex
    real(dp), dimension(6,K)                        :: rhoVertex, rhouVertex, rhovVertex, enerVertex
    real(dp), dimension(6*K)                        :: rhoVertexFlat,rhouVertexFlat,rhovVertexFlat,enerVertexFlat

    ! Unlimited gradients
    real(dp), dimension(3,K)                        :: dVdxFace, dVdyFace, dVdxNeigh, dVdyNeigh
    real(dp), dimension(1,K)                        :: dVdxC0, dVdyC0
    real(dp), dimension(3,K)                        :: dVdxC, dVdyC

    ! Limited gradient
    real(dp), dimension(1,K)                        :: LdVdx, LdVdy

    ! Matrices used in calculation of unlimited gradients
    real(dp), dimension(1,K)                        :: VC0
    real(dp), dimension(3,K)                        :: VC
    real(dp), dimension(6,K)                        :: VVertex

    ! Weights for limited gradient
    real(dp), dimension(3,K)                        :: g, w
    real(dp), parameter                             :: epsilon = 1e-10
  
    ! Taylor expansion
    real(dp), dimension(Np,K,4)                       :: dV, aV

    ! Matrix for All function in negative pressure correction
    real(dp), dimension(Np) :: zeros

    ! Array to hold min values
    real(dp), dimension(K)  :: minArrayDensity, minArrayPressure

    ! Integer arrays to hold elements that need density and pressure correction
    integer, dimension(:), allocatable  :: idRhoNeg, idPNeg


    real(dp), dimension(3*NGauss,K) ::  bu, bv, bp, brho, bT, bener
    real(dp), dimension(Ncub,K)     ::  cu, cv, cp, crho, cT, cener

    PCUnpackMask = .true.
    ! Note: Most arrays are stored as a 1 x n matrix so they can be matrix multiplied,  
    ! or be used with arrays that need to be matrix multiplied

    ! Construct a vector which is 1/3 on nodes which are a vertex, and 0 otherwise
    zeros = 0.0_dp
    
    aveflat = sum(massMatrix,dim=1)/2.0_dp
    ave(1,:) = aveflat

    ! Identity matrix
    call IdentityMatrix(Np,eye)

    ! Matric of ones
    allocate(ones(Np,1))
    ones = 1.0_dp

    ! Calculate the distance of vertices to cell center
    dropAve = eye - matmul(ones,ave)
    dx = matmul(dropAve,x)
    dy = matmul(dropAve,y)
    
    ! Neighbours for the patch
    E1 = EToE(:,1)
    E2 = EToE(:,2)
    E3 = EToE(:,3)

    ! Extract the vertices
    v1 = EToV(:,1)
    v2 = EToV(:,2)
    v3 = EToV(:,3)
    
    xv(1,:) = VX(v1)
    xv(2,:) = VX(v2)
    xv(3,:) = VX(v3)

    yv(1,:) = VY(v1)
    yv(2,:) = VY(v2)
    yv(3,:) = VY(v3)
   
    
    ! Face normals
    fnx(1,:) = yv(2,:) - yv(1,:)
    fnx(2,:) = yv(3,:) - yv(2,:)
    fnx(3,:) = yv(1,:) - yv(3,:)

    fny(1,:) = -1*(xv(2,:) - xv(1,:))
    fny(2,:) = -1*(xv(3,:) - xv(2,:))
    fny(3,:) = -1*(xv(1,:) - xv(3,:))

    
    ! Normalise face normals
    fLen = sqrt(fnx**2.0_dp+fny**2.0_dp)
    fnx = fnx/fLen
    fny = fny/fLen

    ! Calculate element centers by averaging the vertex values
    xc0 = matmul(ave,x)
    yc0 = matmul(ave,y)

    xc(1,:) = xc0(1,E1)
    xc(2,:) = xc0(1,E2)
    xc(3,:) = xc0(1,E3)

    yc(1,:) = yc0(1,E1)
    yc(2,:) = yc0(1,E2)
    yc(3,:) = yc0(1,E3)
    
   
    ! The area of each element is matmul(ave,J)*2
    ! This is because the straight sided triangle has area 2, 
    ! and so the real area is the Jacobian multiplied by this
    A0 = matmul(ave,J)*2.0_dp/3.0_dp
    A(1,:) = A0(1,:) + A0(1,E1)
    A(2,:) = A0(1,:) + A0(1,E2)
    A(3,:) = A0(1,:) + A0(1,E3)
    

    allocate(id1(0),id2(0),id3(0))
   
    ! Find boundary faces
    do i = 1,size(BCType(:,1))
        if (BCType(i,1).ne.0) then
            id1 = [id1, i]
        end if
        if (BCType(i,2).ne.0) then
            id2 = [id2, i]
        end if
        if (BCType(i,3).ne.0) then
            id3 = [id3, i]
        end if
    end do
   
    ! Location of center of ghost elements
    H1 = 2.0_dp*(A0(1,id1)/fLen(1,id1))
    H2 = 2.0_dp*(A0(1,id2)/fLen(2,id2))
    H3 = 2.0_dp*(A0(1,id3)/fLen(3,id3))

    xc(1,id1) = xc(1,id1) + 2.0_dp*fnx(1,id1)*H1
    xc(2,id2) = xc(2,id2) + 2.0_dp*fnx(2,id2)*H2
    xc(3,id3) = xc(3,id3) + 2.0_dp*fnx(3,id3)*H3

    yc(1,id1) = yc(1,id1) + 2.0_dp*fny(1,id1)*H1
    yc(2,id2) = yc(2,id2) + 2.0_dp*fny(2,id2)*H2
    yc(3,id3) = yc(3,id3) + 2.0_dp*fny(3,id3)*H3
    
   
  
    
    ! Extract conserved variables
    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)

    cu = matmul(cubV,rhou/rho)
    bu = matmul(gaussinterp,rhou/rho)

    cv = matmul(cubV,rhov/rho)
    bv = matmul(gaussinterp,rhov/rho)

    cener = matmul(cubV,ener)
    bener = matmul(gaussinterp,ener)

    crho = matmul(cubV,rho)
    brho = matmul(gaussinterp,rho)

    cT = (gammaGas - 1.0_dp)*(cener - 0.5_dp*(cu**2.0_dp + cv**2.0_dp))
    bT = (gammaGas - 1.0_dp)*(bener - 0.5_dp*(bu**2.0_dp + bv**2.0_dp))

    
    
    ! Construct flat arrays for later use
    rhoflat  = pack(rho,.true.)
    rhouflat = pack(rhou,.true.)
    rhovflat = pack(rhov,.true.)
    enerflat = pack(ener,.true.)

    ! Cell averages of conserved variables
    rhoC  = matmul(ave,rho)
    rhouC = matmul(ave,rhou)
    rhovC = matmul(ave,rhov)
    enerC = matmul(ave,ener)
    
    ! Place cell average information at all nodal points
    averRho  = matmul(ones,rhoC)
    averRhou = matmul(ones,rhouC)
    averRhov = matmul(ones,rhovC)
    averEner = matmul(ones,enerC)


    ! Primitive variables at element centers
    PC0(:,:,1) = rhoC
    PC0(:,:,2) = rhouC/rhoC
    PC0(:,:,3) = rhovC/rhoc
    PC0(:,:,4) = (gammaGas - 1.0_dp)*(enerC - 0.5_dp*(rhouC**2.0_dp+rhovC**2.0_dp)/rhoC)

    ! print*,enerc(1,2463),rhouC(1,2463),rhovC(1,2463),rhoC(1,2463)
    ! print*,PC0(1,2463,4)
    
    ! An array for primitive variable neighbours
    ! Currently in conservative form to be used with BC's
    ! Will be converted after BC's applied
    PC(1,:,1) = rhoC(1,EToE(:,1))
    PC(2,:,1) = rhoC(1,EToE(:,2))
    PC(3,:,1) = rhoC(1,EToE(:,3))

    PC(1,:,2) = rhouC(1,EToE(:,1))
    PC(2,:,2) = rhouC(1,EToE(:,2))
    PC(3,:,2) = rhouC(1,EToE(:,3))

    PC(1,:,3) = rhovC(1,EToE(:,1))
    PC(2,:,3) = rhovC(1,EToE(:,2))
    PC(3,:,3) = rhovC(1,EToE(:,3))

    PC(1,:,4) = enerC(1,EToE(:,1))
    PC(2,:,4) = enerC(1,EToE(:,2))
    PC(3,:,4) = enerC(1,EToE(:,3))
    
  
    T = (  (gammaGas-1.0_dp)*( PC(:,:,4) - 0.5_dp* (PC(:,:,2)**2.0_dp + PC(:,:,3)**2.0_dp)/PC(:,:,1) )  ) /GasConstant/PC(:,:,1)
    ! Find and classify boundary faces

    allocate(idW(0),idMW(0),idI(0),idO(0),idN(0),idD(0),idS(0),idB(0))
    BCTypeFlat = pack(transpose(BCType),.true.)
    
    ! This do loop performs the same function as find in Matlab
    ! Need a flattened BCType array to perform in Fortran
    do i = 1,size(BCTypeFlat)
        if(BCTypeFlat(i).eq.Inflow) then
            idI = [idI,i]
        end if

        if(BCTypeFlat(i).eq.Outflow) then
            idO = [idO,i]
        end if

        if(BCTypeFlat(i).eq.Wall) then
            idW = [idW,i]
        end if

        if(BCTypeFlat(i).eq.Moving) then
            idMW = [idMW,i]
        end if

        if(BCTypeFlat(i).eq.Slip) then
            idS = [idS,i]
        end if

        if(BCTypeFlat(i).eq.Dirichlet) then
            idD = [idD,i]
        end if

        if(BCTypeFlat(i).eq.Neumann) then
            idN = [idN,i]
        end if

        if(BCTypeFlat(i).ne.0) then
            idB = [idB, i]
        end if
    end do

 
   


    ! Apply BC's
    select case(BCScheme)
        case(0)
            call ZeitounTestCaseNoSlipBCLimiter(PC,T,idW)
        case(1)
            call ZeitounTestCaseConstantSlipBCLimiter(PC,T,idW)
        case(2)
            call ZeitounTestCaseSlipBCLimiter(PC, T, cu, bu, cT, bT, brho, idW)
        case(3)
            call ShearFlow2DBC(PC,idB,xc,yc)
        case default
             print*,"Warning: Invalid boundary conditions have been chosen"
    end select

    
    ! Convert to primitive variables
    PC(:,:,2) = PC(:,:,2)/PC(:,:,1)
    PC(:,:,3) = PC(:,:,3)/PC(:,:,1)
    PC(:,:,4) = (gammaGas-1.0_dp)*( PC(:,:,4)-0.5_dp*PC(:,:,1)*(PC(:,:,2)**2.0_dp + PC(:,:,3)**2.0_dp) )

   

    ! Extract the vertex nodes from the face nodes
    idV = [1, Nfp, Nfp+1, 2*Nfp, 2*Nfp+1, 3*Nfp]
   
    ! For each element, pick out vertex nodes from vmapM and vmapP
    ! To form a vmap of only vertex nodes
    do i = 1,K
        vmapMVertex(6*i-5:6*i) = vmapM([3*Nfp*(i-1) + 1 , 3*Nfp*(i-1) + Nfp , 3*Nfp*(i-1) + Nfp + 1, 3*Nfp*(i-1) + 2*Nfp, 3*Nfp*(i-1) + 3*Nfp, 3*Nfp*(i-1) +2*Nfp+1])
        vmapPVertex(6*i-5:6*i) = vmapP([3*Nfp*(i-1) + 1 , 3*Nfp*(i-1) + Nfp , 3*Nfp*(i-1) + Nfp + 1, 3*Nfp*(i-1) + 2*Nfp, 3*Nfp*(i-1) + 3*Nfp, 3*Nfp*(i-1) +2*Nfp+1])
    end do
  
   

    ! Take an average of exterior and interior nodes to calculate the value of values at vertices
    rhoVertexFlat  = (rhoflat(vmapMVertex)  + rhoflat(vmapPVertex ))/2.0_dp
    rhouVertexFlat = (rhouflat(vmapMVertex) + rhouflat(vmapPVertex))/2.0_dp
    rhovVertexFlat = (rhovflat(vmapMVertex) + rhovflat(vmapPVertex))/2.0_dp
    enerVertexFlat = (enerflat(vmapMVertex) + enerflat(vmapPVertex))/2.0_dp

    ! Reform matrices
    allocate(unPackMask(6,K))
    unPackMask = .true.

    rhoVertex  = unpack(rhoVertexFlat,unPackMask,rhoVertex)
    rhouVertex = unpack(rhouVertexFlat,unPackMask,rhouVertex)
    rhovVertex = unpack(rhovVertexFlat,unPackMask,rhovVertex)
    enerVertex = unpack(enerVertexFlat,unPackMask,enerVertex)

    

    ! Construct an array of primitive variables at vertices
    PVVertex(:,:,1) = rhoVertex
    PVVertex(:,:,2) = rhouVertex/rhoVertex
    PVVertex(:,:,3) = rhovVertex/rhoVertex
    PVVertex(:,:,4) = (gammaGas -1.0_dp)*(enerVertex - 0.5_dp*rhoVertex*(PVVertex(:,:,2)**2.0_dp + PVVertex(:,:,3)**2.0_dp))

   
    
    ! Use do loop to calculate unlimited gradient
    do i = 1,4
        VC0     = PC0(:,:,i)
        VC      = PC(:,:,i)
        VVertex = PVVertex(:,:,i)
  
        ! Calculate face gradients
        dVdxFace(1,:) =  0.5_dp * ( (VC(1,:)-VC0(1,:) )*(yv(2,:)-yv(1,:)) + (VVertex(1,:)-VVertex(2,:))*(yc(1,:)-yc0(1,:)) )/A(1,:)
        dVdyFace(1,:) = -0.5_dp * ( (VC(1,:)-VC0(1,:) )*(xv(2,:)-xv(1,:)) + (VVertex(1,:)-VVertex(2,:))*(xc(1,:)-xc0(1,:)) )/A(1,:)

        dVdxFace(2,:) =  0.5_dp * ( (VC(2,:)-VC0(1,:) )*(yv(3,:)-yv(2,:)) + (VVertex(3,:)-VVertex(4,:))*(yc(2,:)-yc0(1,:)) )/A(2,:)
        dVdyFace(2,:) = -0.5_dp * ( (VC(2,:)-VC0(1,:) )*(xv(3,:)-xv(2,:)) + (VVertex(3,:)-VVertex(4,:))*(xc(2,:)-xc0(1,:)) )/A(2,:)

        dVdxFace(3,:) =  0.5_dp * ( (VC(3,:)-VC0(1,:) )*(yv(1,:)-yv(3,:)) + (VVertex(5,:)-VVertex(6,:))*(yc(3,:)-yc0(1,:)) )/A(3,:)
        dVdyFace(3,:) = -0.5_dp * ( (VC(3,:)-VC0(1,:) )*(xv(1,:)-xv(3,:)) + (VVertex(5,:)-VVertex(6,:))*(xc(3,:)-xc0(1,:)) )/A(3,:)

        
        ! Area weighted averaged to calculate cell gradient
        dVdxC0(1,:) = ( A(1,:)*dVdxFace(1,:) + A(2,:)*dVdxFace(2,:) + A(3,:)*dVdxFace(3,:) )/(A(1,:)+A(2,:)+A(3,:))

        dVdyC0(1,:) = ( A(1,:)*dVdyFace(1,:) + A(2,:)*dVdyFace(2,:) + A(3,:)*dVdyFace(3,:) )/(A(1,:)+A(2,:)+A(3,:))
        
        ! Neighour gradients
        dVdxC(1,:) = dVdxc0(1,E1)
        dVdyC(1,:) = dVdyc0(1,E1)

        dVdxC(2,:) = dVdxc0(1,E2)
        dVdyC(2,:) = dVdyc0(1,E2)

        dVdxC(3,:) = dVdxc0(1,E3)
        dVdyC(3,:) = dVdyc0(1,E3)
        
        ! Use face gradients for ghost elements
        dVdxC(1,id1) = dVdxFace(1,id1)
        dVdyC(1,id1) = dVdyFace(1,id1)

        dVdxC(2,id2) = dVdxFace(1,id2)
        dVdyC(2,id2) = dVdyFace(1,id2)

        dVdxC(3,id3) = dVdxFace(1,id3)
        dVdyC(3,id3) = dVdyFace(1,id3)

        
        ! Calculate weights
        
        g(1,:) = dVdxC(1,:)**2.0_dp + dVdyC(1,:)**2.0_dp
        g(2,:) = dVdxC(2,:)**2.0_dp + dVdyC(2,:)**2.0_dp
        g(3,:) = dVdxC(3,:)**2.0_dp + dVdyC(3,:)**2.0_dp

        w(1,:) = (g(2,:)*g(3,:)+epsilon)/(g(1,:)**2.0_dp+g(2,:)**2.0_dp+g(3,:)**2.0_dp + 3.0_dp*epsilon)
        w(2,:) = (g(1,:)*g(3,:)+epsilon)/(g(1,:)**2.0_dp+g(2,:)**2.0_dp+g(3,:)**2.0_dp + 3.0_dp*epsilon)
        w(3,:) = (g(2,:)*g(1,:)+epsilon)/(g(1,:)**2.0_dp+g(2,:)**2.0_dp+g(3,:)**2.0_dp + 3.0_dp*epsilon)

        ! Calculate limited gradient
        LdVdx(1,:) = w(1,:)*dVdxC(1,:) + w(2,:)*dVdxC(2,:) + w(3,:)*dVdxC(3,:)
        LdVdy(1,:) = w(1,:)*dVdyC(1,:) + w(2,:)*dVdyC(2,:) + w(3,:)*dVdyC(3,:)

        
        ! Calculate the limited correction to average cell value
        dV(:,:,i) = dx*matmul(ones,LdvdX) + dy*matmul(ones,LdVdy)

        ! Average vell value placed at all nodal values
        aV(:,:,i) = matmul(ones,VC0)
    end do
   
    ! Limited density
    rho = aV(:,:,1) + dV(:,:,1)
   
    
    ! Check for negative densities
    allocate(idRhoNeg(0),idPNeg(0))
    ! Find min value in each element
    minArrayDensity = minval(rho,dim=1)
    ! Find location of elements less than some tolerance
    do i = 1,K
        if (minArrayDensity(i).lt.0) then
            idRhoNeg = [idRhoNeg, i]
        end if  
    end do
    
    ! Correct by halving gradient in problem elements
    ! And then check for problem elements
    ! Repeat until no more problem elements
   
    do while(size(idRhoNeg).gt.0)
        
        dV(:,idRhoNeg,1) = 0.5_dp*dV(:,idRhoNeg,1)
        rho = aV(:,:,1) + dV(:,:,1)

        minArrayDensity = minval(rho,dim=1)

        deallocate(idRhoNeg)
        allocate(idRhoNeg(0))
        do i = 1,K
            if (minArrayDensity(i).lt.0) then
                idRhoNeg = [idRhoNeg, i]
            end if
         end do
        
         
    end do

    

    ! Limited momentum, rhou = averRhou + d(rho*u) = averRhou + drho*u + du*rho
    rhou = averRhou + dV(:,:,1)*aV(:,:,2) + dV(:,:,2)*aV(:,:,1)
    rhov = averRhov + dV(:,:,1)*aV(:,:,3) + dV(:,:,3)*aV(:,:,1)
    
    ener = averEner + 1.0_dp/(gammaGas -1.0_dp)*dV(:,:,4) + 0.5_dp*dV(:,:,1)*(aV(:,:,2)**2.0_dp + aV(:,:,3)**2.0_dp) + aV(:,:,1)*(aV(:,:,2)*dV(:,:,2) + aV(:,:,3)*dV(:,:,3))
    ! Update Limiter2D output argument
    Q(:,:,1) = rho
    Q(:,:,2) = rhou
    Q(:,:,3) = rhov
    Q(:,:,4) = ener

    
    deallocate(idRhoNeg,idPNeg)
    deallocate(BCTypeFlat)
    deallocate(idW, idMW, idI, idO, idS, idD, idN, idB)
    deallocate(H1,H2,H3)
    deallocate(id1,id2,id3)
    deallocate(unPackMask)
    deallocate(ones)
end subroutine Limiter2D

subroutine Print2DMatrix(M)
    ! Prints a 2D matrix
    ! Arguments:    M       -   Matrix to be printed
    !               shapeM  -   Shape of M (input as shape(M)) 

    ! Author: Sam MacPherson
    implicit none

    real(dp), dimension(:,:), intent(in)    :: M
    integer :: i,j

    do i = 1,size(M(:,1))
        print*, (M(i,j), j=1,size(M(1,:)))
    end do


end subroutine Print2DMatrix

subroutine Residual(Q2,Q1,dt,R)
    ! Subroutine to calculate the Residuals using the same process as fluent
    ! Arguments:    Q2          -   Latest vector of conserved variable
    !               Q1          -   Previous vector of conserved variables
    !               dt          -   Timestep
    !               timestep    -   The timestep
    !               R           -   Vector of Residuals

    ! Author: Sam Macpherson

    implicit none
    ! Define arguments
    real(dp), dimension(:,:,:), intent(in)  :: Q2, Q1
    real(dp), dimension(Np,K,4)             :: RMatrix
    real(dp), intent(in)                    :: dt
    real(dp), dimension(4), intent(out)     :: R

    integer :: i

    RMatrix = (Q2 - Q1)/dt
    RMatrix = RMatrix**2
    do i = 1,4
        R(i) = sqrt(sum(RMatrix(:,:,i))/Np/K)
    end do
    
end subroutine Residual

subroutine IdentityMatrix(N,eye)
    ! Computes identity matrix of shape NxN
    ! Arguments:    N   -   Size of matrix
    !               eye -   Matrix to be returned

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    integer, intent(in)                     ::  N
    real(dp), dimension(N,N), intent(out)   ::  eye

    integer ::  i

    eye = 0.0_dp

    do i = 1,N
        eye(i,i) = 1.0_dp
    end do
end subroutine IdentityMatrix



end Module