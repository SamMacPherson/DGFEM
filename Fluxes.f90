Module Fluxes
use GlobalVariables
use MatrixOperations
use PrimConTransformation
contains


subroutine LaxFriedrichFlux(FM,FP,gQM,gQP,flux)
    ! Function to calculate the Lax-Friedrich flux
    ! Arguments:    FM      - Interior flux (Left state) at Gauss nodes
    !               FP      - Exterior flux (Right state) at Gauss nodes
    !               gQM     - Conserved vector at interior (left) Gauss nodes
    !               gWP     - Conserved vector at exterior (right) Gauss nodes
    !               flux    - Flux to be returned

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(in)              ::  FM, FP
    real(dp), dimension(:,:), allocatable               ::  grho, gu, gv, gp
    real(dp), dimension(:,:,:), intent(inout)              ::  gQM, gQP
    real(dp), dimension(:,:,:), intent(out)             ::  flux
   

    real(dp), dimension(:,:), allocatable               ::  fluxTemp



    ! Define useful lambda matrices
    real(dp), dimension(:,:), allocatable   ::  glambda
    real(dp), dimension(NGauss,3*K)         ::  glambdafaces
    real(dp), dimension(:), allocatable     ::  glambdaflat    

    ! Masks for unpacking variables
    logical, dimension(3*NGauss,K)   ::  unPackMask
    logical, dimension(NGauss,3*K)   ::  unPackMaskglambdaFaces


    ! Loop counters
    integer :: i,j

    real:: start, finish


    grho = gQM(:,:,1)
    gu   = gQM(:,:,2)/grho
    gv   = gQM(:,:,3)/grho
    gp   = (gammaGas-1.0_dp)*(gQM(:,:,4)-0.5_dp*grho*(gu**2.0_dp+gv**2.0_dp))
    
    unPackMask = .true.
    unPackMaskglambdaFaces = .true.
    

    ! Calculate lambda for Lax-Friedrich scheme
   
    glambda = sqrt(gu**2.0_dp + gv**2.0_dp) + sqrt(abs(gammaGas*gp/grho))
    glambdaflat = pack(glambda,.true.)
    
    
    ! Loop through and find maximum between the interior and exterior Gauss nodes.
    do i = 1,size(glambdaflat)
        glambdaflat(i) = maxval( [glambdaflat(gaussmapM(i)),glambdaflat(gaussmapP(i))]    )
    end do

    
    ! Reform the matrix to Ngauss x 3*K shape, ie each column represents a face
    glambdafaces = unpack(glambdaflat,unPackMaskglambdaFaces,glambdafaces)

    
    ! Loop through each face, and set the value at each Gauss node to the max value on its face
    do i = 1,size(glambdafaces(1,:))
        glambdafaces(:,i) = maxval(glambdafaces(:,i))
    end do
    
    ! Reform glambda back to NGauss*3 x K shape, ie each column represents an element
    ! Pack and unpack only work between vectors and arrays, so need to flatten before reforming again
    glambdaflat = pack(glambdafaces,.true.)
  
    glambda     = unpack(glambdaflat,unPackMask,glambda)

    do i = 1,4
        ! Lax-Friedrich flux
        flux(:,:,i) = (FM(:,:,i)+FP(:,:,i))/2.0_dp + glambda*(gQM(:,:,i) - gQP(:,:,i))/2.0_dp

        ! Convert back to polynomial reconstruction nodes
        !flux(:,:,i) = matmul(transpose(gaussinterp),gaussW*fluxTemp)
    end do
 
  
    
    deallocate(glambda, glambdaflat)
   
end subroutine LaxFriedrichFlux

subroutine HLLFlux(FM,FP,gQM,gQP,flux)
    ! Function to calculate HLL flux
    ! Arguments:    FM      - Interior flux (Left state) at Gauss nodes
    !               FP      - Exterior flux (Right state) at Gauss nodes
    !               gQM     - Conserved vector at interior (left) Gauss nodes
    !               gQP     - Conserved vector at exterior (right) Gauss nodes
    !               flux    - Flux to be returned

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(in)      ::  FM, FP
    real(dp), dimension(:,:,:), intent(inout)   ::  gQM, gQP
    real(dp), dimension(:,:,:), intent(out)     ::  flux


    real(dp), dimension(3*NGauss,K,4)   ::  fluxRotated
    real(dp), dimension(3*NGauss,K,4)   ::  fluxTemp
    real(dp), dimension(3*NGauss,K)     ::  grhouM, grhouP, grhovM, grhovP
    real(dp), dimension(3*NGauss,K)     ::  H, HM, HP, HTilde, rhoM, rhoP, uM, uP, vM
    real(dp), dimension(3*NGauss,K)     ::  vP, pM, pP, utilde, atilde, enerM, enerP, cM, cP
    real(dp), dimension(3*NGauss,K)     ::  ubar, dbar
    real(dp), dimension(3*NGauss,K)     ::  SM, SP, eta
    real(dp), dimension(3*NGauss,K)     ::  minSP0, min0SM, t1, t2, t3

    integer ::  i,j,n

 

    grhouM = gQM(:,:,2)
    grhouP = gQP(:,:,2)

    grhovM = gQM(:,:,3)
    grhovP = gQP(:,:,3)

    

    ! Transform to coordinate system normal and tangent to the element boundary
    gQM(:,:,2) =  gaussnx*grhouM  + gaussny*grhovM
    gQM(:,:,3) = -gaussny*grhouM  + gaussnx*grhovM

    gQP(:,:,2) =  gaussnx*grhouP  + gaussny*grhovP
    gQP(:,:,3) = -gaussny*grhouP  + gaussnx*grhovP

    enerM = gQM(:,:,4)
    enerP = gQP(:,:,4)

    ! Convert to primitive variables
    call ConToPrim(gQM,rhoM,uM,vM,pM)
    call ConToPrim(gQP,rhoP,uP,vP,pP)

 

    call SoundSpeed(pM,rhoM,cM)
    call SoundSpeed(pP,rhoP,cP)
    HM = (enerM + pM)/rhoM
    HP = (enerP + pP)/rhoP


    ! Roe averaged enthalpy
    HTilde = ( sqrt(rhoM)*HM + sqrt(rhoP)*HP ) / (sqrt(rhoM) + sqrt(rhoP))

    ! Roe averaged velocity
    utilde = ( sqrt(rhoM)*uM + sqrt(rhoP)*uP ) / (sqrt(rhoM) + sqrt(rhoP))

    
   
    ! Sound speed estimate
    atilde = sqrt( (gammaGas-1.0_dp)*(HTilde-0.5_dp*utilde**2.0_dp) )

    eta = 0.5_dp*( sqrt(rhoM)*sqrt(rhoP) )/ ( sqrt(rhoM) + sqrt(rhoP) )**2.0_dp

    dbar = ( sqrt(rhoM)*cM**2.0_dp + sqrt(rhoP)*cP**2.0_dp ) / (sqrt(rhoM) + sqrt(rhoP)) &
         +   eta*(uP-uM)**2.0_dp
    dbar = sqrt(dbar)
    
    ubar = 0.5_dp*(uP+uM)

    ! Wave speed estimate
    SM = utilde - atilde
    SP = utilde + atilde
    
    ! Apply HLL flux in a matrix like way
    do i = 1,3*NGauss
        do j = 1,K

            minSP0(i,j) = minval([SP(i,j), 0.0_dp ])
            min0SM(i,J) = minval([0.0_dp, SM(i,j)])
        end do
    end do

    ! Apply HLL flux in a matrix like way
    t1 = (minSP0 - min0SM)/(SP - SM)
    t2 = 1.0_dp - t1
    t3 = (SP*abs(SM) - SM*abs(SP))/(2*(SP-SM))
 
    ! Apply HLL flux in a matrix like way
    gQM(:,:,2) = grhouM
    gQM(:,:,3) = grhovM

    gQP(:,:,2) = grhouP
    gQP(:,:,3) = grhovP

    do i = 1,4
        flux(:,:,i) = t1*FP(:,:,i) + t2*FM(:,:,i) + t3*(gQM(:,:,i)-gQP(:,:,i))
        
    end do
 
   

  
end subroutine HLLFlux

subroutine HLLCFlux(FM,FP,gQM,gQP,flux)
    ! Function to calculate HLL flux
    ! Arguments:    FM      - Interior flux (Left state) at Gauss nodes
    !               FP      - Exterior flux (Right state) at Gauss nodes
    !               gQM     - Conserved vector at interior (left) Gauss nodes
    !               gQP     - Conserved vector at exterior (right) Gauss nodes
    !               flux    - Flux to be returned

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(in)      ::  FM, FP
    real(dp), dimension(:,:,:), intent(inout)   ::  gQM, gQP
    real(dp), dimension(:,:,:), intent(out)     ::  flux


    real(dp), dimension(3*NGauss,K,4)   ::  fluxRotated
    real(dp), dimension(3*NGauss,K,4)   ::  fluxTemp, UstarM, UstarP, UstarMRot, UstarPRot
    real(dp), dimension(3*NGauss,K)     ::  grhouM, grhouP, grhovM, grhovP
    real(dp), dimension(3*NGauss,K)     ::  H, HM, HP, HTilde, rhoM, rhoP, uM, uP, vM
    real(dp), dimension(3*NGauss,K)     ::  vP, pM, pP, utilde, atilde, enerM, enerP, cM, cP
    real(dp), dimension(3*NGauss,K)     ::  SM, SP, SStar
    


    integer ::  i,j,n



    grhouM = gQM(:,:,2)
    grhouP = gQP(:,:,2)

    grhovM = gQM(:,:,3)
    grhovP = gQP(:,:,3)

    

    ! Transform to coordinate system normal and tangent to the element boundary
    gQM(:,:,2) =  gaussnx*grhouM  + gaussny*grhovM
    gQM(:,:,3) = -gaussny*grhouM  + gaussnx*grhovM

    gQP(:,:,2) =  gaussnx*grhouP  + gaussny*grhovP
    gQP(:,:,3) = -gaussny*grhouP  + gaussnx*grhovP

    enerM = gQM(:,:,4)
    enerP = gQP(:,:,4)

    ! Convert to primitive variables
    call ConToPrim(gQM,rhoM,uM,vM,pM)
    call ConToPrim(gQP,rhoP,uP,vP,pP)

 

    call SoundSpeed(pM,rhoM,cM)
    call SoundSpeed(pP,rhoP,cP)
    HM = (enerM + pM)/rhoM
    HP = (enerP + pP)/rhoP


    ! Roe averaged enthalpy
    HTilde = ( sqrt(rhoM)*HM + sqrt(rhoP)*HP ) / (sqrt(rhoM) + sqrt(rhoP))

    ! Roe averaged velocity
    utilde = ( sqrt(rhoM)*uM + sqrt(rhoP)*uP ) / (sqrt(rhoM) + sqrt(rhoP))

    
   
    ! Sound speed estimate
    atilde = sqrt( (gammaGas-1.0_dp)*(HTilde-0.5_dp*utilde**2.0_dp) )

    ! Wave speed estimate
    SM = utilde - atilde
    SP = utilde + atilde
    SStar = ( pP - pM + rhoM*uM*(SM - uM) - rhoP*uP*(SP - uP)  ) / ( rhoM*(SM-uM) - rhoP*(SP - uP)  )



    UstarM(:,:,1)       = rhoM*(SM-uM)/(SM-SStar)
    UstarMRot(:,:,2)    = rhoM*(SM-uM)/(SM-SStar)*SStar
    UstarMRot(:,:,3)    = rhoM*(SM-uM)/(SM-SStar)*vM
    UstarM(:,:,4)       = rhoM*(SM-uM)/(SM-SStar)*( enerM/rhoM + (SStar - uM)*(SStar + pM/(rhoM*(SM-uM))  )  )

    UstarP(:,:,1)       = rhoP*(SP-uP)/(SP-SStar)
    UstarPRot(:,:,2)    = rhoP*(SP-uP)/(SP-SStar)*SStar
    UstarPRot(:,:,3)    = rhoP*(SP-uP)/(SP-SStar)*vP
    UstarP(:,:,4)       = rhoP*(SP-uP)/(SP-SStar)*( enerP/rhoP + (SStar - uP)*(SStar + pP/(rhoP*(SP-uP))  )  )


    UstarM(:,:,2) = gaussnx*UstarMRot(:,:,2) - gaussny*UstarMRot(:,:,3)
    UstarM(:,:,3) = gaussny*UstarMRot(:,:,2) + gaussnx*UstarMRot(:,:,3)

    UstarP(:,:,2) = gaussnx*UstarPRot(:,:,2) - gaussny*UstarPRot(:,:,3)
    UstarP(:,:,3) = gaussny*UstarPRot(:,:,2) + gaussnx*UstarPRot(:,:,3)
 
    gQM(:,:,2) = grhouM
    gQM(:,:,3) = grhovM

    gQP(:,:,2) = grhouP
    gQP(:,:,3) = grhovP

    do n = 1,4
        do i = 1,3*NGauss
            do j = 1,K

                if      (0.0_dp.le.SM(i,j))                 then

                    flux(i,j,n) = FM(i,j,n)

                elseif  (SM(i,j).le.0.0_dp.le.SStar(i,j))   then

                    flux(i,j,n) = FM(i,j,n) + SM(i,j)*(UstarM(i,j,n) - gQM(i,j,n))

                elseif  (SStar(i,j).le.0.0_dp.le.SP(i,j))   then

                    flux(i,j,n) = FP(i,j,n) + SP(i,j)*(UstarP(i,j,n) - gQP(i,j,n))

                elseif  (0.0_dp.ge.SP(i,j))                 then
                    
                    flux(i,j,n) = FP(i,j,n)

                end if
            end do
        end do
    end do
    
   

  
end subroutine HLLCFlux

subroutine CenteredFlux(FM,FP,flux)
    ! Function to calculate the flux for viscous terms
    ! Arguments:    FM      - Interior flux (Left state) at Gauss nodes
    !               FP      - Exterior flux (Right state) at Gauss nodes
    !               flux    - Flux to be returned

    ! Author: Sam MacPherson
    implicit none

    real(dp), dimension(:,:,:), intent(in)              ::  FM, FP
    real(dp), dimension(:,:,:), intent(out)             ::  flux


    flux = 0.5_dp*(FM+FP)
end subroutine CenteredFlux









end Module