Module PrimConTransformation
use GlobalVariables
contains

subroutine ConToPrim(Q,rho,u,v,p)
    ! Extracts the primitive variables from the vector of conserved variables
    ! Arguments:    Q       -   Vector of conserved variables
    !               rho     -   Matrix of densities
    !               u       -   Matrix of u velocities
    !               v       -   Matrix of v velocities
    !               p       -   Matrix of pressures

    ! Author: Sam MacPherson

    implicit none

    ! Define arguments

    real(dp), dimension(:,:,:), intent(in)  :: Q
    real(dp), dimension(:,:),   intent(out) :: rho,u,v,p

    
    ! Field variables
    real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1))) :: rhou,rhov,ener

    rho  = Q(:,:,1)
    rhou = Q(:,:,2)
    rhov = Q(:,:,3)
    ener = Q(:,:,4)
   
    u = rhou/rho
    v = rhov/rho
    p = (gammaGas-1.0_dp)*(ener-0.5_dp*rho*(u**2.0_dp+v**2.0_dp))
end subroutine ConToPrim

subroutine PrimToCon(Q,rho,u,v,p)
    ! constructs the conservative field vector from the primitive variables
    ! Arguments:    Q       -   Vector of conserved variables
    !               rho     -   Matrix of densities
    !               u       -   Matrix of u velocities
    !               v       -   Matrix of v velocities
    !               p       -   Matrix of pressures

    ! Author: Sam MacPherson

    
    implicit none

    ! Define arguments
    real(dp), dimension(:,:,:), intent(out) :: Q
    real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1))), intent(in)            :: rho,u,v,p

    Q(:,:,1) = rho
    Q(:,:,2) = rho*u
    Q(:,:,3) = rho*v
    Q(:,:,4) = p/(gammaGas-1.0_dp) +0.5_dp*rho*(u**2.0_dp+v**2.0_dp)


end subroutine PrimToCon

subroutine GetTemperature(p,rho,T)
    ! Calculates the temperature using ideal gas law
    ! Arguments:    p   -   Matrix of pressure values
    !               rho -   Matrix of density values
    !               T   -   Matrix of temperature values

    ! Author: Sam MacPherson

    implicit none

    real(dp), dimension(:,:), intent(in)                            :: p, rho
    real(dp), dimension(size(p(:,1)),size(p(1,:)) ), intent(out)    :: T    

    T = p/(rho*GasConstant)

end subroutine GetTemperature
    
subroutine Sutherland(T,mu)
    ! Calculates the dynamic viscosity according to Sutherlands law
    ! Arguments:    T   -   Matrix of temperatures
    !               mu  -   Dynamic viscosity

    ! Author: Sam MacPherson

    implicit none
    
    ! Define arguments
    real(dp), dimension(:,:), intent(in)            :: T
    real(dp), dimension(size(T(:,1)),size(T(1,:)))  :: mu

    mu = mu0*(T/T0)**(1.5_dp)*(T0 + S)/(T+S)
end subroutine Sutherland

subroutine SoundSpeed(p,rho,c)
    ! Calculates the sound speed at every location in the grid
    ! Arguments:    p       -   Matrix of pressure values
    !               rho     -   Matrix of density values
    !               c       -   Speed of sound

    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:,:), intent(in)                        :: p, rho
    real(dp), dimension(size(p(:,1)),size(p(1,:))), intent(out) :: c

    c = sqrt(abs(gammaGas*p/rho))

end subroutine SoundSpeed



end Module