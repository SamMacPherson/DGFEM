Module Geometric
use GlobalVariables
use MatrixOperations

contains

subroutine Startup()
    ! Author: Sam MacPherson

    real(dp), dimension(size(Nodex))                        ::  r, s
    real(dp), dimension(size(r),(int(N)+2)*(int(N)+1)/2)    ::  Dr, Ds
    real(dp), dimension(Np,K)                               ::  rx, sx, ry, sy

    
    call xytors(Nodex,Nodey,r,s)
    
    allocate(J(Np,K),V(Np,Np))
    
    call Vandermonde2D(r,s,N,V)
    
    call Dmatrices2D(r,s,N,V,Dr,Ds)
    

    call GeometricFactors2D(x,y,Dr,Ds,rx,sx,ry,sy,J)
    
   
end subroutine Startup

subroutine xytors(x,y,r,s)
    ! Tranforms to r,s co-ordinates on straight sided triangle from x,y geometric co-ordinates
    ! Author: Sam MacPherson

    implicit none

    ! Define arguments
    real(dp), dimension(:), intent(in)  ::  x,y
    real(dp), dimension(:), intent(out) ::  r,s

    real(dp), dimension(size(x))        ::  L1, L2, L3

    L1 = (sqrt(3.0_dp)*y+1.0_dp)/3.0_dp
    L2 = (-3.0_dp*x - sqrt(3.0_dp)*y + 2.0_dp)/6.0_dp
    L3 = ( 3.0_dp*x - sqrt(3.0_dp)*y + 2.0_dp)/6.0_dp

    r = -L2 + L3 - L1
    s = -L2 - L3 + L1

end subroutine xytors

subroutine rstoab(r,s,a,b)
    ! Transforms from straight sided to equilateral triangle nodes
    ! Author: Sam MacPherson

    implicit none

    ! Define arguments
    real(dp), dimension(:), intent(in)    ::  r, s
    real(dp), dimension(:), intent(out)   ::  a, b

    integer ::  i

    a = 0
 
    do i = 1,Np
        if(s(i).ne.1) then
            a(i) = 2.0_dp*(1.0_dp + r(i))/(1.0_dp-s(i)) - 1.0_dp
        else
            a(i) = -1
        end if
    end do

    b = s

end subroutine rstoab

subroutine JacobiP(x,alpha,beta,NJac,P)
    ! Author: Sam MacPherson

    implicit none

    ! Define arguments
    real(dp), dimension(:), intent(in)  ::  x
    real(dp), intent(in)                ::  alpha, beta, NJac
    real(dp), dimension(:), intent(out) ::  P

    real(dp), dimension(Njac+1,size(x)) ::  PL
    real(dp)    ::  gamma0, gamma1, aold, anew, bnew, h1
    intrinsic gamma

    integer ::  i
    

    PL = 0

    gamma0 = 2.0_dp**(alpha+beta+1.0_dp)/(alpha+beta+1.0_dp)*gamma(alpha+1.0_dp)*gamma(beta+1.0_dp)/gamma(alpha+beta+1.0_dp)

    PL(1,:) = 1.0_dp/sqrt(gamma0)

    if (Njac.eq.0) then
        P = PL(1,:)
        return
    end if

    gamma1 = (alpha+1.0_dp)*(beta+1.0_dp)/(alpha+beta+3.0_dp)*gamma0
    PL(2,:) = ((alpha+beta+2.0_dp)*x/2.0_dp + (alpha-beta)/2.0_dp )/sqrt(gamma1)

    if (Njac.eq.1) then
        P = PL(2,:)
        return
    end if

    aold = 2.0_dp/(2.0_dp+alpha+beta)*sqrt((alpha+1.0_dp)*(beta+1.0_dp)/(alpha+beta+3.0_dp))

    do i = 1,Njac-1
        h1 = 2.0_dp*real(i)+alpha+beta
        anew = 2.0_dp/(h1+2.0_dp)*sqrt( (real(i)+1.0_dp)*(real(i)+1.0_dp+alpha+beta)*(real(i)+1.0_dp+alpha)*(real(i)+1.0_dp+beta)/(h1+1.0_dp)/(h1+3.0_dp))
        bnew = - (alpha**2.0_dp-beta**2.0_dp)/h1/(h1+2.0_dp)

        PL(i+2,:) = 1.0_dp/anew*( -aold*PL(i,:) + (x-bnew)*PL(i+1,:))
        aold =anew
    end do
    P = PL(Njac+1,:)

end subroutine JacobiP

subroutine Simplex2DP(a,b,i,j,P)
    ! Author: Sam MacPherson

    implicit none

    ! Define arguments
    real(dp), intent(in)                ::  i, j
    real(dp), dimension(:), intent(in)  ::  a, b
    real(dp), dimension(:), intent(out) ::  P
    real(dp), dimension(size(P))        ::  h1, h2

    call JacobiP(a,0.0_dp,0.0_dp,i,h1)
    call JacobiP(b, 2.0_dp*i + 1.0_dp, 0.0_dp, j, h2)

    P = sqrt(2.0_dp)*h1*h2*(1-b)**i
end subroutine Simplex2DP

subroutine Vandermonde2D(r,s,N,V2D)
    ! Calculates the Vandermonde matrix
    ! Author: Sam MacPherson

    implicit none

    ! Define arguments
    real(dp), dimension(:), intent(in)      ::  r, s
    real(dp), intent(in)                    ::  N
    real(dp), dimension(:,:), intent(out)   ::  V2D
    real(dp), dimension(size(r))            ::  a, b

    integer ::  sk, i, j

    V2D = 0.0_dp

    call rstoab(r,s,a,b)
   
    sk = 1
    do i = 0,int(N)
        do j = 0,int(N) - i
            
            call Simplex2DP(a,b,real(i),real(j),V2D(:,sk))
            sk = sk + 1
        end do
    end do

end subroutine Vandermonde2D

subroutine GradJacobiP(r,alpha,beta,NJac,deltaP)
    ! Author: Sam MacPherson

    implicit none

    ! Define arguments
    real(dp), dimension(:), intent(in)  ::  r
    real(dp), intent(in)                ::  alpha, beta, NJac
    real(dp), dimension(size(r)), intent(out) ::  deltaP

    deltaP = 0.0_dp

    if (int(N).eq.0) then
        deltaP = 0.0_dp
    else
        call JacobiP(r, alpha + 1.0_dp, beta + 1.0_dp,Njac-1.0_dp,deltaP)
        deltaP = deltaP*sqrt(NJac*(NJac+alpha+beta+1.0_dp))
    end if
    

    
end subroutine 

subroutine GradSimplex2DP(a,b,id,jd,dmodedr,dmodeds)
    ! Author: Sam MacPherson
    implicit none

    ! Define arguments
    real(dp), dimension(:), intent(in)          ::  a, b
    real(dp), intent(in)                        ::  id, jd
    real(dp), dimension(size(a)), intent(out)   ::  dmodedr, dmodeds
    real(dp), dimension(size(a))                ::  fa, dfa, gb, dgb, tmp

    call JacobiP(a,0.0_dp,0.0_dp,id,fa)
    call GradJacobiP(a,0.0_dp,0.0_dp,id,dfa)

    call JacobiP(b,2.0_dp*id+1.0_dp,0.0_dp,jd,gb)
    call GradJacobiP(b,2.0_dp*id+1.0_dp,0.0_dp,jd,dgb)

    dmodedr = dfa*gb
    if (id.gt.0) then
        dmodedr = dmodedr * ((0.5_dp*(1.0_dp-b))**(id-1.0_dp))
    end if

    dmodeds = dfa*(gb*(0.5_dp*(1+a)))
    if (id.gt.0) then
        dmodeds = dmodeds * ((0.5_dp*(1.0_dp-b))**(id-1.0_dp))
    end if

    tmp = dgb*((0.5_dp*(1.0_dp-b))**id)

    if(id.gt.0) then
        tmp = tmp-0.5_dp*id*gb*((0.5_dp*(1.0_dp-b))**(id-1.0_dp))
    end if

    dmodeds = dmodeds+fa*tmp

    dmodedr = 2**(id+0.5_dp)*dmodedr 
    dmodeds = 2**(id+0.5_dp)*dmodeds

end subroutine GradSimplex2DP

subroutine GradVandermonde2D(r,s,N,V2Dr,V2Ds)
    ! Author: Sam MacPherson

    ! Define arguments
    real(dp), dimension(:), intent(in)      ::  r, s
    real(dp), intent(in)                    ::  N
    real(dp), dimension(:,:), intent(out)   ::  V2Dr, V2Ds

    real(dp), dimension(size(r))            ::  a, b
    integer ::  sk, i, j

    V2Dr = 0.0_dp
    V2Ds = 0.0_dp

    sk = 1
    do i = 0,int(N)
        do j = 0,int(N)-i
            call GradSimplex2DP(a,b,real(i),real(j),V2Dr(:,sk),V2Ds(:,sk))
            sk = sk + 1
        end do
    end do
end subroutine GradVandermonde2D

subroutine Dmatrices2D(r,s,N,V2D,Dr,Ds)
    ! Author: Sam MacPherson

    ! Define arguments
    real(dp), dimension(:), intent(in)        ::  r, s
    real(dp), intent(in)                      ::  N
    real(dp), dimension(:,:), intent(in)      ::  V2D
    real(dp), dimension(Np,Np)                ::  invV2D
    real(dp), dimension(Np,Np), intent(out)   ::  Dr, Ds

    real(dp), dimension(Np,Np)                ::  Vr, Vs
  
    invV2D = inv(V2D)
   
    call GradVandermonde2D(r,s,N,Vr,Vs)

    
    Dr = matmul(Vr,invV2d)
    Ds = matmul(Vs,invV2d)
    


end subroutine Dmatrices2D

subroutine GeometricFactors2D(x,y,Dr,Ds,rx,sx,ry,sy,J)
    ! Author Sam MacPherson

    ! Define arguments
    real(dp), dimension(:,:), intent(in)                        ::  x, y
    real(dp), dimension(size(x(:,1)),size(x(:,1))), intent(in)  ::  Dr, Ds
    real(dp), dimension(size(x(:,1)),size(x(1,:))), intent(out) ::  rx, sx, ry, sy, J  

    real(dp), dimension(size(x(:,1)),size(x(1,:)))   ::  xr, xs, yr, ys

    xr = matmul(Dr,x)
    xs = matmul(Ds,x)

    yr = matmul(Dr,y)
    ys = matmul(Ds,y)

    J = -xs*yr + xr*ys

    rx =  ys/J
    sx = -yr/J
    ry = -xs/J
    sy =  xr/J
    


end subroutine GeometricFactors2D



end Module