Program Main
! Script that solves the 2D CNS equations using a DG-FEM
! Method uses the Hesthaven and Warbuton DG-FEM code as a reference
! The pre-processing Matlab script needs to first be run

! Author: Sam MacPherson


! External modules
use Load
use CNSSubroutines
use GlobalVariables
use BCandIC
use PrimConTransformation
use Geometric
implicit none

! Define conserved vectors
real(dp), dimension(:,:,:), allocatable :: Q1,Q2,Q3, rhsQ, ki, pi
integer, dimension(3)                   :: shapeQ
! Primitive variables
real(dp), dimension(:,:), allocatable   :: rho, uvel, vvel, p
! Residual
real(dp), dimension(4)  :: R
! Define time parameters
real(dp)    :: currentTime
integer     :: timestep
logical     :: lastTimeStep, autoTimeStep, residualCheck, residualReached

integer :: i,a,b,c

real    :: start, finish

! 4th order LSERK Coefficients
aL = [0., -567301805773./1357537059087., -2404267990393./2016746695238., -3550918686646./2091501179385., -1275806237668./842570457699. ]
bL = [1432997174477./9575080441755., 5161836677717./13612068292357., 1720146321549./2090206949498., 3134564353537./4481467310338., 2277821191437./14882151754819.]
! cL = [0, 1432997174477/9575080441755, 2526269341429/6820363962896, 2006345519317/3224310063776, 2802321613138/2924317926251]

! aL = [0., 0.481231743, -1.049562606, -1.602529574, -1.778267193]
! bL = [0.097618354, 0.412253292, 0.440216964, 1.426311463, 0.197876053]

lastTimeStep = .false.
autoTimeStep = .false.
residualCheck = .false.
residualReached = .false.


! Read all matrices from the Matlab preprocessing script
call ReadAllMatrices()
print*, "Read matrices into Fortran"
print*, " "
Ncub = size(cubV(:,1))
! Read in .ini file

call ReadIniFile()

call Startup()

! Set up Q vectors, use x just for shape
shapeQ(1) = size(x(:,1))
shapeQ(2) = size(x(1,:))
shapeQ(3) = 4
allocate( Q1( shapeQ(1),shapeQ(2),shapeQ(3)),Q2(shapeQ(1),shapeQ(2),shapeQ(3)))
allocate(rho(shapeQ(1),shapeQ(2)),uvel(shapeQ(1),shapeQ(2)),vvel(shapeQ(1),shapeQ(2)),p(shapeQ(1),shapeQ(2)))

Q1(:,:,1) = x
Q1(:,:,2) = x
Q1(:,:,3) = x
Q1(:,:,4) = x
Q2 = Q1
Q3 = Q1
rhsQ = Q1


! Initial conditions
select case(testProblem)
    case(0)
        call ZeitounTestCaseIC(Q1)
    case(1)
        call ShearFlow2DIC(Q1)
    case default
        print*, "Warning: Invalid test problem selected"
end select

ki = Q1
pi = Q1

currentTime = 0
timestep = 0

! Output initial conditions
call ConvertToParaView(Q1,timestep,.false.)


if(dt.eq.0) then
    ! Compute initial timestep only if not specified in setup file
    call CNSdt2D(Q1,dt)
    dt = CFL*dt
    ! Set flag so code knows it needs to determine time step
    autoTimeStep = .true.
end if

if (totalTime.eq.0) then
    residualCheck = .true.
end if



mainLoop: do 
    ! Change dt if its the last timestep and residuals aren't being used to stop the simulation
    if (residualCheck.eq..false.) then
    if ( currentTime + dt > totalTime ) then
        dt = totalTime - currentTime   
        lastTimeStep = .true.
    end if
    end if

    select case (RKOrder)
        case(2)
            
            !print*,"dt",dt
            !print*,Q1(1,2463,4)
            ! 2nd order 2 stage SSP
            call CNS2DRHS(Q1,rhsQ)
            
            !print*,rhsQ(1,2463,4)
            Q2 = Q1 + dt*rhsQ
             
            
            !print*,Q2(1,2463,:)
            ! Apply limitier
            if(useLimiter.eq.1) then
                call Limiter2D(Q2)
            end if
            !print*,Q2(1,2463,:)
            
            
            call CNS2DRHS(Q2,rhsQ)
            !print*,rhsQ(1,2463,:)
            Q2 = 0.5*(Q1 + Q2 + dt*rhsQ)
            !print*,Q2(1,2463,:)
            ! Apply limitier
            if(useLimiter.eq.1) then
                call Limiter2D(Q2)
            end if
            !print*,Q2(1,2463,4)
            
        case(3)
            call cpu_time(start)
            ! 3rd Order 3 stage SSP
            
            call CNS2DRHS(Q1,rhsQ)
            Q2 = Q1 + dt*rhsQ
            
            ! Apply limiter
            if(useLimiter.eq.1) then
                call Limiter2D(Q2)
            end if
            
           
            call CNS2DRHS(Q2,rhsQ)
            Q3 = 0.25_dp*(3.0_dp*Q1 + Q2 + dt*rhsQ)
            
            ! Apply limiter
            if(useLimiter.eq.1) then
                call Limiter2D(Q3)
            end if
           
            call CNS2DRHS(Q3,rhsQ)
            Q2 = 1.0_dp/3.0_dp*(Q1 + 2.0_dp*Q3 + 2.0_dp*dt*rhsQ)
            
            ! Apply limiter
            if(useLimiter.eq.1) then
                call Limiter2D(Q2)
            end if
            call cpu_time(finish)
        case(4)
            ! 4th Order LSERK
            ! Initialise arrays to use, for k just need to initialise something with correct dimensions
           
           
            do i = 1,5
              
                call CNS2DRHS(pi,rhsQ)
               
                ki = aL(i)*ki + dt*rhsQ
                
                pi = pi + bL(i)*ki
                
                ! Apply limitier
                if(useLimiter.eq.1) then
                    call Limiter2D(pi)
                end if
              
               
                
            end do

            
                
            Q2 = pi

        case default
            print*, "No RK selected, choose one of the implemented schemes"
            exit mainLoop
    end select





    ! Update time level
    currentTime = currentTime + dt
    timestep = timestep + 1

    call Residual(Q2,Q1,dt,R)

    ! Give update on simulation to terminal,
    ! output depends on whether residuals are being used to stop simulation
    if (Modulo(timestep,outputFrequency).eq.0) then

        ! Save data to file
        call ConvertToParaView(Q2,timestep,.false.)

        if (residualCheck.eq..true.) then
            print*," "
            print*," "
            print*,"Density residual = ",       R(1)
            print*,"X-momentum residual = ",    R(2)
            print*,"Y-momentum residual = ",    R(3)
            print*,"Energy residual = ",        R(4)

        else
            print*," "
            print*," "
            print*,currentTime/totalTime*100,"% of simulation has run"
            print*,"Timestep time: ", finish - start
        end if
    end if


    ! Decide whether simulation has finished or not
    ! Depends if residuals or time is the determining factor
    if (residualCheck.eq..true.) then

        ! Check if residuals have been met
        residualReached = .true.
        do i = 1,4
            if(R(i).gt.residualStop) then
                residualReached = .false.
            end if
        end do
        
        ! If residual is reached then save data and exit simulation
        if (residualReached) then
            call ConvertToParaView(Q2,timestep,.false.)
            call ConvertToParaView(Q2,timestep,.true.)

            print*," "
            print*," "
            print*,"Density residual = ",       R(1)
            print*,"X-momentum residual = ",    R(2)
            print*,"Y-momentum residual = ",    R(3)
            print*,"Energy residual = ",        R(4)

            EXIT
        end if

    else
        ! This section runs if time is used to deduce the end of simulation
        if (currentTime.eq.totalTime) then
        !if (timestep.eq.1) then
            print*," "
            print*," "
            print*, "Final time reached"
            call ConvertToParaView(Q2,timestep,.false.)
            call ConvertToParaView(Q2,timestep,.true.)
            EXIT
        end if

    end if
    

    ! Calculate new timestep only if not specified in setup file
    if (autoTimeStep) then
        call CNSdt2D(Q2,dt)
        dt = CFL*dt
    end if


    
    Q1 = Q2

end do mainLoop

end Program Main