Module Load
use GlobalVariables
use PrimConTransformation
implicit none

Contains

    subroutine ReadMatrixFromText(M,filename)
        ! Subroutine which reads in a matrix that has been saved in text form
        ! Arguments:    M           -   Matrix for data to be written in to, size nrows x ncols
        !               filename    -   Name of the file (with file extension) 

        ! Author: Sam MacPherson
        implicit none

        ! Counters
        integer :: i, j
        
        real(dp), dimension(:,:), allocatable, intent(out)  :: M
        integer                                             :: nrows, ncols
        character (*), intent(in)                           :: filename


        open(12,file = filename)
        
            ! Read in size of matrix
            read(12,*) nrows, ncols

            ! Allocate matrix
            Allocate( M(nrows,ncols) )

            ! Read in matrix
            read(12,*) ((M(i,j), j=1,ncols), i=1,nrows)

        close(12)
   
    end subroutine ReadMatrixFromText

    subroutine ReadIntegerMatrixFromText(M,filename)
        ! Subroutine which reads in a matrix that has been saved in text form
        ! Arguments:    M           -   Matrix for data to be written in to, size nrows x ncols
        !               filename    -   Name of the file (with file extension) 

        ! Author: Sam MacPherson
        implicit none

        ! Counters
        integer :: i, j
        
        integer, dimension(:,:), allocatable, intent(out)  :: M
        integer                                             :: nrows, ncols
        character (*), intent(in)                           :: filename


        open(12,file = filename)
        
            ! Read in size of matrix
            read(12,*) nrows, ncols

            ! Allocate matrix
            Allocate( M(nrows,ncols) )

            ! Read in matrix
            read(12,*) ((M(i,j), j=1,ncols), i=1,nrows)

        close(12)
   
    end subroutine ReadIntegerMatrixFromText

    subroutine ReadAllMatrices()
        ! Reads all matrices

        ! Author: Sam MacPherson

        ! Read in matrices from preprocessing

        call ReadMatrixFromText(x,"Matrices/x.dat")
        call ReadMatrixFromText(y,"Matrices/y.dat")

        call ReadMatrixFromText(NodexMatrix,"Matrices/Nodex.dat")
        call ReadMatrixFromText(NodeyMatrix,"Matrices/Nodey.dat")

        call ReadMatrixFromText(cubDr,"Matrices/cubDr.dat")
        call ReadMatrixFromText(cubDs,"Matrices/cubDs.dat")

        call ReadMatrixFromText(cubV,"Matrices/cubV.dat")
        call ReadMatrixFromText(cubW,"Matrices/cubW.dat")
        call ReadMatrixFromText(Fscale,"Matrices/Fscale.dat")
        call ReadIntegerMatrixFromText(vmapMMatrix,"Matrices/vmapM.dat")
        call ReadIntegerMatrixFromText(vmapPMatrix,"Matrices/vmapP.dat")

        call ReadMatrixFromText(cubrx,"Matrices/cubrx.dat")
        call ReadMatrixFromText(cubry,"Matrices/cubry.dat")
        call ReadMatrixFromText(cubsx,"Matrices/cubsx.dat")
        call ReadMatrixFromText(cubsy,"Matrices/cubsy.dat")

        call ReadMatrixFromText(gaussinterp,"Matrices/gaussinterp.dat")
        call ReadMatrixFromText(gaussW,"Matrices/gaussW.dat")

        
        call ReadIntegerMatrixFromText(gaussmapMMatrix,"Matrices/gaussmapM.dat")
        call ReadIntegerMatrixFromText(gaussmapPMatrix,"Matrices/gaussmapP.dat")

        call ReadMatrixFromText(gaussnx,"Matrices/gaussnx.dat")
        call ReadMatrixFromText(gaussny,"Matrices/gaussny.dat")

        !call ReadMatrixFromText(V,"Matrices/V.dat")
        call ReadMatrixFromText(massMatrix,"Matrices/MassMatrix.dat")
        !call ReadMatrixFromText(J,"Matrices/J.dat")

        call ReadIntegerMatrixFromText(vmapIMatrix,"Matrices/vmapI.dat")
        call ReadIntegerMatrixFromText(vmapOMatrix,"Matrices/vmapO.dat")
        call ReadIntegerMatrixFromText(vmapWMatrix,"Matrices/vmapW.dat")
        call ReadIntegerMatrixFromText(vmapMWMatrix,"Matrices/vmapMW.dat")
        call ReadIntegerMatrixFromText(vmapSMatrix,"Matrices/vmapS.dat")
        call ReadIntegerMatrixFromText(vmapNMatrix,"Matrices/vmapN.dat")
        call ReadIntegerMatrixFromText(vmapDMatrix,"Matrices/vmapD.dat")

        call ReadIntegerMatrixFromText(gaussmapBMatrix,"Matrices/gaussmapB.dat")
        call ReadIntegerMatrixFromText(gaussmapIMatrix,"Matrices/gaussmapI.dat")
        call ReadIntegerMatrixFromText(gaussmapOMatrix,"Matrices/gaussmapO.dat")
        call ReadIntegerMatrixFromText(gaussmapWMatrix,"Matrices/gaussmapW.dat")
        call ReadIntegerMatrixFromText(gaussmapMWMatrix,"Matrices/gaussmapMW.dat")
        call ReadIntegerMatrixFromText(gaussmapSMatrix,"Matrices/gaussmapS.dat")
        call ReadIntegerMatrixFromText(gaussmapNMatrix,"Matrices/gaussmapN.dat")
        call ReadIntegerMatrixFromText(gaussmapDMatrix,"Matrices/gaussmapD.dat")

        call ReadIntegerMatrixFromText(EToE,"Matrices/EToE.dat")
        call ReadIntegerMatrixFromText(EToV,"Matrices/EToV.dat")
        call ReadIntegerMatrixFromText(BCType,"Matrices/BCType.dat")

        call ReadMatrixFromText(VXMatrix,"Matrices/VX.dat")
        call ReadMatrixFromText(VYMatrix,"Matrices/VY.dat")
        
        ! Turns maps into vectors
        vmapM = pack(vmapMMatrix,.true.)
        vmapP = pack(vmapPMatrix,.true.)
        vmapI = pack(vmapIMatrix,.true.)
        vmapO = pack(vmapOMatrix,.true.)
        vmapW = pack(vmapWMatrix,.true.)
        vmapMW = pack(vmapMWMatrix,.true.)
        vmapS = pack(vmapSMatrix,.true.)
        vmapN = pack(vmapNMatrix,.true.)
        vmapD = pack(vmapDMatrix,.true.)

        gaussmapM = pack(gaussmapMMatrix,.true.)
        gaussmapP = pack(gaussmapPMatrix,.true.)

        gaussmapB = pack(gaussmapBMatrix,.true.)
        gaussmapI = pack(gaussmapIMatrix,.true.)
        gaussmapO = pack(gaussmapOMatrix,.true.)
        gaussmapW = pack(gaussmapWMatrix,.true.)
        gaussmapMW = pack(gaussmapMWMatrix,.true.)
        gaussmapS = pack(gaussmapSMatrix,.true.)
        gaussmapN = pack(gaussmapNMatrix,.true.)
        gaussmapD = pack(gaussmapDMatrix,.true.)

        VX = pack(VXMatrix,.true.)
        VY = pack(VYMatrix,.true.)

        Nodex = pack(NodexMatrix,.true.)
        Nodey = pack(NodeyMatrix,.true.)

        ! Number of elements
        K = size(x(1,:))

        ! Used for dividing by J
        divJ = 1/J
    end subroutine ReadAllMatrices

    subroutine ConvertToParaView(Q,timestep,isLastTimeStep)
        ! Converts the file data to csv format readable by ParaView and then saves it to file
        ! Arguments:    Q       -   Conserved vector

        ! Author: Sam MacPherson

        implicit none

        ! Define Arguments
        real(dp), dimension(:,:,:), intent(in)  :: Q
        integer, intent(in)                     :: timestep
        logical, intent(in)                     :: isLastTimeStep
        ! Primitive variables       
        real(dp), dimension(size(Q(:,1,1)), size(Q(1,:,1))  )    :: rho, uvel, vvel, p, T

        ! Result matrix
        real(dp), dimension(:,:), allocatable   :: result

        ! String which is appended on the end of filename
        character (6)   :: timestr
        character (50)   :: filename
        character (10)  :: fmt
        ! Loop counters
        integer :: i,j,row

        ! Shape
        integer, dimension(3)   :: shapeQ

        write (timestr,"(I0.6)") timestep

        if (isLastTimeStep) then
            filename = "Results/results.csv.last"
        else
            filename = "Results/results.csv."//trim(timestr)
        end if

        filename = trim(filename)
        ! Shape of Q
        shapeQ = shape(Q)

        ! Extract primitive variables
        call ConToPrim(Q,rho,uvel,vvel,p)
        call GetTemperature(p,rho,T)
        allocate(result(2*size(rho),9))

        row = 1
        do i = 1,shapeQ(1)
            do j = 1,shapeQ(2)
                !x,y,z=0,rho,u,v,p,vel_mag
                result(row,1) = x(i,j)
                result(row,2) = y(i,j)
                result(row,3) = 0
                result(row,4) = rho(i,j)
                result(row,5) = uvel(i,j)
                result(row,6) = vvel(i,j)
                result(row,7) = p(i,j)
                result(row,8) = sqrt(uvel(i,j)**2+vvel(i,j)**2)
                result(row,9) = T(i,j)

                row = row + 1
                !x,y,z=1,rho,u,v,p
                result(row,1) = x(i,j)
                result(row,2) = y(i,j)
                result(row,3) = 1
                result(row,4) = rho(i,j)
                result(row,5) = uvel(i,j)
                result(row,6) = vvel(i,j)
                result(row,7) = p(i,j)
                result(row,8) = sqrt(uvel(i,j)**2+vvel(i,j)**2)
                result(row,9) = T(i,j)

                row = row + 1
            end do
        end do
       
        open(12,file = filename)
            write(12,*) "x",",", "y",",", "z",",", "rho",",", "u",",", "v",",", "p", ",", "vel_mag", ",","T"
            do i = 1, size(result(:,1))
                write(12,'(*(E22.15 : ", "))') (result(i,j), j=1,9)
            end do
        close(12)

    end subroutine ConvertToParaview

    subroutine ReadIniFile()
        ! Read the settings and parameter file
        ! Also set the cubature and gauss order

        ! Author: Sam MacPherson

        open(12,file = "SettingsAndParameters.ini")
            read(12,*) N
            read(12,*) gammaGas
            read(12,*)
            read(12,*) totalTime
            read(12,*) dt
            read(12,*) CFL
            read(12,*) outputFrequency
            read(12,*) RKOrder
            read(12,*) fluxScheme
            read(12,*) residualStop
            read(12,*) GasConstant
            read(12,*) mu0
            read(12,*) T0
            read(12,*) S
            read(12,*) constMu
            read(12,*) thermalConductivity
            read(12,*) useLimiter
            read(12,*) BCScheme
            read(12,*) testProblem
        close(12)

        ! Set cubature and Gauss order
        Np = (N+1)*(N+2)/2
        NGauss = ceiling(3*(N+1)/2)
        Nfp = N + 1

        ! Decide whether to use Sutherland law, or use constant viscosity.
        if (constMu.eq.0) then
            useSutherland = .true.
        else
            useSutherland = .false.
        end if
        
    end subroutine ReadIniFile





end Module