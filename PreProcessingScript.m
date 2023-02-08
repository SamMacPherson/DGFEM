% Script that computes all pre-processing for a 2D Compressible Navier-Stokes DG-FEM
% This script uses the Hesthaven and Warbuton Nodal DG-FEM code for the pre-processing

% Author: Sam MacPherson


% Add Folders to Matlabs path
pw = pwd;
addpath([pw,'/MatlabCodes']);
addpath([pw,'/Grids']);
disp("Codes added to path");
disp(" ");


% Namespace with all the variables
Globals2D;

% Read in setup file
fid = fopen('SettingsAndParameters.ini');

N = fscanf(fid,'%g');

% Skip to mesh line
tline = fgets(fid);
tline = fgets(fid);

% Need all the alphabet as it tells fscanf the allowable chars, it stops when it hits a blank space
meshFileName = fscanf(fid,'%[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ._1234567890]');

disp(meshFileName)

% Set polynomial order
CubatureOrder = 3*(N+1);
NGauss = ceil(3*(N+1)/2);

% Read in mesh
% Mesh should be in Gambit .neu format
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(meshFileName);
disp("Mesh read in");
disp(" ");

% Initialize solver and construct grid and metric
StartUp2D;
disp("Constructed matrices and geometric factors");
disp(" ");





% Build maps for BC's
BuildBCMaps2D;
disp("Built BC maps");
disp(" ");



% Build cubature information
cub = CubatureVolumeMesh2D(CubatureOrder);

% Build Gauss node data
gauss = GaussFaceMesh2D(NGauss);
disp("Generated cubature and Gauss information");
disp(" ");

% If the following are empty then Fortran will throw an error when loading in the matrix
% So if they're empty, give them a dummy matrix that will never be used
if(isempty(vmapI))
    vmapI = [1 1];
end

if(isempty(vmapO))
    vmapO = [1 1];
end

if(isempty(vmapW))
    vmapW = [1 1];
end

if(isempty(vmapMW))
    vmapMW = [1 1];
end

if(isempty(vmapS))
    vmapS = [1 1];
end

if(isempty(vmapN))
    vmapN = [1 1];
end

if(isempty(vmapD))
    vmapD = [1 1];
end



if(isempty(gauss.mapI))
    gauss.mapI = [1 1];
end

if(isempty(gauss.mapO))
    gauss.mapO = [1 1];
end

if(isempty(gauss.mapW))
    gauss.mapW = [1 1];
end

if(isempty(gauss.mapMW))
    gauss.mapMW = [1 1];
end

if(isempty(gauss.mapS))
    gauss.mapS = [1 1];
end

if(isempty(gauss.mapN))
    gauss.mapN = [1 1];
end

if(isempty(gauss.mapD))
    gauss.mapD = [1 1];
end


% Save matrices to be read by Fortran
SaveMatrix(x,"Matrices/x");
SaveMatrix(y,"Matrices/y");

SaveMatrix(VX,"Matrices/VX");
SaveMatrix(VY,"Matrices/VY");

SaveMatrix(BCType,"Matrices/BCType");

SaveMatrix(V,"Matrices/V");
SaveMatrix(MassMatrix,"Matrices/MassMatrix");
SaveMatrix(J,"Matrices/J");
SaveMatrix(Fscale,"Matrices/Fscale");
SaveMatrix(vmapM,"Matrices/vmapM");
SaveMatrix(vmapP,"Matrices/vmapP");
SaveMatrix(cub.Dr,"Matrices/cubDr");
SaveMatrix(cub.Ds,"Matrices/cubDs");

SaveMatrix(cub.rx,"Matrices/cubrx");
SaveMatrix(cub.ry,"Matrices/cubry");
SaveMatrix(cub.sx,"Matrices/cubsx");
SaveMatrix(cub.sy,"Matrices/cubsy");

SaveMatrix(cub.V,"Matrices/cubV");
SaveMatrix(cub.W,"Matrices/cubW");

SaveMatrix(gauss.interp,"Matrices/gaussinterp");

SaveMatrix(gauss.nx,"Matrices/gaussnx");
SaveMatrix(gauss.ny,"Matrices/gaussny");

SaveMatrix(gauss.mapB,"Matrices/gaussmapB");
SaveMatrix(gauss.mapM,"Matrices/gaussmapM");
SaveMatrix(gauss.mapP,"Matrices/gaussmapP");

SaveMatrix(gauss.W,"Matrices/gaussW");

SaveMatrix(vmapI,"Matrices/vmapI");
SaveMatrix(vmapO,"Matrices/vmapO");
SaveMatrix(vmapW,"Matrices/vmapW");
SaveMatrix(vmapMW,"Matrices/vmapMW");
SaveMatrix(vmapS,"Matrices/vmapS");
SaveMatrix(vmapN,"Matrices/vmapN");
SaveMatrix(vmapD,"Matrices/vmapD");

SaveMatrix(gauss.mapB,"Matrices/gaussmapB");
SaveMatrix(gauss.mapI,"Matrices/gaussmapI");
SaveMatrix(gauss.mapO,"Matrices/gaussmapO");
SaveMatrix(gauss.mapW,"Matrices/gaussmapW");
SaveMatrix(gauss.mapMW,"Matrices/gaussmapMW");
SaveMatrix(gauss.mapS,"Matrices/gaussmapS");
SaveMatrix(gauss.mapN,"Matrices/gaussmapN");
SaveMatrix(gauss.mapD,"Matrices/gaussmapD");

SaveMatrix(EToE,"Matrices/EToE");
SaveMatrix(EToV,"Matrices/EToV");

disp("Matrices saved to be imported to Fortran");
disp(" ");


