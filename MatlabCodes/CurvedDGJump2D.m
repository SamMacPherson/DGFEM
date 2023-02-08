function [jumpU] = CurvedDGJump2D(gU, gmapD, bcU)

% function [jumpU] = CurvedDGJump2D(gU, gmapD, bcU)
% purpose: compute discontinuous Galerkin jump applied
%          to a field given at cubature and Gauss nodes

Globals2D;

% surface traces at Gauss nodes
gUM = gU(gauss.mapM); gUP = gU(gauss.mapP); gUP(gmapD) = bcU(gmapD);

% compute jump term and lift to triangle interiors
fx = gUM - gUP;
fx(1,2463)/2;
jumpU = -gauss.interp'*(gauss.W.*fx);

% multiply straight sided triangles by inverse mass matrix




return
