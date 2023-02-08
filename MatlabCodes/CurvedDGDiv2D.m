function [divU] = CurvedDGDiv2D(cU, cV, gU, gV, gmapN, bcNdotU)
format long
% function [divU] = CurvedDGDiv2D(cU, cV, gU, gV, gmapN, bcNdotU)
% Purpose: compute the divergence of a vectorial function given at cubature
%          and surface Gauss nodes

Globals2D;
a = cub.W.*(cub.rx.*cU+cub.ry.*cV);

cub.Dr(:,1)
cU(:,2463)
% volume terms: U
divU = (cub.Dr')*(cub.W.*(cub.rx.*cU+cub.ry.*cV)) + ...
          (cub.Ds')*(cub.W.*(cub.sx.*cU+cub.sy.*cV));

divU(1,2463);
% surface traces at Gauss nodes
gUM = gU(gauss.mapM); gUP = gU(gauss.mapP); 
gVM = gV(gauss.mapM); gVP = gV(gauss.mapP);

% normal fluxes
gFxM = gauss.nx.*gUM + gauss.ny.*gVM;
gFxP = gauss.nx.*gUP + gauss.ny.*gVP;

% Apply boundary conditions
gFxP(gmapN) = bcNdotU(gmapN);

% add flux terms to divergence
divU = divU - gauss.interp'*(gauss.W.*(gFxM + gFxP))/2;

a = gauss.interp'*(gauss.W.*(gFxM + gFxP))/2;
a(1,2463);

% multiply straight sided triangles by inverse mass matrix




% Correct sign
divU = -divU;

return
