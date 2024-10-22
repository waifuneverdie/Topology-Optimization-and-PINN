%%%% A 110 LINE TOPOLOGY OPTIMIZATION CODE WITH HEAVISIDE FILTERING Nov, 2010%%%%
function top111(nelx,nely,volfrac,penal, penal_c, rmin, rmin_c,lambda_e, ft)
%% Define new parameters
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE1 = lambda_e *1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
KE2 = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
edofVec = reshape(2*nodenrs(1:end-1, 1:end-1) + 1, (nelx)*(nely), 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1], (nelx)*(nely), 1);

iK1 = reshape(kron(edofMat,ones(8,1))',64*(nelx)*(nely),1);
jK1 = reshape(kron(edofMat,ones(1,8))',64*(nelx)*(nely),1);
iK2 = reshape(kron(edofMat,ones(8,1))',64*(nelx)*(nely),1);
jK2 = reshape(kron(edofMat,ones(1,8))',64*(nelx)*(nely),1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1/2,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH1 = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH1 = ones(size(iH1));
sH1 = zeros(size(iH1));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH1(k) = e1;
        jH1(k) = e2;
        sH1(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H1 = sparse(iH1,jH1,sH1);
Hs1 = sum(H1,2);

iH2 = ones(nelx*nely*(2*(ceil(rmin_c)-1)+1)^2,1);
jH2 = ones(size(iH2));
sH2 = zeros(size(iH2));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin_c)-1),1):min(i1+(ceil(rmin_c)-1),nelx)
      for j2 = max(j1-(ceil(rmin_c)-1),1):min(j1+(ceil(rmin_c)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH2(k) = e1;
        jH2(k) = e2;
        sH2(k) = max(0,rmin_c-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H2 = sparse(iH2,jH2,sH2);
Hs2 = sum(H2,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
beta = 1;
eta = 0.1;
if ft == 1 || ft == 2
  xPhys = x;
elseif ft == 3
  xTilde = x;
  xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
  
  x2 = xPhys;
  x2Tilde = x2;
  x2Phys = (tanh(beta*eta)+tanh(beta*(x2Tilde-eta))) / (tanh(beta*eta)+tanh(beta*(1-eta)));
  
  xSubs = xPhys;
  xCoat = x2Phys .* (1 - xPhys);
end
loopbeta = 0;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01 %&& loop ~=205
  loopbeta = loopbeta+1;
  loop = loop+1;
  %% FE-ANALYSIS
  sK1 = reshape(KE1(:)*(Emin+xSubs(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  sK2 = reshape(KE2(:)*(Emin+xCoat(:)'.^penal_c*(E0-Emin)),64*nelx*nely,1);
  K1 = sparse(iK1,jK1,sK1); K1 = (K1+K1')/2;
  K2 = sparse(iK2,jK2,sK2); K2 = (K2+K2')/2;
  K = K1+K2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce1 = reshape(sum((U(edofMat)*KE1).*U(edofMat),2),nely,nelx);
  c1 = sum(sum((Emin+xSubs.^penal*(E0-Emin)).*ce1));
  dc1 = -penal*(E0-Emin)*xSubs.^(penal-1).*ce1;
  dv1 = ones(nely,nelx);

  ce2 = reshape(sum((U(edofMat)*KE2).*U(edofMat),2),nely,nelx);
  c2 = sum(sum((Emin+xCoat.^penal_c*(E0-Emin)).*ce2));
  dc2 = -penal_c*(E0-Emin)*xCoat.^(penal_c-1).*ce2;
  dv2 = ones(nely,nelx);

  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc1(:) = H1*(x(:).*dc1(:))./Hs1./max(1e-3,x(:));
  elseif ft == 2
    dc1(:) = H1*(dc1(:)./Hs1);
    dv1(:) = H1*(dv1(:)./Hs1);
  elseif ft == 3
    dx1 = beta*exp(-beta*xTilde)+exp(-beta);
    dc1(:) = H1*(dc1(:).*dx1(:)./Hs1);
    dv1(:) = H1*(dv1(:).*dx1(:)./Hs1);

    dx2 = beta*(1-tanh(beta*(x2Tilde-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));
    dc2_a(:) = (1-xPhys(:)).*H2*(dc2(:).*dx2(:)./Hs2);
    dc2_b(:) = -x2Phys(:) .* dc2(:);
    dc2(:) = dc2_a(:) + dc2_b(:);
    dv2_a(:) = (1-xPhys(:)).*H2*(dv2(:).*dx2(:)./Hs2);
    dv2_b(:) = -x2Phys(:) .* dv2(:);
    dv2(:) = dv2_a(:) + dv2_b(:);

    dx2 = beta*exp(-beta*xTilde)+exp(-beta);
    dc2(:) = H1*(dc2(:).*dx2(:)./Hs1);
    dv2(:) = H1*(dv2(:).*dx2(:)./Hs1);

    dc = dc1 + dc2;
    dv = dv1;

  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H1*xnew(:))./Hs1;
    elseif ft == 3
      xTilde(:) = (H1*xnew(:))./Hs1;
      xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);

      x2 = xPhys;
      x2Tilde = (H2*x2(:))./Hs2; 
      x2Phys = (tanh(beta*eta)+tanh(beta*(x2Tilde-eta))) / (tanh(beta*eta)+tanh(beta*(1-eta)));
      
      xSubs = xPhys;
      xCoat(:) = x2Phys(:) .* (1 - xPhys(:));
    end
    if sum(xSubs(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c1, ...
    mean(xSubs(:)),change);
%% PLOT DENSITIES

% Maximize figure window
set(gcf, 'Position', get(0, 'Screensize'));

% Clear the current axes
cla;

% Plot xCoat with full opacity
colormap(gray);
h = imagesc(1-xCoat);
set(h, 'AlphaData', 1);  % Full opacity for xCoat
caxis([0 1]); 
axis equal; 
axis off; 
hold on;

% Overlay xSubs with transparency
l = imagesc(1-xSubs); 
set(l, 'AlphaData', 0.25);  % Set transparency for overlay

title('xSubs and xCoat overlayed');
drawnow;
% colormap(gray); imagesc(1-xCoat); caxis([0 1]); axis equal; axis off; drawnow;
  %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
  if ft == 3 && beta < 32 && (loopbeta >= 50 || change <= 0.01)
    beta = 2*beta;
    loopbeta = 0;
    change = 1;
    eta = eta;
    fprintf('Parameter beta increased to %g.\n',beta);
  end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
