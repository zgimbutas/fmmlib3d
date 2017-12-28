% Debugging output to screen and fort.13
fmm3dprini(6,13)

% Source triangle
VF = [0 0 0;-1 1 0;-0.2 -0.2 0.2];
N = [1 0 0];
C = mean(VF');
% Target grid
[X,Y] = meshgrid(linspace(-1,1),linspace(-1,1));
Z = zeros(size(X));
% Target positions
Q = [X(:) Y(:) Z(:)];

iprec = 3;
UFMM = lfmm3dtria(iprec,1,VF,N',C',0,0,1,1,N',0,0,size(Q,1),Q',1,0);
UDIR = l3dtriadirect   (1,VF,N',C',0,0,1,1,N',0,0,size(Q,1),Q',1,0);

UDIFF = real(UFMM.pottarg - UDIR.pottarg);
UREL = abs(UDIFF)./abs(UDIR.pottarg);
err_abs = max(abs(UDIFF))
err_rel = max(abs(UREL))

err = reshape(UDIR.pottarg - UFMM.pottarg, 100, 100);

err = reshape(UFMM.pottarg, 100, 100);
imagesc(log10(abs(err)))
colorbar