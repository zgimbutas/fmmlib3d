% Debugging output to screen and fort.13
fmm3dprini(6,13)

verts = [0 0 0;-1 1 0;-0.2 -0.2 0.2];
ifaces = [1 2 3]';
[verts,ifaces]=atrirefine(verts,ifaces);
[verts,ifaces]=atrirefine(verts,ifaces);
[verts,ifaces]=atrirefine(verts,ifaces);
[verts,ifaces]=atrirefine(verts,ifaces);

[triaflat,trianorm,triaarea,centroid,triatang1,triatang2] = ...
     atriproc(verts,ifaces);

VF = triaflat;
N = trianorm;
C = centroid;

% Target grid
[X,Y] = meshgrid(linspace(-1,1),linspace(-1,1));
Z = zeros(size(X));
% Target positions
Q = [X(:) Y(:) Z(:)];

nsource = size(C,2);
charge = zeros(nsource);
dipstr = ones(nsource);

iprec = 3;
UFMM = lfmm3dtria(iprec,nsource,VF,N,C,0,0,1,dipstr,N,0,0,size(Q,1),Q',1,0);
UDIR = l3dtriadirect   (nsource,VF,N,C,0,0,1,dipstr,N,0,0,size(Q,1),Q',1,0);

UDIFF = real(UFMM.pottarg - UDIR.pottarg);
UREL = abs(UDIFF)./abs(UDIR.pottarg);
err_abs = max(abs(UDIFF))
err_rel = max(abs(UREL))

err = reshape(UDIR.pottarg - UFMM.pottarg, 100, 100);

err = reshape(UFMM.pottarg, 100, 100);
imagesc(log10(abs(err)))
colorbar