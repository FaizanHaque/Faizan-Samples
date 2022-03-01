function ops = p_diff_operators(mesh)
%
% CDDiscretize creates vectors which will be used to create spherical
% derivative operators in cartesian coordinates
%
% Author: 
%   - Syed Faizanul Haque, Francois Primeau
%
% Use as:
%
% Inputs:
%
%   mesh
nX = mesh.nX; nY = mesh.nY; nZ = mesh.nZ;
dX = mesh.dX; dY = mesh.dY; dZt = mesh.dZt;dZw = mesh.dZw;
dA = mesh.dA;
d0 = @(x) spdiags(x(:),0,length(x(:)),length(x(:)));

%Setup the multiplying factors for the operators
II = zeros (nY,nX,nZ);
II(:) = 1:nY*nX*nZ;

%The lower case are the indices
iE = II (:, [ 2:end,1],:);
iW = II (:, [ end, 1:end-1],:);
iN = II ( [ 2:end,1],:,:);
iS = II ( [ end, 1:end-1],:,:);
iD = II (:,:, [ 2:end,1]);
iU = II (:,:, [ end, 1:end-1]);

%Create a sparse identity matrix (only stores the nonzero values)
I = speye(nX*nY*nZ); %

%
IE= I(iE(:),:);
IW= I(iW(:),:);
IN= I(iN(:),:);
IS= I(iS(:),:);
IU= I(iU(:),:);
ID= I(iD(:),:);


%Averaging operators
QxAQy = (I+IE+IS+IS*IE)/4;
QyAQx = (I+IN+IW+IW*IN)/4;
QyAu =(I+IN)/2;
QxAu = (I+IE)/2;

%Difference operators
% (i,j,k)--> (y,x,z)
BDi = (I - IS); BDY = d0(1./dY(:))*BDi;
FDi = (IN - I); FDY = d0(1./dY(:))*FDi;
BDj = (I - IW); BDX = d0(1./dX(:))*BDj;
FDj = (IE - I); FDX = d0(1./dX(:))*FDj;
FDk = (IU - I); FDZ = d0(1./dZt(:))*FDk;
BDk = (I - ID); BDZ = d0(1./dZw(:))*BDk;
% Derivitive operators
DIV2D = d0(1./dA)*[BDj*d0(dY),BDi*d0(dX)];
GRAD2D = -DIV2D';
DEL2H = DIV2D*GRAD2D; % horizontal Laplacian
D2Z = FDZ*BDZ;
% glue everything together
ops.BDX = BDX;
ops.FDX = FDX;
ops.BDY = BDY;
ops.FDY = FDY;
ops.FDZ = FDZ;
ops.BDZ = BDZ;
%
ops.DIV2D = DIV2D;
ops.GRAD2D = GRAD2D;
%
ops.D2Z = D2Z;
ops.DEL2H = DEL2H;
%
ops.QxAQy = QxAQy;
ops.QyAQx = QyAQx;
ops.QxAu = QxAu;
ops.QyAu = QyAu;
