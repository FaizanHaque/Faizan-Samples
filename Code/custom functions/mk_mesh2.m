function mesh = mk_mesh2(d)
Xw = d.lon(d.rj); % longitudes of the w points
Yw = d.lat(d.ri); % latitudes of the w points
Zw = d.zw(1:(d.rk(end)+1)); % depths of the w points
Zt = d.zt(1:(d.rk(end)+1));
%
dlat = Yw(2)-Yw(1);
dlon = Xw(2)-Xw(1);
dzt = d.dzt(1:d.rk(end)+1);
dzw = d.dzw(1:d.rk(end)+1);
%                          
nX = length(Xw); %length of x vector
nY = length(Yw); %length of y vector
nZ = length(Zw); %length of Z vector            
% taking the f-plane approximation
%dX = ones(nY,nX,nZ)*dlon; % meters
%dY = ones(nY,nX,nZ)*dlat; %           
phi0 = pi*mean(Yw)/180;

%dX = ones(nY,nX,nZ)*dlon*111.699*1000*cos(phi0); % meters
%dY = ones(nY,nX,nZ)*dlat*111.699*1000; %                                     
dX = ones(nY,nX,nZ)*dlon*111.699*1000.*cos(pi*Yw/180); % allowing the dX to vary with latitude
%keyboard
dY = ones(nY,nX,nZ)*dlat*111.699*1000; %                                     


% dX = ones(nY,nX,nZ)*dlon*a*cos(phi0); % meters
% dY = ones(nY,nX,nZ)*dlat*a; %                                     
dZw = zeros(1,1,nZ); 
dZt = zeros(1,1,nZ); 
dZt(1,1,:) = dzt(1:d.rk(end)+1);
dZw(1,1,:) = dzw(1:d.rk(end)+1);
dZt = dZt(ones(nY,1),ones(nX,1),:); % copy to all horizontal point
dZw = dZw(ones(nY,1),ones(nX,1),:); % copy to all horizontal point
dV = dY.*dX.*dZt;
dA = dX.*dY;
%
[mesh.X,mesh.Y,mesh.Zw] = meshgrid(Xw,Yw,Zw);
[mesh.X,mesh.Y,mesh.Zt] = meshgrid(Xw,Yw,Zt);

mesh.phi0 = phi0;
mesh.dZt = dZt;
mesh.dZw = dZw;
mesh.dX = dX;
mesh.dY = dY;
mesh.dA = dA;
mesh.dV = dV;
mesh.nX = nX;
mesh.nY = nY;
mesh.nZ = nZ;
mesh.dlat = dlat;
mesh.dlon = dlon;
mesh.dzt = dzt;
mesh.dzw = dzw;