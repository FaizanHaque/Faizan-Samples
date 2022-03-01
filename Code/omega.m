%
% Load in the data
%

%fname = '../Data/dataset-armor-3d-rep-monthly_1553733806594_Seo2005GS_dlMar27.nc'
%fname = '../Data/dataset-armor-3d-rep-weekly_1553712920759.nc';
fname = '../../../Data/armor-3d-rep-monthly_sep_2005_gulfstream_pascual_et_al_2015.nc';
%fname = '../../Data/armor-3d-rep-monthly_sep_2005_gulfstream_pascual_et_al_2015.nc';
data = p_get_data(fname);
u0 = data.u(:,:,:,1);
v0 = data.v(:,:,:,1);
T0 = data.temp(:,:,:,1);
S0 = data.salt(:,:,:,1);
z0 = -data.depth;

%
% Define the size and location of the domain
%
% [i,j,k] <-> [Y,X,Z]
rj = 1:136; % range of grid boxes in X direction 
ri = 1:49;% range of grid boxes in Y direction
rk = 1:20; % range of grid boxes in Z direction
fNo=605
%61 15:120,9:40
% make a new mesh
%

nz = 33;
[zt,zw,dzt,dzw] = zgrid(nz);

domain.lat = data.lat;   domain.lon = data.lon(6:end);
domain.dlat = data.dlat; domain.dlon = data.dlon;
domain.depth = data.depth;
zw = zw(2:end-1)'; zt = zt(2:end-1)'; dzw = dzw(2:end-1)'; dzt = dzt(2:end-1)';
domain.zt = zt; domain.zw = zw; domain.dzt = dzt; domain.dzw = dzw;
domain.ri = ri; domain.rj = rj; domain.rk = rk;
mesh = mk_mesh2(domain);
% interpolate onto a new mesh
[u, v, T, S ] = regridZData(u0, v0,T0,S0,z0,zt,ri,rj,rk);
u = u(ri,rj,1:rk(end)+1);
v = v(ri,rj,1:rk(end)+1);

%
% buoyancy frequency
%
[LON,LAT,DEPTH] = meshgrid(domain.lon,domain.lat,zw);

P = 0*mesh.Zw;
P(:) = gsw_p_from_z(mesh.Zw(:),mesh.Y(:));
N2 = 0*P;

%i = 24;
%j = 37;
% [n2, p_mid] = gsw_Nsquared(squeeze(S(i,j,1:21)),...
%                                    squeeze(T(i,j,1:21)),...
%                                    squeeze(P(i,j,1:21)),...
%                                    squeeze(mesh.Y(i,j,1:21)));

for j = 1:size(P,2)
    for i = 1:size(P,1)
        [n2, p_mid] = gsw_Nsquared(squeeze(S(i,j,1:21)),...
                                   squeeze(T(i,j,1:21)),...
                                   squeeze(P(i,j,1:21)),...
                                   squeeze(mesh.Y(i,j,1:21)));
            N2(i,j,:) = interp1(p_mid,n2,squeeze(P(i,j,1:21)));
    end
end
%N2 = N2(:,:,1:20);
% Coriolis frequency
oneday = 24*60^2; %seconds
Omega = 2*pi/oneday; %rotation/sec
f = 0*N2;
for l = 1:size(mesh.Y,1)
%f(l,:,:) = 2*Omega*sin(pi/180*mesh.Y(l,:,1:21));

    f(l,:,:) = 2*Omega*sin(mesh.phi0);
end
f2 = f.*f;


% Make the finite difference operators
%
ops = p_diff_operators(mesh);

%
D2Z = ops.D2Z;                          
DEL2H = ops.DEL2H;
BDX = ops.BDX;
FDX = ops.FDX;
BDY = ops.BDY;
FDY = ops.FDY;
FDZ = ops.FDZ;
BDZ = ops.BDZ;
DIV2D = ops.DIV2D;
%
QxAQy = ops.QxAQy; 
QyAQx = ops.QyAQx; 
QxAu = ops.QxAu; 
QyAu = ops.QyAu; 

% Make the elliptic operator for the Omega equation
%
%keyboard
d0 = @(x) spdiags(x(:),0,length(x(:)),length(x(:)));
%L = f^2*D2Z + DEL2H*d0(N2(:));
L = 0*D2Z;
L(:) = d0(f2(:))*D2Z + DEL2H*d0(N2(:));
%L(:) = d0(f2(:))*D2Z + d0(N2(:))*DEL2H;
    
%N = 1e-5;

%L = f^2*D2Z + (N^2)*DEL2H;

% impose lateral and vertical boundary conditions
% B=1 for boundary points
% B=0 for interior points
B = 0*mesh.dX;
B(1,:,:) = 1; B(end,:,:) = 1;  % y boundaries
B(:,:,1) = 1; B(:,:,end) = 1;  % z boundaries
B(:,1,:) = 1; B(:,end,:) = 1;  % x boundaries
ib = find(B(:)); 
L(ib,:) = 0;
L(ib,ib) = speye(length(ib));




dudx =0*mesh.X; dvdx = 0*mesh.X;
dudy =0*mesh.X; dvdy = 0*mesh.X;
dudzI =0*mesh.X; dvdzI = 0*mesh.X;

%
dudx(:) = FDX*u(:); dvdx(:) = FDX*v(:);
dudy(:) = FDY*u(:); dvdy(:) = FDY*v(:);
dudz(:) = FDZ*u(:); dvdz(:) = FDZ*v(:);
%dudzI = permute(dudzI,[3 1 2]);

%dvdzI = permute(dvdzI,[3 1 2]);

%dudz = 0*dudzI;dvdz = dvdzI*0; 

%for j = rj
 %   for i = ri
  %      dudz(:,i,j) = interp1(zt(1:21),dudzI(:,i,j),zw(1:21));
   %     dvdz(:,i,j) = interp1(zt(1:21),dvdzI(:,i,j),zw(1:21));
    %end    
%end
%dudz = permute(dudz ,[2 3 1]);
%dvdz = permute(dvdz ,[2 3 1]);



x_new = mesh.X + data.dlon/2;% - data.lon(1);
y_new = mesh.Y + data.dlat/2;% - data.lon(1);

Dudz = 0*u;
Dvdz = 0*u;
Dudx = 0*u;
Dvdx = 0*u;
Dudy = 0*u;
Dvdy = 0*u;

for j = 1:length(rj)
   for i = 1:length(ri)
       Dvdz(i,j,:) = interp1(squeeze(mesh.Zt(i,j,:)),squeeze(dvdz(i,j,:)),squeeze(mesh.Zw(i,j,:)));
       Dudz(i,j,:) = interp1(squeeze(mesh.Zt(i,j,:)),squeeze(dudz(i,j,:)),squeeze(mesh.Zw(i,j,:)));
   end
%    for k = 1:rk(end)+1
%        Dvdy(:,j,k) = interp1(squeeze(y_new(:,j,k)),squeeze(dvdy(:,j,k)),squeeze(mesh.Y(:,j,k)));
%        Dudy(:,j,k) = interp1(squeeze(y_new(:,j,k)),squeeze(dudy(:,j,k)),squeeze(mesh.Y(:,j,k)));
%    end    
end
dudz = Dudz;
dvdz = Dvdz;
% for k = 1:rk(end)+1
%     for i = ri
%         Dvdx(i,:,k) = interp1(squeeze(x_new(i,:,k)),squeeze(dvdx(i,:,k)),squeeze(mesh.X(i,:,k)));
%         Dudx(i,:,k) = interp1(squeeze(x_new(i,:,k)),squeeze(dudx(i,:,k)),squeeze(mesh.X(i,:,k)));
%     end
% end


% Vertical Velocity- Original Version
kBot = 19;
Qx = 0*Dudx;
Qy = 0*Dudx;
QX = 0*Dudx;
QY = 0*Dudx;

DivQ = 0*Dudx;
% Ref. Pascual et al. (2015)

Qx(:) =  2*f(:).*(dvdx(:).*(QxAu*dudz(:)) + (QxAQy*dvdy(:)).*(QxAu*dvdz(:)));
Qy(:) = -2*f(:).*((QyAQx*dudx(:)).*(QyAu*dudz(:)) + dudy(:).*(QyAu*dvdz(:))); 

%QX(:) =  2*f(:).*(Dvdx(:).*Dudz(:) + Dvdy(:).*Dvdz(:));
%QY(:) = -2*f(:).*(Dudx(:).*Dudz(:) + Dudy(:).*Dvdz(:)); 
Q = [Qx(:); Qy(:)];
DivQ(:) = DIV2D*Q;

%Mask out points shallower than 1000m
[sX,sY,sZ] = find(isnan(DivQ(:,:,kBot)));
for i = 1:length(sX)
DivQ(sX(i),sY(i),1:kBot) = -99999; %Mark every layer from surface to z=-1000m for the above mentioned points
end
is = find(DivQ == -99999);
DivQ(is) = 0;
L(is,:) = 0;
L(is,is) = speye(length(is));

%Find all nans in DivQ (land/missing data)and set the w=0 there
in = find(isnan(DivQ));
DivQ(in) = 0;
L(in,:) = 0;
L(in,in) = speye(length(in));

% impose the boundary condition on w
DivQ(ib) = 0;
w = 0*DivQ;
w(:) = L\DivQ(:);
%put back land points and the shallow ocean floor points
w(in) = nan;
w(is) = nan;
%u(is) = nan;
u(in) = nan;
%v(is) = nan;
v(in) = nan;
% for i = 1:length(mX)
% w(sX(i),sY(i),1:19) = nan;
% end
%w = w*100*86400;%convert m/s to cm/day

%keyboard

%latitudonally averaged profile
%uLat = zeroes(1,size(u,2),size(u,3));
%wLat = 0*uLat;
uLat = mean(u(:,:,:),1);
uLat = permute(uLat,[3 2 1]);
wLat = mean(w(:,:,:),1);
wLat = permute(wLat,[3 2 1]);
% figure(8),
% qLat=quiver(domain.lon(1:end-1),zw(1:21),uLat,wLat,'k');
% qLat.MaxHeadSize=1;
% xticks([-72 -66 -60 -54 -48 -42 ])
% xticklabels({'72°','66°', '60°', '54°' ,'48°', '42°' })


%Estimated vertival velocity: W~Ro*V*H/L=V^2*H/fL^2
%set L~100km , H~1km, f~0.855*10^-4 s^-1, v~data~ sqrt(u^2+v^2)
%pg 43 eq 2.109 in Jim McWilliams book
L = 100000;%m
H= 1000;%m
f= 0.855*10^-4;%s^-1
V_Ro = 0*v;
V_Ro(:)=sqrt(u(:).^2+v(:).^2);
W_Ro = (V_Ro.^2)*H/(f*L^2); 
caxis([-200 200]);


d = 8; % level 8 corresponds to 100 m depth




figure(fNo)
%subplot(2,1,1)
%for d = 2:30    
%clevs = [-10:1:-2,2:1:10];
clevs = -30:.1:30;
%contourf(mesh.X(:,:,1)-360, mesh.Y(:,:,1), HEIGHT(:,:,1), 'edgecolor', 'none');

[C_cont,h_cont] = contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (w(:,:,d))*86400,clevs, 'edgecolor', 'none');
%h_cont = contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (w(:,:,d))*86400, 'edgecolor', 'none');
clabel(C_cont,h_cont);
colormap(redblue)

%caxis([-10 10]);
%colormap(jet(256));
hold on
qScale = 2.5;
q = quiver(domain.lon(rj)-360, data.lat(ri), u(:,:,d), v(:,:,d),qScale,'k');
%q.Color = 'black';  
dNo = num2str(d);
dVal = num2str(DEPTH(1,1,d));
    title(['W at  z =' dVal 'm ']);
    xlabel('Longitude');
    ylabel('Latitiude');
    hcb = colorbar;
    set(get(hcb,'Title'),'String',' W [m/day]');
xticks([-72 -66 -60 -54 -48 -42 ])
xticklabels({'72°','66°', '60°', '54°' ,'48°', '42°' })
    wMax = max(max(w(:,:,d)))*86400;
    caxis([-wMax wMax]);
   dim = [.2 .5 .3 .3];
str = ['Max w = ', num2str(wMax), 'm/day'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
figure(7)

%quiver3(LON(1:10:end-1,1:10:end-1,3:2:21)-360,LAT(1:10:end-1,1:10:end-1,3:2:21),DEPTH(1:10:end-1,1:10:end-1,3:2:21),u(1:10:end,1:10:end,3:2:21),v(1:10:end,1:10:end,3:2:21),w(1:10:end,1:10:end,3:2:21),01,'k');
quiver3(LON(1:10:end-1,1:10:end-1,5:2:21)-360,LAT(1:10:end-1,1:10:end-1,5:2:21),DEPTH(1:10:end-1,1:10:end-1,5:2:21),0*u(1:10:end,1:10:end,5:2:21),0*v(1:10:end,1:10:end,5:2:21),w(1:10:end,1:10:end,5:2:21),.3,'k');
quiver3(LON(1:3:end-1,1:3:end-1,1:2:21)-360,LAT(1:3:end-1,1:3:end-1,1:2:21),DEPTH(1:3:end-1,1:3:end-1,1:2:21),0*u(1:3:end,1:3:end,1:2:21),0*v(1:3:end,1:3:end,1:2:21),w(1:3:end,1:3:end,1:2:21)*86400,1,'k');

rotate3d on
xlabel('Longitude')
xticks([-72 -66 -60 -54 -48 -42 ])
xticklabels({'72°','66°', '60°', '54°' ,'48°', '42°' })
ylabel('Latitude')
yticks(32:2:44);
yticklabels({'32°','34°', '36°', '38°' ,'40°', '42°','44°' })
zlabel('Depth')

% figure(9)
% %subplot(2,1,1)
% %for d = 2:30    
% %clevs = [-10:1:-1,1:1:10];
% %contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (w(:,:,d))*86400,clevs, 'edgecolor', 'none');
% contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (W_Ro(:,:,d))*86400, 'edgecolor', 'none');
% 
% %colormap(jet(256));
% colormap(redblue)
% hold on
% qScale = 2.5;
% q = quiver(domain.lon(rj)-360, data.lat(ri), u(:,:,d), v(:,:,d),qScale,'k');
% %q.Color = 'black';  
% dNo = num2str(d);
% dVal = num2str(DEPTH(1,1,d));
%     title(['W (not interpolated)  at  ' dVal 'm depth']);
%     xlabel('Longitude');
%     ylabel('Latitiude');
%     hcb = colorbar;
%     set(get(hcb,'Title'),'String',' W [m/day]');
% xticks([-72 -66 -60 -54 -48 -42 ])
% xticklabels({'72°','66°', '60°', '54°' ,'48°', '42°' })
    
    % subplot(2,1,2)
% %for d = 2:30    
% 
% contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (wTemp(:,:,d))*86400, 'edgecolor', 'none');
% %colormap(jet(256));
% colormap(redblue)
% hold on
% %q = quiver(domain.lon(rj)-360, data.lat(ri), u(:,:,d), v(:,:,d),'k');
% %q.Color = 'black';  
% dNo = num2str(d);
% dVal = num2str(DEPTH(1,1,d));
%     title(['W (interpolated)  at  ' dVal 'm depth']);
%     xlabel('Longitude');
%     ylabel('Latitiude');
%     hcb = colorbar;
%     set(get(hcb,'Title'),'String',' W [m/day]');
% 
    %caxis([-1 1]);

% 
% figure(5)
% 
% contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (w(:,:,d)-wTemp(:,:,d))*86400, 'edgecolor', 'none');
% %colormap(jet(256));
% colormap(redblue)
% hold on
% %q = quiver(domain.lon(rj)-360, data.lat(ri), u(:,:,d), v(:,:,d),'k');
% %q.Color = 'black';  
% dNo = num2str(d);
% dVal = num2str(DEPTH(1,1,d));
%     title(['W( not interpolated)  - W (interpolated)  at  ' dVal 'm depth']);
%     xlabel('Longitude');
%     ylabel('Latitiude');
%     hcb = colorbar;
%     set(get(hcb,'Title'),'String',' W [m/day]');
% %caxis([-1 1]);
% 
% 
% figure(12)
% %subplot(2,1,1)
% 
% contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (w(:,:,d))*86400, 'edgecolor', 'none');
% %colormap(jet(256));
% colormap(redblue)
% hold on
% %q = quiver(domain.lon(rj)-360, data.lat(ri), u(:,:,d), v(:,:,d),'k');
% %q.Color = 'black';  
% dNo = num2str(d);
% dVal = num2str(DEPTH(1,1,d));
%     title(['W (constant N) at  ' dVal 'm depth']);
%     xlabel('Longitude');
%     ylabel('Latitiude');
%     hcb = colorbar;
%     set(get(hcb,'Title'),'String',' W [m/day]');
% caxis([-2*10^-16 2*10^-16])
%     subplot(2,1,2)
% %for d = 2:30    
% 
% contourf(mesh.X(:,:,d)-360, mesh.Y(:,:,d), (Qy(:,:,d)), 'edgecolor', 'none');
% %colormap(jet(256));
% colormap(redblue)
% hold on
% %q = quiver(domain.lon(rj)-360, data.lat(ri), u(:,:,d), v(:,:,d),'k');
% %q.Color = 'black';  
% dNo = num2str(d);
% dVal = num2str(DEPTH(1,1,d));
%     title(['Qy at  ' dVal 'm depth']);
%     xlabel('Longitude');
%     ylabel('Latitiude');
%     hcb = colorbar;
%     set(get(hcb,'Title'),'String','  [1/s^3]');
% caxis([-4*10^-12 4*10^-12])
