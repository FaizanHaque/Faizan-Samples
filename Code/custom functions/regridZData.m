function [u v T S ] = regridZData(u0, v0,T0,S0,z0,zt,ri,rj,rk)
S0 = permute(S0,[3 1 2]);
T0 = permute(T0,[3 1 2]);
v0 = permute(v0,[3 1 2]);
u0 = permute(u0,[3 1 2]);
S = 0*S0;        
T = 0*T0;        
u = 0*u0;        
v = 0*v0;        



for j = rj
    for i = ri
        S(:,i,j) = interp1(zt,S0(:,i,j),z0);
        T(:,i,j) = interp1(zt,T0(:,i,j),z0);
        u(:,i,j) = interp1(zt,u0(:,i,j),z0);
        v(:,i,j) = interp1(zt,v0(:,i,j),z0);

    end    
end
S = permute(S,[2 3 1]);
T = permute(T,[2 3 1]);
u = permute(u,[2 3 1]);
v = permute(v,[2 3 1]);

