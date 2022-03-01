function data = p_get_data(fname);
%data.u = ncread(fname, 'zvelocity'); %(m/s)
%data.v = ncread(fname, 'mvelocity'); %(m/s)
data.u = ncread(fname, 'ugo'); %(m/s)
data.v = ncread(fname, 'vgo'); %(m/s)
data.salt = double(ncread(fname,'so'));% practical salinity unit(0.001)
data.temp = double(ncread(fname,'to')); %temperature (C*)

data.lat = double(ncread(fname, 'latitude')); 
data.lon = double(ncread(fname, 'longitude')); 
data.depth = double(ncread(fname, 'depth'));
% data.salt = double(ncread(fname,'salinity'));
% data.temp = double(ncread(fname,'temperature'));
data.dlat = data.lat(2)-data.lat(1); % assume constant resolution
data.dlon = data.lon(2)-data.lon(1); % assume constant resolution

data.v = permute(data.v,[2 1 3 4]);
data.u = permute(data.u,[2 1 3 4]);
data.salt = permute(data.salt,[2 1 3 4]);
data.temp = permute(data.temp,[2 1 3 4]);



