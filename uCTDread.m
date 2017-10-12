function [data,mdata]= uCTDread(filename,lagtime)
%
% [data,mdata]= uCTDread(filename)
% [data,mdata]= uCTDread(filename,lagtime)
%
% Reads data from underway CTD files .asc, isolates downcast,
% and transforms them in MATLAB format after computing absolute salinity
% and potential density anomaly.
% It aligns temperature and conductivity measurements before calculating
% salinity by applying a negative time lag to temperature measurements.
%
% The format for the coordinates in the uCTD file should be:
% <degrees> <minutes> <direction>, as in the example below
% *Lat 23 27.197 N
% *Lon 157 51.439 W
% Notice that longitude in the metadata is positive to the West
%
% REQUIREMENTS: 
%   Gibbs-Seawater (GSW) Oceanographic Toolbox, http://www.TEOS-10.org
%   MATLAB R2013b or above (for table format)
%
% OUTPUT:
%   filename: underway CTD file .asc
%   lagtime: lag in seconds for temperature data (default 0.09)
%
% OUTPUT:
%   data: table with columns depth, press, sa, t, ct, sig
%   mdata: metadata structure for the profile
%
% VERSION HISTORY
%   v0.11: adds temperature time-lag correction to compute salinity
% 
% Benedetto Barone - Jul 18, 2017 - Revision 0.11

if nargin < 2
    lagtime = 0.09;
end

f = fopen(filename);

% Read metadata
A = fgetl(f);
while strncmp(A,'*scan#',6) == 0
    A = fgetl(f);
    if strncmp(A,'*Lat',4) == 1
        dlat = textscan(A,'*Lat %f %f %c');
        mdata.lat = dlat{1} + dlat{2}/60;
        if strcmp(dlat{3},'S'),data.lat = -data.lat ,end
    elseif strncmp(A,'*Lon',4) == 1
        dlon = textscan(A,'*Lon %f %f %c');
        mdata.lon = -dlon{1} -dlon{2}/60;
        if strcmp(dlon{3},'E'),data.lon = -data.lon ,end
    elseif strncmp(A,'*Cast',5) == 1
        ddate = textscan(A,'*Cast %*f %s %s %s %s %*f %*f %*s');
        mdata.date = datenum([char(ddate{1}) ' ' char(ddate{2}) ' ' char(ddate{3}) ' ' char(ddate{4})]);
    elseif strncmp(A,'*SerialNumber=7',15) == 1
        dserial = textscan(A,'*SerialNumber= %f');
        mdata.serialn = dserial{1};
    end
end

% Read CTD observations
ddata = textscan(f,'%f %f %f %f','CollectOutput',1);
ddata = ddata{1};
ddata(ddata(:,4)<0,:) = []; % Remove data above sea surface
[~,ind_max] = max(ddata(:,4));
ddata(ind_max:end,:) = [];% Remove upcast
mdata.ind_d = [0; ddata(2:end,4)-ddata(1:end-1,4)];
tm = ddata(:,1)*0.0625;
tm_t = tm - lagtime;
uncorr_t = ddata(:,3); % unaligned temperature
t = interp1(tm_t,ddata(:,3),tm,'linear','extrap'); % align temperature
p = ddata(:,4);
depth = - gsw_z_from_p(p,mdata.lat);
sp = gsw_SP_from_C(ddata(:,2)*10,t,p);
sa = gsw_SA_from_SP(sp,p,mdata.lon,mdata.lat);
uncorr_sp = gsw_SP_from_C(ddata(:,2)*10,uncorr_t,p);
uncorr_sa = gsw_SA_from_SP(uncorr_sp,p,mdata.lon,mdata.lat);
ct = gsw_CT_from_t(sa,t,p);
sig = gsw_sigma0(sa,ct);
data = table(p,depth,t,sp,ct,sa,sig);
data.Properties.VariableUnits = {'db','m','deg C','-','deg C','g kg-1','kg m-3'};
fclose(f);

% Plot results
subplot(1,3,1)
plot(data.t,data.depth,'k-')
set(gca,'ydir','rev')
xlabel('T (deg C)'),ylabel('Depth (m)')
subplot(1,3,2)
plot(uncorr_sa,data.depth,'r-',data.sa,data.depth,'k-')
set(gca,'ydir','rev')
xlim([34 36])
xlabel('S_A (g Kg^{-1})'),ylabel('Depth (m)')
legend('no correction','lag correction')
legend('boxoff')
subplot(1,3,3)
plot(data.sig,data.depth,'k-')
set(gca,'ydir','rev')
xlim([22 28])
xlabel('sigma0 (Kg m^{-3})'),ylabel('Depth (m)')

end