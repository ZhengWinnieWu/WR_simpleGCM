%% interpolate MERRA observation temperature to the grid points of model and calculate Teq for next iteration using new combine script and pressure levels in MERRA
clear all; close all; 
% open model file
filename11='$filename'           % e.g. 'ltm.p.00001_00500.atmos_monthly.nc' 

% read model data full field
lon11=ncread(filename11,'lon');
lonb11=ncread(filename11,'lonb');
lat11=ncread(filename11,'lat');
latb11=ncread(filename11,'latb');
level11=ncread(filename11,'level');
time11=ncread(filename11,'time');
t11=ncread(filename11,'temp');

% open file for Teq
filename12='$filename'           % e.g. 'teq_n-1.nc';
 
time12=ncread(filename12,'time');
teq_old=ncread(filename12,'teq');

% if Teq has the same pressure level order as T, skip this step
% Teq from low level to high level the same order as T 
teq=zeros(length(lon11),length(lat11),length(level11),length(time11));
for izteq=1:length(level11)
  teq(:,:,izteq,:)=teq_old(:,:,length(level11)+1-izteq,:);
end

missing_value=-1e+10;
missing_value_obs=1000;              % depend on what Reanalysis data is
teq_n=zeros(length(lon11),length(lat11),length(level11),length(time11)); 

% interpolation Reananlysis data to the model grid, if the reanalysis data is the same as the model grid, skip this step
fn=['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'];
for i=1:12
    
    % open file
    filename1=['/uufs/chpc.utah.edu/common/home/reichler-group3/group/REANALYSIS/ATMOS/MERRA/daily/T/mon_cli/' strtrim(fn(i,:)) '.nc'];
    
    % read data
    lon1=ncread(filename1,'lon');
    lat1=ncread(filename1,'lat');
    lev1=ncread(filename1,'lev');
    time1=ncread(filename1,'time');
    t1=ncread(filename1,'T');
    
    % interpolation
    for ix=1:length(lon1)
        if lon1(ix)<0 
            lon1(ix)=lon1(ix)+360.0; 
        end 
        if lon1(ix)-360.0>-1e-10 
            lon1(ix)=0.0; 
        end 
    end 
    lon1;
    
    logic_missing=t1>500.0;
    t1(logic_missing)=missing_value_obs;
    
    tlon=zeros(length(lon11),length(lat1),length(lev1));
    for izo=1:length(lev1)
        for iyo=1:length(lat1)
            for ixm=1:length(lon11)
                for ixo=1:length(lon1)-1
                    if (lon1(ixo)-lon11(ixm) <= 0 && lon1(ixo+1)-lon11(ixm) >= 0) 
                        if (t1(ixo,iyo,izo)~=missing_value_obs && t1(ixo+1,iyo,izo)~=missing_value_obs) 
                            tlon(ixm,iyo,izo)=t1(ixo,iyo,izo)*(lon1(ixo+1)-lon11(ixm))+t1(ixo+1,iyo,izo)*(lon11(ixm)-lon1(ixo)); 
                            tlon(ixm,iyo,izo)=tlon(ixm,iyo,izo)/(lon1(ixo+1)-lon1(ixo)); 
                        elseif (t1(ixo,iyo,izo)==missing_value_obs && t1(ixo+1,iyo,izo)~=missing_value_obs) 
                            tlon(ixm,iyo,izo)=t1(ixo+1,iyo,izo); 
                        elseif (t1(ixo,iyo,izo)~=missing_value_obs && t1(ixo+1,iyo,izo)==missing_value_obs)
                            tlon(ixm,iyo,izo)=t1(ixo,iyo,izo); 
                        else 
                            tlon(ixm,iyo,izo)=missing_value_obs; 
                        end 
                    end 
                end 
            end 
        end 
    end 
    
    tlat=zeros(length(lon11),length(lat11),length(lev1)); 
    for izo=1:length(lev1)
        for ixm=1:length(lon11)
            for iym=1:length(lat11)
                for iyo=1:length(lat1)-1
                    if (lat1(iyo)-lat11(iym) <= 0 && lat1(iyo+1)-lat11(iym) >= 0) 
                        if (tlon(ixm,iyo,izo)~=missing_value_obs && tlon(ixm,iyo+1,izo)~=missing_value_obs) 
                            tlat(ixm,iym,izo)=tlon(ixm,iyo,izo)*(lat1(iyo+1)-lat11(iym))+tlon(ixm,iyo+1,izo)*(lat11(iym)-lat1(iyo)); 
                            tlat(ixm,iym,izo)=tlat(ixm,iym,izo)/(lat1(iyo+1)-lat1(iyo)); 
                        elseif (tlon(ixm,iyo,izo)==missing_value_obs && tlon(ixm,iyo+1,izo)~=missing_value_obs) 
                            tlat(ixm,iym,izo)=tlon(ixm,iyo+1,izo); 
                        elseif (tlon(ixm,iyo,izo)~=missing_value_obs && tlon(ixm,iyo+1,izo)==missing_value_obs) 
                            tlat(ixm,iym,izo)=tlon(ixm,iyo,izo); 
                        else 
                            tlat(ixm,iym,izo)=missing_value_obs; 
                        end 
                    end 
                end 
            end 
        end 
    end 
    
    t_n=zeros(length(lon11),length(lat11),length(level11));
    for iym=1:length(lat11)
        for ixm=1:length(lon11)
            for izm=1:length(level11)-3
                for izo=1:length(lev1)-1
                    if (lev1(izo)-level11(izm) >= 0 && lev1(izo+1)-level11(izm) <= 7.4506e-08) 
                        if (tlat(ixm,iym,izo)~=missing_value_obs && tlat(ixm,iym,izo+1)~=missing_value_obs) 
                            t_n(ixm,iym,izm)=tlat(ixm,iym,izo)*(level11(izm)-lev1(izo+1))+tlat(ixm,iym,izo+1)*(lev1(izo)-level11(izm)); 
                            t_n(ixm,iym,izm)=t_n(ixm,iym,izm)/(lev1(izo)-lev1(izo+1)); 
                        elseif (tlat(ixm,iym,izo)==missing_value_obs && tlat(ixm,iym,izo+1)~=missing_value_obs) 
                            t_n(ixm,iym,izm)=tlat(ixm,iym,izo+1); 
                        elseif (tlat(ixm,iym,izo)~=missing_value_obs && tlat(ixm,iym,izo+1)==missing_value_obs) 
                            t_n(ixm,iym,izm)=tlat(ixm,iym,izo); 
                        else 
                            t_n(ixm,iym,izm)=missing_value_obs; 
                        end 
                    end 
                end 
            end 
        end 
    end 
   
    % if Reanalysis data has data above 0.1 hPa, skip this step
    % calculation of temperature of above 0.1hPa using linear regression
    for iy=1:length(lat11)
        for ix=1:length(lon11)
            y=squeeze(tlat(ix,iy,40:42));
            x=lev1(40:42);
            p=polyfit(x,y,1);
            for iz=1:3
                tlat_upper(ix,iy,iz)=p(1)*level11(iz+31)+p(2);
            end 
        end 
    end 
    
    t_n(:,:,32:34)=tlat_upper(:,:,1:3);
    
    % new teq by adding up temp, teq and observation temp
    for ix=1:length(lon11)
        for iy=1:length(lat11)
            for iz=1:length(level11)
                if (teq(ix,iy,iz,i)==missing_value) || (t11(ix,iy,iz,i)==missing_value || t_n(ix,iy,iz)==missing_value_obs) 
                    teq_n(ix,iy,iz,i)=missing_value; 
                else
                    teq_n(ix,iy,iz,i)=teq(ix,iy,iz,i)-2.0/3.0*(t11(ix,iy,iz,i)-t_n(ix,iy,iz)); 
                end 
            end 
        end 
    end 
    
    i
end 

% calculation of zonal mean of new teq 
teq_ave=zeros(length(lat11),length(level11),length(time11));
for it=1:12
    for iz=length(level11):-1:1           % length(level11)-1:-1:1
        for iy=1:length(lat11)
            nval=0.0;
            teq_sum=0.0;
            for ix=1:length(lon11)
                if ((teq_n(ix,iy,iz,it)~=missing_value) || ((teq_n(ix,iy,iz,it)<430.0) && (teq_n(ix,iy,iz,it)>100.0)))
                    nval=nval+1;
                    teq_sum=teq_sum+teq_n(ix,iy,iz,it); 
                end 
            end 
            if nval>0.0 
                teq_ave(iy,iz,it)=teq_sum/nval;
            else 
                teq_ave(iy,iz,it)=teq_ave(iy,iz+1,it); 
            end 
        end 
    end 
end 

% this is not necessary, and it is used to replace the missing values in lower levels, can use different ways to treat the missing values
% replace missing values with near pressure level values
for it=1:12
    for iz=1:length(level11)-1
        for iy=1:length(lat11)
            for ix=1:length(lon11)
                nexval=0.0;
                if abs(teq_n(ix,iy,iz,it)-missing_value)<0.0001 
                    mt=0;
                    %if teq_n(ix,iy,iz+1,it)~=missing_value 
                    %    teq_n(ix,iy,iz,it)=teq_n(ix,iy,iz+1,it); 
                    %else 
                    %    teq_n(ix,iy,iz,it)=teq_ave(iy,iz+1,it);
                    %end 
                    for izt=iz+1:length(level11) 
                        if abs(teq_n(ix,iy,izt,it)-missing_value)>0.0001 
                            mt=1; 
                        end 
                        if mt==1 
                            izt; 
                            break 
                        end 
                    end 
                    teq_n(ix,iy,iz,it)=teq_n(ix,iy,izt,it);
                end
                %if (teq_n(ix,iy,iz,it)>380.0 || teq_n(ix,iy,iz,it)<100.0) 
                %if (teq_n(ix,iy,iz,it)<100.0)     
                %  nexval=nexval+1;  
                %  teq_n(ix,iy,iz,it)=teq_ave(iy,iz,it); 
                %end 
            end 
            nexval;
        end 
    end 
    it; 
end 

teq_n(find(teq_n < 100.0))=100.0; 
teq_n(find(teq_n > 400.0))=400.0; 

% this is not necessary
% replace levels above 1hPa with its own Teq of 1hPa
ilev=find(level11==1.0);
for iz=1:length(level11)
  if level11(iz)<1.0 
    %teq_n(:,:,iz,:)=teq_n(:,:,ilev,:);  
    %teq_n(:,:,iz,:)=teq_n(:,:,ilev,:);  
  end 
end 

% re-order the pressure level order, not necessary
% from high level to low level 
teqn_rev=zeros(length(lon11),length(lat11),length(level11),length(time11));
pfull=zeros(length(level11),1);                % new pfull
for i=1:length(level11)
    teqn_rev(:,:,i,:)=teq_n(:,:,length(level11)+1-i,:); 
    pfull(i)=level11(length(level11)+1-i); 
end 

% new phalf corresponding to new pfull
phalf = zeros(length(pfull)+1,1);
phalf(2:end-1) = pfull(2:end) - diff(pfull)/2;
phalf(end) = 1e3; 

% write data into new file
filename='$newfilename'          % e.g. 'teq_n.nc'; 
ncid=netcdf.create(filename,'64BIT_OFFSET');

% define dimensions and variables 
lond_id=netcdf.defDim(ncid,'lon',length(lon11)); 
lonbd_id=netcdf.defDim(ncid,'lonb',length(lonb11));
latd_id=netcdf.defDim(ncid,'lat',length(lat11)); 
latbd_id=netcdf.defDim(ncid,'latb',length(latb11));
pfd_id=netcdf.defDim(ncid,'pfull',length(pfull)); 
phd_id=netcdf.defDim(ncid,'phalf',length(phalf));
td_id=netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED')); 

lon_id=netcdf.defVar(ncid,'lon','double',lond_id);
netcdf.putAtt(ncid,lon_id,'long_name','longitude');
netcdf.putAtt(ncid,lon_id,'units','degrees_E');
netcdf.putAtt(ncid,lon_id,'cartesian_axis','X');
netcdf.putAtt(ncid,lon_id,'edges','lonb');
lonb_id = netcdf.defVar(ncid,'lonb','double',lonbd_id);
netcdf.putAtt(ncid,lonb_id,'long_name','longitude_edges');
netcdf.putAtt(ncid,lonb_id,'units','degrees_E');
netcdf.putAtt(ncid,lonb_id,'cartesian_axis','X');
lat_id=netcdf.defVar(ncid,'lat','double',latd_id);
netcdf.putAtt(ncid,lat_id,'long_name','latitude');
netcdf.putAtt(ncid,lat_id,'units','degrees_N');
netcdf.putAtt(ncid,lat_id,'cartesian_axis','Y');
netcdf.putAtt(ncid,lat_id,'edges','latb');
latb_id=netcdf.defVar(ncid,'latb','double',latbd_id);
netcdf.putAtt(ncid,latb_id,'long_name','latitude_edges');
netcdf.putAtt(ncid,latb_id,'units','degrees_N');
netcdf.putAtt(ncid,latb_id,'cartesian_axis','Y');
pf_id=netcdf.defVar(ncid,'pfull','double',pfd_id);
netcdf.putAtt(ncid,pf_id,'long_name','approx full pressure level');
netcdf.putAtt(ncid,pf_id,'units','hPa');
netcdf.putAtt(ncid,pf_id,'cartesian_axis','Z');
netcdf.putAtt(ncid,pf_id,'positive','down');
netcdf.putAtt(ncid,pf_id,'edges','phalf');
ph_id=netcdf.defVar(ncid,'phalf','double',phd_id);
netcdf.putAtt(ncid,ph_id,'long_name','approx half pressure level');
netcdf.putAtt(ncid,ph_id,'units','hPa');
netcdf.putAtt(ncid,ph_id,'cartesian_axis','Z');
t_id=netcdf.defVar(ncid,'time','double',td_id);
netcdf.putAtt(ncid,t_id,'long_name','time');
netcdf.putAtt(ncid,t_id,'cartesian_axis','T');
netcdf.putAtt(ncid,t_id,'units','days since 0000-00-00 00:00:00');
netcdf.putAtt(ncid,t_id,'cell_methods','time:mean'); 

teq_id=netcdf.defVar(ncid,'temp','double',[lond_id,latd_id,pfd_id,td_id]);
netcdf.putAtt(ncid,teq_id,'long_name','radiative relaxation');
netcdf.putAtt(ncid,teq_id,'units','deg_K');
netcdf.putAtt(ncid,teq_id,'valid_range',[0,500]);
netcdf.putAtt(ncid,teq_id,'missing_value',-1.e+10);  

netcdf.endDef(ncid);
% write fields
netcdf.putVar(ncid,lon_id,lon11);
netcdf.putVar(ncid,lonb_id,lonb11);
netcdf.putVar(ncid,lat_id,lat11);
netcdf.putVar(ncid,latb_id,latb11);
netcdf.putVar(ncid,pf_id,pfull);
netcdf.putVar(ncid,ph_id,phalf);
netcdf.putVar(ncid,t_id,0,length(time12),time12);
netcdf.putVar(ncid,teq_id,[0,0,0,0],[128,64,34,12],teqn_rev);
netcdf.close(ncid)
