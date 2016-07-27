function field=get_missing_val(lon,lat,field,missvalue,default)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function field=get_missing_val(lon,lat,field,missvalue,ro,default)
%
%  pierrick 2001
%
%  perform an objective analysis to fill
%  the missing points of an horizontal gridded slice
%
%
%  input: 
%    lon      : longitude
%    lat      : latitude
%    field    : input 2D field
%    missvalue: value of the bad points (e.g. -99.999)
%    ro       : oa decorrelation scale
%    default  : default value given if there is only missing data
%
%  output:
%    field    : output 2D field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
  missvalue=NaN;
  default=0;
elseif nargin<5
  default=0;
end

%
%  get a masking matrix and the good data matrix
%
if isnan(missvalue)
  ismask=isnan(field);
else
  ismask=(field==missvalue);
end
isdata=1-ismask;
[M,L]=size(field);

if sum(size(lon))==(length(squeeze(lon))+1)
  [lon,lat]=meshgrid(lon,lat);
end
%
% test if there are any data
%
if (sum(sum(isdata))==0) 
%  disp('no data')
  field=zeros(M,L)+default;
  return
elseif (sum(sum(isdata))<6) 
  default=min(field(isdata==1));
  disp('no enough data to fill missing values')
  disp([' ... using default value:',num2str(default)])
  field=zeros(M,L)+default;
  interp_flag=0;
  return
end
if (sum(sum(ismask))==0) 
%  disp('no mask')
  return
end

%---------------------------------------------------------------
% Extrapolation using nearest values
%--------------------------------------------------------------

field(ismask)=griddata(lon(~ismask),lat(~ismask),field(~ismask),...
                       lon(ismask),lat(ismask),'nearest');
return
