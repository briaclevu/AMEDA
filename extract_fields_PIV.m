function extract_fields_PIV
%extract_fields.m
%
%  Concatene velocities field from a list of vectors of u and v in m/s
%
%  For a description of the input parameters see param_eddy_tracking.m.
%
%  Nov 2015 Briac Le Vu
%
%=========================

param_eddy_tracking_PIV

global path_in

clear files list result
% get velocities field
[~,result] = system(['dir ',path_in,'*/*txt'])
files = textscan( result(1:end), '%s', 'delimiter', '\n' )
list = files{1}

% initialize mat
x = zeros(350,347);
y = x;
mask = x;
u = nan(350,347,length(list));
v = u;

% boucle sur les fichiers 'list'
for k=1:length(list)

    file = list{k}

    if exist(file,'file')
        
        % lecture de l'entete du fichier j et des données
        fid = fopen(file);
        dataE = textscan(fid,'%s',1,'Delimiter','"');
        dataH = textscan(fid,'%s%s%s%s%s%s',1,'Delimiter',' ');
        dataM = textscan(fid,'%s%s%s%s','Delimiter','\t');
        fclose(fid);

        % fill the matrices (lon,lat,time)
        if k==1
            x = reshape(str2double(strrep(dataM{1}, ',', '.')),size(x));
            y = reshape(str2double(strrep(dataM{2}, ',', '.')),size(x));
        end
        u(:,:,k) = reshape(str2double(strrep(dataM{3}, ',', '.')),size(x));
        v(:,:,k) = reshape(str2double(strrep(dataM{4}, ',', '.')),size(x));
        mask(u(:,:,k).*v(:,:,k)~=0)=1;
    end
end
x = x*1e-6; % from mm to km
y = y*1e-6; % from mm to km

% determine grid spacing xd in km
Dx = diff(x,1,1);
dx = mean(Dx(:)); % in km

Dy = diff(y,1,2);
dy = mean(Dy(:)); % in km

xd = (dx+dy)/2;

% store data in 'mat' file velocities in m/s and grid in km
save([path_in,'EXP19_new'],'x','y','mask','dx','dy','xd','u','v')

