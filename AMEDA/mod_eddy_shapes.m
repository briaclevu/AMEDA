function mod_eddy_shapes(rd,update)
%  mod_eddy_shapes(rd,update) computes the shapes of
%  the eddies identified by the centers detected with mod_eddy_centers.m
%  and saves them in [path_out,'eddy_shapes'];
%
%  - rd defines the area where the streamfunction (PSI) is firstly computed;
%  - update is a flag allowing to update an existing tracking (update is also 
%    = number of steps backward to consider)
%
%  (For a description of the input parameters see param_eddy_tracking.m)
%
%  Eddy centers are updated with field about single and deouble eddies:
%
%  - centers2(t).lon1,lat1,lon2,lon2(n) : coordinates of the centers
%
%  Eddy shapes are saved in the structure array [path_out,'eddy_shapes']:
%  shapes1 for the shapes2 with one center and shapes2 for double eddies
% 
%  - lonlat : lonlat is a [2xm double] cell array, containing
%    the longitude and latitude position of the m vertices that define the
%    boundaries of the n-th eddy of the t-th step
%  - velmax : velmax is the maximum mean velocity along the contour lonlat
%  - tau : tau is the minimal eddy turnover time inside
%    the eddy bourndaries
%  - deta : deta is the difference between the ssh on the
%    limiting contour and the ssh in the extremum in the middle (max or min)
%  - rmax : rmax is the averaged radius from the eddy centre to
%    the boundary vertices (see 'rayon_moyen' for more details)
%  - aire : aire is the area delimited by the eddy boundaries
%    (see 'rayon_moyen' for more details)
%
%  Another file which contains information on the process to  compute eddy 
%  shapes is saved as [path_out,'warnings_shapes']:
%
%  - warn_shapes(t).no_curve(n): 1 if no closed contour of PSI was found
%    around the eddy center;
%  - warn_shapes(t).fac(n): number of times the area where the shape is
%    computed was enlarged;
%  - warn_shapes(t).calcul_curve(n): 1 if the closed contour as been obtain
%    through the computatio of psi field thanks to 'compute_psi' in the
%    case ssh do not allowed to close contour around the center
%  - warn_shapes(t).large_curve(n): 1 if no closed contour around the
%    center with increasing velocity across. Shape is defined by just
%    the larges closed contour;
%  - warn_shapes(t).too_large2(n): 1 if double eddy shapes is too large
%
%  (n is the total number of eddy for a given step; t the total number of
%  steps)
%
%  The function returns also a log file [path_out,'log_eddy_shapes.txt']
%  which contains additional information on the process to compute eddy
%  shapes.
%  (NOTE: if the file already exist, the new log file will be append to it)
%  
%  'eddy_dim' is used to compute eddy shapes. Check the documentation in 
%  eddy_dim.m for further details.
%
%-------------------------
%   Ver. 3.2 Apr.2015 Briac Le Vu
%   Ver. 3.1 2014 LMD
%   Ver. 2.1 Oct.2012
%   Ver. 2.0 Jan.2012
%   Ver. 1.3 Apr.2011
%   Ver. 1.2 May.2010
%   Ver. 1.1 Dec.2009
%   Authors: Francesco Nencioli, francesco.nencioli@univ-amu.fr
%            Charles Dong, cdong@atmos.ucla.edu
%-------------------------
%
% Copyright (C) 2009-2012 Francesco Nencioli and Charles Dong
%
% This file is part of the Vector-Geometry Eddy Detection Algorithm.
%
% The Vector-Geometry Eddy Detection Algorithm is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
% 
% The Vector-Geometry Eddy Detection Algorithm is distributed in the 
% hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Vector-Geometry Eddy Detection Algorithm.  
% If not, see <http://www.gnu.org/licenses/>.
%
%=========================

global path_out
global runname
global type_detection
global streamlines
global daystreamfunction

% load the computed field in mod_fields
load([path_out,'fields_',runname]);
% begin the log file
diary([path_out,'log_eddy_shapes_',runname,'.txt']);
% load eddy centers
load([path_out,'eddy_centers_',runname]);

step=size(u,3);

%----------------------------------------------
% preallocate shape and warning array
global extended_diags
if update && exist([path_out,'eddy_shapes_',runname,'.mat'],'file')
    load([path_out,'eddy_shapes_',runname])
    load([path_out,'warnings_shapes_',runname])
    step0=step-update+1;
    centers2=centers2(1:step0-1);
    shapes1=shapes1(1:step0-1);
    shapes2=shapes2(1:step0-1);
    warn_shapes=warn_shapes(1:step0-1);
    warn_shapes2=warn_shapes2(1:step0-1);
else
% preallocate shape and warning array
    if extended_diags==0 || extended_diags==2
        shapes1 = struct('lonlat',{},'velmax',{},'tau',{},'deta',{},'rmax',{},'aire',{});
        shapes2 = shapes1;
    elseif extended_diags==1
        shapes1 = struct('lonlat',{},'velmax',{},'tau',{},'deta',{},'rmax',{},'aire',{},...
            'xbary',{},'ybary',{},'z1',{},'z2',{},'a',{},'b',{},'alpha',{},'ellip',{},...
            'ke',{},'vort',{},'vortM',{},'OW',{},'LNAM',{});
        shapes2 = shapes1;
    else
        display('Wrong choice of extended_diags option')
        stop
    end
    centers2 = struct('step',{},'type',{},'lon1',{},'lat1',{},'lon2',{},'lat2',{},'dc',{},'ds',{});
    % shapes struct for the possible second shape with 2 centers
    warn_shapes = struct('no_curve',{},'fac',{},'rd',{},'calcul_curve',{},...
	'large_curve1',{},'large_curve2',{},'too_large2',{});
    warn_shapes2 = warn_shapes;
    if streamlines
        profil2 = struct('step',{},'nc',{},'eta',{},'rmoy',{},'vel',{},'tau',{});
        fname = fieldnames(profil2);
    end
    step0=1;
end
 
disp(['Detrmine contour shapes on ',num2str(rd),'X',num2str(rd),' grid for ',runname])

%----------------------------------------------
% Compute eddy shape

% loop through all steps of the time series
for i=step0:step

    disp(['  Step ',num2str(i),' %-------------'])

    uu = squeeze(u(:,:,i));
    vv = squeeze(v(:,:,i));
    if type_detection==2 || type_detection==3
        sshh = squeeze(ssh(:,:,i));
    else
        sshh = [];
    end
    
    if streamlines
        profil2(i).step=centers(i).step;
    end
    centers2(i).step=centers(i).step;
    shapes1(i).velmax=[];
    shapes2(i).velmax=[];
    warn_shapes(i).no_curve=1;
    
    %----------------------------------------------
    % loop through all centers detected for a given step
    for ii=1:length(centers(i).type)

        disp([' === Center ',num2str(ii),' ==='])

        % initialization
        box=1; % flag that indicates that no maximum velocity was found yet
        fac=0; % increase factor for the area where PSI is computed
        tmp_velmax=zeros(1,2); % velocity to be tested every eddy_dim computation
        tmp_allines=0; % streamlines features to be tested every eddy_dim computation

        while box
            % factor which increases the area where psi is computed
            % if eddy dimensions are too big (box=1)
            fac=fac+1; % this determine larger area

            %----------------------------------------------
            % lonlat is the computed shape;
            % the others are flags output by eddy_dim, and saved in 'warnings_shapes';
            [CD,lonlat,allines,velmax,tau,deta,large,warn,box,calcul]=...
            	eddy_dim(uu,vv,sshh,mask,lon,lat,centers(i),fac,rd,ii);

            %----------------------------------------------
            % flags exploitation
            if warn
                disp('   No significant streamlines closed around the center')

            % No maximum velocity find but only an increasing; test a bigger area
            elseif box
                % if eddy velocity is still increasing then
                % temporary save lonlat
                %if velmax(1) > tmp_velmax(1)+0.0001 || velmax(2) > tmp_velmax(2)+0.0001
                if size(allines,1) > size(tmp_allines,1) && ...
                        (velmax(1) > tmp_velmax(1) || velmax(2) > tmp_velmax(2))
                    disp(['   Big eddy: going to fac = ',num2str(fac+1)])
                    % temporary save flags
                    tmp_allines=allines;
                    tmp_velmax=velmax;
                % if no closed curve more intense in the larger area then
                % final eddy shape is the one computed in the smaller area
                else
                    disp(['   No closed or largest curve at fac ',num2str(fac),...
                        ' back to the largest curve at fac ',num2str(fac-1)])
                    % stop enlarging the area
                    box=0;
                    fac=fac-1;
                    % go back to saved eddy_shape
                    allines=tmp_allines;
                    velmax=tmp_velmax;
                end
            end
        end % end while loop

        %----------------------------------------------
        % which kind of eddy found
        if ~warn
            if large(2) == 0
                disp('   Eddy with 2 centers')
            elseif large(1) == 0
                disp('   Eddy with 1 center')
            elseif large(2) == 1
                disp('   Largest with 2 centers')
            else
                disp('   Largest with 1 center')
            end
        end

        %----------------------------------------------------------
        % save dim in a struct array shapes1 the single eddies
        if ~isempty(lonlat{1})
            [rmax,aire,~,~] = rayon_moyen(lonlat{1});
            centers2(i).type(ii)=centers(i).type(ii);
            centers2(i).lon1(ii)=centers(i).lon(ii);
            centers2(i).lat1(ii)=centers(i).lat(ii);
            shapes1(i).lonlat(ii)=lonlat(1);
            shapes1(i).velmax(ii)=velmax(1);
            shapes1(i).tau(ii)=tau(1);
            shapes1(i).deta(ii)=deta(1);
            shapes1(i).rmax(ii)=rmax(1);
            shapes1(i).aire(ii)=aire;
            if extended_diags==1
                [xbary,ybary,z1,z2,a,b,alpha,~] = calcul_barycentre(lonlat{1});
                shapes1(i).xbary(ii)=xbary;
                shapes1(i).ybary(ii)=ybary;
                shapes1(i).z1(ii)=z1;
                shapes1(i).z2(ii)=z2;
                shapes1(i).a(ii)=a;
                shapes1(i).b(ii)=b;
                shapes1(i).alpha(ii)=alpha;
                if a>=b && a~=0
                    shapes1(i).ellip(ii)=1-b/a;
                elseif a<b && b~=0
                    shapes1(i).ellip(ii)=1-a/b;
                else
                    shapes1(i).ellip(ii)=NaN;
                end
                in_eddy=inpolygon(lon,lat,lonlat{1}(1,:),lonlat{1}(2,:));
                ke=detection_fields(i).ke(in_eddy~=0);
                vort=detection_fields(i).vort(in_eddy~=0);
                OW=detection_fields(i).OW(in_eddy~=0);
                LNAM=detection_fields(i).LNAM(in_eddy~=0);
                shapes1(i).ke(ii)=nansum(ke(:));
                shapes1(i).vort(ii)=nanmean(vort(:));
                if nanmean(vort(:)) > 0
                    shapes1(i).vortM(ii)=max(vort);
                elseif nanmean(vort(:)) < 0
                    shapes1(i).vortM(ii)=min(vort);
                else
                    display('vort==0?')
                end
                shapes1(i).OW(ii)=nanmean(OW(:));
                shapes1(i).LNAM(ii)=nanmean(LNAM(:));
            end
            %----------------------------------------------------------
            % save the streamlines features at i=daystreamfunction
            if streamlines
                for j=daystreamfunction
                    if j==i
                        name=fieldnames(profil2);
                        for n=2:length(name)
                            eval(['profil2(i).',name{n},'{ii}=allines(:,n-1)'';'])
                        end
                    end
                end
            end
        else
            name=fieldnames(centers2);
            for n=2:length(name)
                eval(['centers2(i).',name{n},'(ii)=NaN;'])
            end
            shapes1(i).lonlat(ii)={NaN};
            name=fieldnames(shapes1);
            for n=2:length(name)
                eval(['shapes1(i).',name{n},'(ii)=NaN;'])
            end
            if streamlines
                for j=daystreamfunction
                    if j==i
                        name=fieldnames(profil2);
                        for n=2:length(name)
                            eval(['profil2(i).',name{n},'(ii)={NaN};'])
                        end
                    end
                end
            end
        end

        %----------------------------------------------------------
        % save dim in a struct array shapes2 double eddies
        if ~isempty(lonlat{2})
            [rmax,aire,~,~] = rayon_moyen(lonlat{2});
            if CD(1,1)==centers(i).lon(ii) && CD(2,1)==centers(i).lat(ii)
                centers2(i).lon2(ii)=CD(1,2);
                centers2(i).lat2(ii)=CD(2,2);
            elseif CD(1,2)==centers(i).lon(ii) && CD(2,2)==centers(i).lat(ii)
                centers2(i).lon2(ii)=CD(1,1);
                centers2(i).lat2(ii)=CD(2,1);
            end
            centers2(i).dc(ii)=sw_dist(CD(2,:),CD(1,:),'km');
            centers2(i).ds(ii)=NaN;
            shapes2(i).lonlat(ii)=lonlat(2);
            shapes2(i).velmax(ii)=velmax(2);
            shapes2(i).tau(ii)=tau(2);
            shapes2(i).deta(ii)=deta(2);
            shapes2(i).rmax(ii)=rmax(1);
            shapes2(i).aire(ii)=aire;
            if extended_diags==1
                [xbary,ybary,z1,z2,a,b,alpha,~] = calcul_barycentre(lonlat{2});
                shapes2(i).xbary(ii)=xbary;
                shapes2(i).ybary(ii)=ybary;
                shapes2(i).z1(ii)=z1;
                shapes2(i).z2(ii)=z2;
                shapes2(i).a(ii)=a;
                shapes2(i).b(ii)=b;
                shapes2(i).alpha(ii)=alpha;
                if a>=b && a~=0
                    shapes2(i).ellip(ii)=1-b/a;
                elseif a<b && b~=0
                    shapes2(i).ellip(ii)=1-a/b;
                else
                    shapes2(i).ellip(ii)=NaN;
                end
                in_eddy=inpolygon(lon,lat,lonlat{2}(1,:),lonlat{2}(2,:));
                ke=detection_fields(i).ke(in_eddy~=0);
                vort=detection_fields(i).vort(in_eddy~=0);
                OW=detection_fields(i).OW(in_eddy~=0);
                LNAM=detection_fields(i).LNAM(in_eddy~=0);
                shapes2(i).ke(ii)=nansum(ke(:));
                shapes2(i).vort(ii)=nanmean(vort(:));
                if nanmean(vort(:)) > 0
                    shapes2(i).vortM(ii)=max(vort);
                elseif nanmean(vort(:)) < 0
                    shapes2(i).vortM(ii)=min(vort);
                else
                    display('vort==0?')
                end
                shapes2(i).OW(ii)=nanmean(OW(:));
                shapes2(i).LNAM(ii)=nanmean(LNAM(:));
            end
        else
            name=fieldnames(centers2);
            for n=5:length(name)
                eval(['centers2(i).',name{n},'(ii)=NaN;'])
            end
            shapes2(i).lonlat(ii)={NaN};
            name=fieldnames(shapes2);
            for n=2:length(name)
                eval(['shapes2(i).',name{n},'(ii)=NaN;'])
            end
        end

        %----------------------------------------------------------
        % warnings from shape computation
        warn_shapes(i).no_curve(ii)=warn;
        warn_shapes(i).fac(ii)=fac;
        warn_shapes(i).rd(ii)=rd;
        warn_shapes(i).calcul_curve(ii)=calcul;
        warn_shapes(i).large_curve1(ii)=large(1);
        warn_shapes(i).large_curve2(ii)=large(2);
        warn_shapes(i).too_large2(ii)=0;

    end % ii=last eddy

    if isempty(centers2(i).type)
        %----------------------------------------------------------
        % Initialize warn_shape for eddies with shapes1
        warn_shapes(i).no_curve = [];
        warn_shapes2(i) = warn_shapes(i);
    else
        %----------------------------------------------------------
        % Initialize warn_shape for eddies with shapes1
        warn_shapes2(i) = warn_shapes(i);
        
        %----------------------------------------------------------
        % remove the shapes and centers where no closed contour was found
        replace=find(isnan(shapes1(i).velmax));
        name=fieldnames(centers2);
        for n=2:length(name)
            eval(['centers2(i).',name{n},'(replace)=[];'])
        end
        name=fieldnames(shapes1);
        for n=1:length(name)
            eval(['shapes1(i).',name{n},'(replace)=[];'])
        end
        name=fieldnames(shapes2);
        for n=1:length(name)
            eval(['shapes2(i).',name{n},'(replace)=[];'])
        end
        name=fieldnames(warn_shapes2);
        for n=1:length(name)
            eval(['warn_shapes2(i).',name{n},'(replace)=[];'])
        end
        if streamlines
            for j=daystreamfunction
                if j==i
                    name=fieldnames(profil2);
                    for n=2:length(name)
                        eval(['profil2(i).',name{n},'(replace)=[];'])
                    end
                end
            end
        end

        %----------------------------------------------------------
        % remove shapes2 too big:
        % if 2 shapes1 exist
        %   remove shapes 2 if double eddy too weak
        %	calcul ds = distance between 2 shapes1 (rmax1 & rmax2)
        %   remove shapes2 with ds > 3/2 (rmax1 + rmax2)
        %   if 2 shapes2 exist remove the calculated and the weakest
        % else
        %	remove shapes2 with dc > 5 * rmax
        % end
        disp(' ')
        for ii=1:length(shapes2(i).velmax)
            % test shapes2(ii) exists
            if ~isnan(shapes2(i).velmax(ii))
                mv=0;
                ll1=shapes1(i).lonlat{ii}; % shapes1(ii)
                % find indice of the other center
                ind=find(centers2(i).lon1==centers2(i).lon2(ii) &...
                centers2(i).lat1==centers2(i).lat2(ii));
                % test shapes1(ind) exists
                if ~isempty(ind)
                    ll2=shapes1(i).lonlat{ind}; % shapes1(ind)
                    % calcul distance betxeen shapes1(ii) and shapes1(ind)
                    centers2(i).ds(ii)=min_dist_shapes(ll1,ll2);
                    % remove shapes2(ii) when shapes1(ind) gots higher velocity
                    if shapes2(i).velmax(ii) < shapes1(i).velmax(ind)
                        disp(['   Double eddy ',num2str(ii),' too weak around second shape !!!'])
                        mv=1;
                    % remove shapes2(ii) when shapes1 far one from the other
                    elseif centers2(i).ds(ii) > 3/2*(shapes1(i).rmax(ii)+shapes1(i).rmax(ind))
                        disp(['   Double eddy ',num2str(ii),' too large space between 2 shapes !!!'])
                        mv=1;
                    end
                    % test shapes2(ind) exist
                    if ~isnan(shapes2(i).velmax(ind))
                        % test shapes2 both calculated or both not calculated
                        if warn_shapes2(i).calcul_curve(ii)==warn_shapes2(i).calcul_curve(ind)
                            % remove shapes2 if weakest
                            if shapes2(i).velmax(ii) < shapes2(i).velmax(ind)
                                mv=1;
                            elseif shapes2(i).velmax(ii) == shapes2(i).velmax(ind)
                                if abs(shapes2(i).deta(ii)) < abs(shapes2(i).deta(ind))
                                    mv=1;
                                end
                            end
                            disp(['   Double eddy ',num2str(ii),' weaker than eddy ',num2str(ind),' !!!'])
                        else
                            % remove shapes2(ii) if calculated
                            if warn_shapes2(i).calcul_curve(ii)==1
                                mv=1;
                            end
                            disp(['   Calculated double eddy ',num2str(ii),' removed  !!!'])
                        end
                    end
                % if shapes1(ind) doesn't exist remove shapes2(ii) with centers
                % very far from each ones and record small shapes2 as shapes1
                else
                    if centers2(i).dc(ii) > 5*shapes1(i).rmax(ii)
                        disp(['   Double eddy ',num2str(ii),' too large distance between 2 centers !!!'])
                        mv=1;
                    end
                    if shapes2(i).aire(ii) < 2*shapes1(i).aire(ii)
                        disp(['   Small double eddies ',num2str(ii),' replace by single eddy !!!'])
                        name=fieldnames(centers2);
                        for n=5:length(name)
                            eval(['centers2(i).',name{n},'(ii)=NaN;'])
                        end
                        name=fieldnames(shapes1);
                        for n=1:length(name)
                            eval(['shapes1(i).',name{n},'(ii)=shapes2(i).',name{n},'(ii);'])
                        end
                        warn_shapes2(i).large_curve1(ii)=warn_shapes2(i).large_curve2(ii);
                        mv=1;
                    end
                end
                if mv
                    shapes2(i).lonlat(ii)={NaN};
                    name=fieldnames(shapes2);
                    for n=2:length(name)
                        eval(['shapes2(i).',name{n},'(ii)=NaN;'])
                    end
                    warn_shapes2(i).too_large2(ii)=1;
                end
            end
        end

        disp(' ')
        
    end

end % i=step

%----------------------------------------
% save warnings, shapes and their centers in structure array
save([path_out,'eddy_centers_',runname],'centers2','-append')
if streamlines
    save([path_out,'eddy_shapes_',runname],'shapes1','shapes2','profil2')
else
    save([path_out,'eddy_shapes_',runname],'shapes1','shapes2')
end
save([path_out,'warnings_shapes_',runname],'warn_shapes','warn_shapes2')

% close log file
diary off

