function mod_eddy_centers(K,rd,update)
%  mode_eddy_centers(K,rd,update) detect the potential eddy centers  
%  present in the domain, for each time step of the time 
%  series of the 2-D velocity field defined by u and v, using the 
%  LNAM and Okubo-Weiss field calculated in mod_fields.m
%  - K is the abs(LNAM) threshold to delimit contour around each center
%  - rd is used to define smaller areas around each center and find which
%    ones are included (and alone) in a closed contour
%  - update is a flag allowing to update an existing tracking
%   (update is also = number of time steps backward to consider)
%
%  Potentiel centers selected are the max(|LNAM(OW<0)>K|) with at least one
%  closed contour of ssh (or psi by using compute_psi) around them
%
%  For a description of the input parameters see param_eddy_tracking.m.
%
%  Eddy centers are saved in the structure array [path_out,'eddy_centers']:
%
%  - centers(t).step : step when the eddy was detected
%  - centers(t).type(n) : eddy type (-1 => cyclonic; 1 => anticyclonic)
%  - centers(t).lat(n) : eddy center latitude
%  - centers(t).lon(n) : eddy center longitude
%  - centers(t).j(n) : eddy center row index
%  - centers(t).i(n) : eddy center column index
%  
%  - centers0 : max of |LNAM(OW<0)>K|
%  
%  (t is the time index; n is the number of eddies detected at t)
%
%-------------------------
%   Ver. 3.2 Apr.2015 Briac Le Vu
%   Ver. 3.1 2013 LMD
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
global H

% load the computed field in mod_fields
load([path_out,'fields_inter_',runname]);
step = size(u,3);

% to calculate psi extrapole u and v to 0 in the land
u(isnan(u))=0;
v(isnan(v))=0;

% preallocate centers array if doesn't exist
if update && exist([path_out,'eddy_centers_',runname,'.mat'],'file')
    load([path_out,'eddy_centers_',runname]);
    step0=step-update+1;
else
    centers0 = struct('step',{},'type',{},'lon',{},'lat',{}, ...
        'j',{},'i',{});
    centers = centers0;
    step0=1;
end

%---------------------------------------------
disp(['Find potential centers from step ',num2str(step0),' to ',num2str(step)])

% cycle through time steps
for i=step0:step

    disp(['  Search centers step ',num2str(i)])

    % eddy centers for a given step k
    centers0(i).step=i;

    % LNAM n OW criteria define contour including potential centers
    OW = detection_fields(i).LOW; % Okubo Weiss
    LNAM = detection_fields(i).LNAM; % LNAM
    LOW = abs(LNAM);
    LOW(OW>=0 | isnan(OW)) = 0;
    
    % test if the grid is regular or not
    if min(lon(1,:)==lon(end,:)) && min(lat(:,1)==lat(:,end))
        CS=contourc(lon(1,:),lat(:,1),LOW,[K K]);
    else
        figure('visible','off')
        CS=contour(lon,lat,LOW,[K K]);
    end

    % Initialization
    j=1; % first coordinates of the contour scan
    k=1; % first contour

    % scan each LNAM contour
    while j < size(CS,2)
        n = CS(2,j); % number of coordinates for the contour(j)
        xv = CS(1,j+1:j+n); % x values serie for the contour(j) coordinates
        yv = CS(2,j+1:j+n); % y values serie for the contour(j) coordinates
        % validate only bigger contour
        if n > 6
            % make a mask outside the contour
            in = inpolygon(lon,lat,xv,yv);
            maskin=mask;
            maskin(~in)=0;
            maskin(in)=1;
            Lm=LNAM.*maskin;
            if sum(maskin(:))>0
            % L maximum value inside the contour(j)  
                LC = Lm(abs(Lm)==max(abs(Lm(:))));
            % save coordinates of the L maximum inside the contour(j)
                if mask(find(Lm==LC,1))==1
                    centers0(i).type(k)=-sign(LC(1));
                    centers0(i).lon(k)=lon(find(Lm==LC,1));
                    centers0(i).lat(k)=lat(find(Lm==LC,1));
                    [centers0(i).j(k),centers0(i).i(k)]=find(Lm==LC,1);
                % increment the counter
                    k = k + 1; % next contour
                end
            end
        end
    % increment the counter
        j = j + n + 1; % series of coordinates of the next contour 
    end
    
    %----------------------------------------------
    % remove center with no close curve around only one center
    disp('    Remove centers without closed streamlines')

    % all centers coordinates for a given step i
    centers(i).step=i;
    centers_lon=centers0(i).lon;
    centers_lat=centers0(i).lat;

    type1=nan(1,length(centers_lon));
    lon1=nan(1,length(centers_lon));
    lat1=nan(1,length(centers_lon));
    j1=nan(1,length(centers_lon));
    i1=nan(1,length(centers_lon));

    %----------------------------------------------
    % compute each centers in a smaller area

    for ii=1:length(centers_lon)
        
        % fix the indice of the main center
        C_I = centers0(i).i(ii);
        C_J = centers0(i).j(ii);
        % center coordinates
        ll_ci=centers0(i).lon(ii);
        ll_cj=centers0(i).lat(ii);

        % resize coordinate and velocity matrix 
        % (making sure not to go outside the domain)
        lt=lat(max(C_J-rd,1):min(C_J+rd,size(lat,1)), ...
            max(C_I-rd,1):min(C_I+rd,size(lat,2)));
        ln=lon(max(C_J-rd,1):min(C_J+rd,size(lon,1)), ...
            max(C_I-rd,1):min(C_I+rd,size(lon,2)));
        mk=mask(max(C_J-rd,1):min(C_J+rd,size(mask,1)), ...
            max(C_I-rd,1):min(C_I+rd,size(mask,2)));
        vv=v(max(C_J-rd,1):min(C_J+rd,size(v,1)), ...
            max(C_I-rd,1):min(C_I+rd,size(v,2)),i);
        uu=u(max(C_J-rd,1):min(C_J+rd,size(u,1)), ...
            max(C_I-rd,1):min(C_I+rd,size(u,2)),i);
        if type_detection==2 || type_detection==3
            sshh=ssh(max(C_J-rd,1):min(C_J+rd,size(ssh,1)), ...
                max(C_I-rd,1):min(C_I+rd,size(ssh,2)),i);
        end
        
        % indice of the center in the small domain
        [cj,ci] = find(lt==ll_cj & ln==ll_ci);
        
        % indices of all eddy centers in the smaller area
        % (contains at least c_j and c_i)
        cts_j=zeros(length(centers_lat),1);
        cts_i=cts_j;
        for k=1:length(centers_lat)
            try
                [cts_j(k),cts_i(k)]=find(lt==centers_lat(k) & ln==centers_lon(k));
            catch
            end
        end
        % can be empty
        zero_centers=find(cts_j==0 | cts_i==0);
        cts_j(zero_centers)=[];
        cts_i(zero_centers)=[];
        % centers position in the smaller area
        if ~isempty(cts_i)
            ll_ctsi=zeros(length(cts_i),1);
            ll_ctsj=ll_ctsi;
            for k=1:length(cts_i)
                ll_ctsj(k)=lt(cts_j(k),cts_i(k));
                ll_ctsi(k)=ln(cts_j(k),cts_i(k));
            end
        else
            ll_ctsj=[];
            ll_ctsi=[];
        end

        % initialize streamlines to be scanned
        CS1=[];
        CS2=[];

        if type_detection==1 || type_detection==3
        %----------------------------------------------
        % compute the psi field from velocity fields
            psi = compute_psi(ln,lt,mk,uu/100,vv/100,ci,cj);
            %H = double(floor(nanmin(psi(:))*1000):ceil(nanmax(psi(:))*1000))/1000;

            % test if the grid is regular or not
            if min(ln(1,:)==ln(end,:)) && min(lt(:,1)==lt(:,end))
                CS1=contourc(ln(1,:),lt(:,1),psi,H(1:2:end));
            else
                figure('visible','off')
                CS1=contour(ln,lt,psi,H(1:2:end));
            end
        end

        if type_detection==2 || type_detection==3
        %----------------------------------------------
        % compute the psi field from ssh

            psi = squeeze(sshh);

            % test if the grid is regular or not
            if min(ln(1,:)==ln(end,:)) && min(lt(:,1)==lt(:,end))
                CS2=contourc(ln(1,:),lt(:,1),psi,H(1:2:end));
            else
                figure('visible','off')
                CS2=contour(ln,lt,psi,H(1:2:end));
            end
        end

        % concantene all streamlines
        CS=[CS1,CS2];

        %-----------------------------------------------------------
        % rearrange all the contours in C to the structure array 'isolines'
        % sort by maximum latitude. Each element of isolines contains
        % all the vertices of a given contour level of PSI

        % fill the structure 'isolines'
        isolines=struct('x',{},'y',{},'l',{});

        % begin two counters
        k = 1;
        kk = 1;

        while k < size(CS,2)
            npoints = CS(2,k);
            isolines(kk).x = CS(1,k+1:k+npoints); % vertex lon's
            isolines(kk).y = CS(2,k+1:k+npoints); % vertex lat's
            isolines(kk).l = max(CS(2,k+1:k+npoints)); % max lat of a curve
            kk=kk+1;
            k=k+npoints+1;
        end

        % sort the contours according to their maximum latitude; this way the first
        % closed contour across which velocity increases will also be the largest
        % one (it's the one which extend further north).
        [~,order] = sort([isolines(:).l],'ascend');
        isolines=isolines(order);

        %----------------------------------------------
        % scan streamlines and validate centers which are alone as potential centers

        % Initialization
        j=1; % first coordinates of the contour scan

        while j <= length(isolines)

            xdata=isolines(j).x; % vortex lon's
            ydata=isolines(j).y; % vortex lat's

            % conditions to have determine a ptential center
            % 1) a closed contour
            % 2) detected the right eddy center inside the polygon

            if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
                    inpolygon(ll_ci,ll_cj,xdata,ydata)

                % searchs the coordinates of the centers in the contour 
                IN=inpolygon(ll_ctsi,ll_ctsj,xdata,ydata);
                [p,~]=find(IN==1); % index of the coordinates in the contour
                nc=length(p); % number of center in the contour
                
                % only one center include in the streamline
                if nc==1
                    type1(ii)=centers0(i).type(ii);
                    lon1(ii)=ll_ci;
                    lat1(ii)=ll_cj;
                    j1(ii)=C_J;
                    i1(ii)=C_I;
                % the contour contains more than 2 centers
                elseif nc>1
                    j=length(isolines); % stop the scan
                end
            end
            % increment the counter
            j=j+1;
        end

    end % centers ii loop

    %----------------------------------------------
    % select only centers with close contour and no other center inside

    % Initialization
    j=1; % first coordinates of the contour scan

    for k=1:length(type1)
    % remove centers with no close contour
        if ~isnan(type1(k))
            centers(i).type(j)=type1(k);
            centers(i).lon(j)=lon1(k);
            centers(i).lat(j)=lat1(k);
            centers(i).j(j)=j1(k);
            centers(i).i(j)=i1(k);
            j=j+1;
        end
    end
        
end % time i loop

disp(' ')

% save centers in struct array
%----------------------------------------

if update && exist([path_out,'eddy_centers_',runname,'.mat'],'file')
    save([path_out,'eddy_centers_',runname],'centers0','centers','-append')
else
    save([path_out,'eddy_centers_',runname],'centers0','centers')
end

