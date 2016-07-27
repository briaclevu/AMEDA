function mod_eddy_centers(K,update)
%mod_eddy_centers(K {,update} )
%
%  Detect the potential eddy centers present in the domain,
%  for each time step of the time serie of the 2-D velocity fields
%  defined by u and v, using the LNAM and Local Okubo-Weiss fields 
%  calculated by mod_fields.m and stored as {detection_fields(t)} in
%  [path_out,'fields_inter',runname]:
%
%  - K is the abs(LNAM(LOW<0)) threshold to delimit the contour of
%       the potential eddy centers (one per contour)
%  - update is a flag allowing to update an existing detection:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the mod_field serie
%
%  For a description of the input parameters see param_eddy_tracking.m.

%
%  Potentiel centers selected are the max(|LNAM(LOW<0)>K|) with at least
%  one closed contour of ssh (or psi by using compute_psi) around them.
%
%  Potential eddy centers are saved/updated as the structure array
%  {center(t)} in [path_out,'eddy_centers_',runname'] with followings
%  fields (type and coordinates):
%  - centers(t).step : step when the eddy was detected
%  - centers(t).type(n) : eddy type (1 => cyclonic; -1 => anticyclonic)
%  - centers(t).x(n) : eddy center x coordinate
%  - centers(t).y(n) : eddy center y coordinate
%  - centers(t).j(n) : eddy center row index
%  - centers(t).i(n) : eddy center column index
%
%
%  Max of |LNAM(LOW<0)>K| with the same fields as {centers} are
%  also saved in {centers0}
%  
%  (t is the time step index; n is the indice of eddy detected at t)
%
%-------------------------
%   Ver. 3.2 Apr 2015 Briac Le Vu
%   Ver. 3.1 2013 LMD from Nencioli et al. routines
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
global grid_ll
global H

% No update by default
if nargin==2
    update = 0;
end

% load the computed field in mod_fields
load([path_out,'fields_inter',runname]);
stepF = size(u,3);

% to calculate psi extrapole u and v to 0 in the land
u(isnan(u)) = 0;
v(isnan(v)) = 0;

% preallocate centers array if doesn't exist
if update && exist([path_out,'eddy_centers',runname,'.mat'],'file')
    
    load([path_out,'eddy_centers',runname]);
    
    step0 = stepF - update+1;
    
    if centers0(step0-1).step ~= step0-1
        display('Gap with the last recorded step!')
        return
    end
    
    centers0 = centers0(1:step0-1);
    centers  = centers(1:step0-1);
    
else
    centers0 = struct('step',{},'type',{},'x',{},'y',{},'i',{},'j',{});
    centers = centers0;
    
    step0 = 1;
    
end

%---------------------------------------------
disp(['Find potential centers from step ',num2str(step0),' to ',num2str(stepF)])

% cycle through time steps
for i=step0:stepF

    disp(['  Search centers step ',num2str(i)])

    % eddy centers for a given step k
    centers0(i).step = i;

    % LNAM n OW criteria define contour including potential centers
    OW = detection_fields(i).LOW; % Okubo Weiss
    LNAM = detection_fields(i).LNAM; % LNAM
    LOW = abs(LNAM);
    LOW(OW>=0 | isnan(OW)) = 0;
    
    % test if the grid is regular or not
    if min(x(1,:)==x(end,:)) && min(y(:,1)==y(:,end))
        CS = contourc(x(1,:),y(:,1),LOW,[K K]);
    else
        figure('visible','off')
        CS = contour(x,y,LOW,[K K]);
    end

    % Initialization
    j = 1; % first coordinates of the contour scan
    k = 1; % first contour

    % scan each LNAM contour
    while j < size(CS,2)
        
        n = CS(2,j); % number of coordinates for the contour(j)
        xv = CS(1,j+1:j+n); % x values serie for the contour(j) coordinates
        yv = CS(2,j+1:j+n); % y values serie for the contour(j) coordinates
        
        % validate only bigger contour
        if n >= 4
            
            % make a mask outside the contour
            in = inpolygon(x,y,xv,yv);
            maskin = mask;
            maskin(~in) = 0;
            maskin(in)  = 1;
            Lm = LNAM.*maskin;
            
            % L maximum value inside the contour(j)  
            if sum(maskin(:))>0
            
                LC = Lm(abs(Lm)==max(abs(Lm(:))));
                
                % save coordinates of the L maximum inside the contour(j)
                if mask(find(Lm==LC,1))==1
                    
                    centers0(i).type(k) = sign(LC(1));
                    centers0(i).x(k)    = x(find(Lm==LC,1));
                    centers0(i).y(k)    = y(find(Lm==LC,1));
                    [centers0(i).j(k),centers0(i).i(k)] = find(Lm==LC,1);
                    
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
    centers(i).step = i;
    centers_type = centers0(i).type;
    centers_x    = centers0(i).x;
    centers_y    = centers0(i).y;
    centers_i    = centers0(i).i;
    centers_j    = centers0(i).j;

    type1 = nan(1,length(centers_x));
    x1    = type1;
    y1    = type1;
    j1    = type1;
    i1    = type1;

    %----------------------------------------------
    % compute each centers in a smaller area

    for ii=1:length(centers_x)
        
        % fix the indice of the main center
        C_I = centers_i(ii);
        C_J = centers_j(ii);
        % center coordinates
        xy_ci = centers_x(ii);
        xy_cj = centers_y(ii);

        % resize coordinate and velocity matrix 
        % (making sure not to go outside the domain)
        xx = x(max(C_J-round(bx(C_J,C_I)/2),1):min(C_J+round(bx(C_J,C_I)/2),size(x,1)), ...
            max(C_I-round(bx(C_J,C_I)/2),1):min(C_I+round(bx(C_J,C_I)/2),size(x,2)));
        yy = y(max(C_J-round(bx(C_J,C_I)/2),1):min(C_J+round(bx(C_J,C_I)/2),size(y,1)), ...
            max(C_I-round(bx(C_J,C_I)/2),1):min(C_I+round(bx(C_J,C_I)/2),size(y,2)));
        mk = mask(max(C_J-round(bx(C_J,C_I)/2),1):min(C_J+round(bx(C_J,C_I)/2),size(mask,1)), ...
            max(C_I-round(bx(C_J,C_I)/2),1):min(C_I+round(bx(C_J,C_I)/2),size(mask,2)));
        vv = v(max(C_J-round(bx(C_J,C_I)/2),1):min(C_J+round(bx(C_J,C_I)/2),size(v,1)), ...
            max(C_I-round(bx(C_J,C_I)/2),1):min(C_I+round(bx(C_J,C_I)/2),size(v,2)),i);
        uu = u(max(C_J-round(bx(C_J,C_I)/2),1):min(C_J+round(bx(C_J,C_I)/2),size(u,1)), ...
            max(C_I-round(bx(C_J,C_I)/2),1):min(C_I+round(bx(C_J,C_I)/2),size(u,2)),i);
        
        if type_detection==2 || type_detection==3
            sshh = ssh(max(C_J-round(bx(C_J,C_I)/2),1):min(C_J+round(bx(C_J,C_I)/2),size(ssh,1)), ...
                max(C_I-round(bx(C_J,C_I)/2),1):min(C_I+round(bx(C_J,C_I)/2),size(ssh,2)),i);
        end
        
        % indice of the center in the small domain
        [cj,ci] = find(yy==xy_cj & xx==xy_ci);
        
        % indices of all eddy centers in the smaller area
        % (contains at least c_j and c_i)
        cts_j = zeros(length(centers_y),1);
        cts_i = cts_j;
        for k=1:length(centers_y)
            try
                [cts_j(k),cts_i(k)] = find(yy==centers_y(k) & xx==centers_x(k));
            catch
            end
        end
        
        % can be empty
        zero_centers = find(cts_j==0 | cts_i==0);
        cts_j(zero_centers) = [];
        cts_i(zero_centers) = [];
        
        % centers position in the smaller area
        if ~isempty(cts_i)
            xy_ctsi = zeros(length(cts_i),1);
            xy_ctsj = xy_ctsi;
            for k=1:length(cts_i)
                xy_ctsj(k) = yy(cts_j(k),cts_i(k));
                xy_ctsi(k) = xx(cts_j(k),cts_i(k));
            end
        else
            xy_ctsj = [];
            xy_ctsi = [];
        end

        % initialize streamlines to be scanned
        CS1 = [];
        CS2 = [];

        %----------------------------------------------
        % compute the psi field from velocity fields
        if type_detection==1 || type_detection==3
        
            psi = compute_psi(xx,yy,mk,uu/100,vv/100,ci,cj,grid_ll);

            % test if the grid is regular or not
            if min(xx(1,:)==xx(end,:)) && min(yy(:,1)==yy(:,end))
                CS1 = contourc(xx(1,:),yy(:,1),psi,H/2);
            else
                figure('visible','off')
                CS1 = contour(xx,yy,psi,H/2);
            end
        end

        %----------------------------------------------
        % compute the psi field from ssh
        if type_detection==2 || type_detection==3

            psi = squeeze(sshh);

            % test if the grid is regular or not
            if min(xx(1,:)==xx(end,:)) && min(yy(:,1)==yy(:,end))
                CS2 = contourc(xx(1,:),yy(:,1),psi,H/2);
            else
                figure('visible','off')
                CS2 = contour(xx,yy,psi,H/2);
            end
        end

        %----------------------------------------------
        % concantene all streamlines
        CS = [CS1,CS2];

        %-----------------------------------------------------------
        % rearrange all the contours in C to the structure array 'isolines'
        % sort by maximum y coord. Each element of isolines contains
        % all the vertices of a given contour level of PSI

        % fill the structure 'isolines'
        isolines = struct('x',{},'y',{},'l',{});

        % begin two counters
        k = 1;
        kk = 1;

        while k < size(CS,2)
            npoints = CS(2,k);
            lvl(kk) = CS(1,k);
            isolines(kk).x = CS(1,k+1:k+npoints); % vertex x's
            isolines(kk).y = CS(2,k+1:k+npoints); % vertex y's
            isolines(kk).l = max(CS(2,k+1:k+npoints)); % max y of a curve
            kk = kk + 1;
            k = k + npoints + 1;
        end

        % sort the contours according to their maximum y coord; this way the first
        % closed contour across which velocity increases will also be the largest
        % one (it's the one which extend further north).
        [~,order] = sort([isolines(:).l],'ascend');
        isolines = isolines(order);
        % ! Debug ! Test the contour value scanned
        %display([min(diff(lvl)) mean(diff(lvl)) max(diff(lvl))]) % !Debug!

        %----------------------------------------------
        % scan streamlines and validate centers which are alone as potential centers

        % Initialization
        j = 1; % first coordinates of the contour scan

        while j <= length(isolines)

            xdata = isolines(j).x; % vortex x's
            ydata = isolines(j).y; % vortex y's

            % conditions to have determine a ptential center
            % 1) a closed contour
            % 2) detected the right eddy center inside the polygon

            if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
                    inpolygon(xy_ci,xy_cj,xdata,ydata) && length(xdata)>=4

                % searchs the coordinates of the centers in the contour 
                IN = inpolygon(xy_ctsi,xy_ctsj,xdata,ydata);
                [p,~] = find(IN==1); % index of the coordinates in the contour
                nc = length(p); % number of center in the contour
                
                % only one center include in the streamline
                if nc==1
                    
                    type1(ii) = centers_type(ii);
                    x1(ii)    = xy_ci;
                    y1(ii)    = xy_cj;
                    j1(ii)    = C_J;
                    i1(ii)    = C_I;
                    
                % the contour contains more than 1 centers
                elseif nc>1
                    
                    j = length(isolines); % stop the scan
                    
                end
            end
            
            % increment the counter
            j = j + 1;
        end

    end % centers ii loop

    %----------------------------------------------
    % select only centers with close contour and no other center inside

    % Initialisation
    j = 1; % first coordinates of the contour scan

    for k=1:length(type1)
        
        % remove centers with no close contour
        if ~isnan(type1(k))
            
            centers(i).type(j) = type1(k);
            centers(i).x(j)    = x1(k);
            centers(i).y(j)    = y1(k);
            centers(i).j(j)    = j1(k);
            centers(i).i(j)    = i1(k);
            j = j + 1;
            
        end
    end
        
end % time i loop

disp(' ')

%----------------------------------------
% save centers in struct array

if update && exist([path_out,'eddy_centers',runname,'.mat'],'file')
    save([path_out,'eddy_centers',runname],'centers0','centers','-append')
else
    save([path_out,'eddy_centers',runname],'centers0','centers')
end
