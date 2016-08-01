<<<<<<< HEAD
function mod_eddy_centers(K,update)
%mod_eddy_centers(K {,update} )
=======
function [centers0,centers] = mod_eddy_centers(source,stp,fields,resolution)
%[centers0,centers] = mod_eddy_centers(source,stp,fields {,resolution} )
>>>>>>> ameda_v2
%
%  Detect the potential eddy centers present in the domain,
%  for each time step of the time serie of the 2-D velocity fields
%  defined by u and v, using the LNAM and Local Okubo-Weiss fields 
%  calculated by mod_fields.m and stored as {detection_fields(t)} in
%  [path_out,'fields_inter',runname]:
%
% - 'source' is the type of netcdf data (AVISO, NEMO, ROMS,...) that
%   determine which load_field_'source'is  used in the routine
% - 'fields' is the step 'stp' of the detection_fields computed with
%   mod_eddy_fields.m
%
% This routine will use:
%  - K: the abs(LNAM(LOW<0)) threshold to delimit the contour of
%       the potential eddy centers (one per contour)
<<<<<<< HEAD
%  - update is a flag allowing to update an existing detection:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the mod_field serie
=======
%  - bx: used to define smaller areas around each center and find which
%       ones are included (and alone) in a closed streamline. Must take into
%       account the resolution factor!
>>>>>>> ameda_v2
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
%  - centers(t).i(n) : eddy center column index
%  - centers(t).j(n) : eddy center row index
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
%=========================

%---------------------------------------------
% read fields and initialisation
%---------------------------------------------

disp(['Find potential centers step ',num2str(stp)])

% load key_source and parameters (use mod_eddy_params.m first)
%----------------------------------------------
load('param_eddy_tracking')

% replace parameters by arguments
%----------------------------------------
if nargin==4
    resol = resolution;
end
bx = bx*resol;

% load 2D velocity fields (m/s) for step stp
%----------------------------------------
eval(['[x,y,mask,u,v,ssh] = load_fields_',source,'(stp,resol,deg);'])

% to calculate psi extrapole u and v to 0 in the land
u(isnan(u)) = 0;
v(isnan(v)) = 0;

% initialise centers as structure
%----------------------------------------
centers0 = struct('step',nan,'type',[],'x',[],'y',[],'i',[],'j',[]);
centers = centers0;

%---------------------------------------------
% Max LNAM 'centers0' for a given step k
%---------------------------------------------

centers0.step = stp;

% LNAM n OW criteria define contour including potential centers
%---------------------------------------------
OW = fields.LOW; % Okubo Weiss
LNAM = fields.LNAM; % LNAM
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
    if n >= n_min

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

                centers0.type(k) = sign(LC(1));
                centers0.x(k)    = x(find(Lm==LC,1));
                centers0.y(k)    = y(find(Lm==LC,1));
                [centers0.j(k),centers0.i(k)] = find(Lm==LC,1);

                % increment the counter
                k = k + 1; % next contour
            end
        end
    end

<<<<<<< HEAD
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
=======
    % increment the counter
    j = j + n + 1; % series of coordinates of the next contour 
end

disp(['  ',num2str(k-1),' max LNAM found'])

%----------------------------------------------
% Remove 'centers0' with no close curve around only one max LNAM
%----------------------------------------------

disp('  Remove max LNAM without closed streamlines')

% all max LNAM coordinates for the given step stp
centers.step = stp;
centers_type = centers0.type;
centers_x    = centers0.x;
centers_y    = centers0.y;
centers_i    = centers0.i;
centers_j    = centers0.j;

type1 = nan(1,length(centers_x));
x1    = type1;
y1    = type1;
j1    = type1;
i1    = type1;

% compute each centers in a smaller area
%----------------------------------------------
for ii=1:length(centers_x)

    % fix the indice of the main center
    C_I = centers_i(ii);
    C_J = centers_j(ii);
    % center coordinates
    xy_ci = centers_x(ii);
    xy_cj = centers_y(ii);

    % resize coordinate and velocity matrix 
    % (making sure not to go outside the domain)
    xx = x(max(C_J-round(bx/2),1):min(C_J+round(bx/2),size(x,1)), ...
        max(C_I-round(bx/2),1):min(C_I+round(bx/2),size(x,2)));
    yy = y(max(C_J-round(bx/2),1):min(C_J+round(bx/2),size(y,1)), ...
        max(C_I-round(bx/2),1):min(C_I+round(bx/2),size(y,2)));
    mk = mask(max(C_J-round(bx/2),1):min(C_J+round(bx/2),size(mask,1)), ...
        max(C_I-round(bx/2),1):min(C_I+round(bx/2),size(mask,2)));
    vv = v(max(C_J-round(bx/2),1):min(C_J+round(bx/2),size(v,1)), ...
        max(C_I-round(bx/2),1):min(C_I+round(bx/2),size(v,2)));
    uu = u(max(C_J-round(bx/2),1):min(C_J+round(bx/2),size(u,1)), ...
        max(C_I-round(bx/2),1):min(C_I+round(bx/2),size(u,2)));

    if type_detection>=2
        sshh = ssh(max(C_J-round(bx/2),1):min(C_J+round(bx/2),size(ssh,1)), ...
            max(C_I-round(bx/2),1):min(C_I+round(bx/2),size(ssh,2)));
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
>>>>>>> ameda_v2
        end
    end

    % can be empty
    zero_centers = find(cts_j==0 | cts_i==0);
    cts_i(zero_centers) = [];
    cts_j(zero_centers) = [];

    % centers position in the smaller area
    if ~isempty(cts_i)
        xy_ctsi = zeros(length(cts_i),1);
        xy_ctsj = xy_ctsi;
        for k=1:length(cts_i)
            xy_ctsi(k) = xx(cts_j(k),cts_i(k));
            xy_ctsj(k) = yy(cts_j(k),cts_i(k));
        end
    else
        xy_ctsi = [];
        xy_ctsj = [];
    end

    % initialize streamlines to be scanned
    CS1 = [];
    CS2 = [];

    % compute the psi field from velocity fields
    %----------------------------------------------
    if type_detection==1 || type_detection==3

        psi = compute_psi(xx,yy,mk,uu/100,vv/100,ci,cj,grid_ll);

<<<<<<< HEAD
        %----------------------------------------------
        % concantene all streamlines
        CS = [CS1,CS2];

        % rearrange all the contours in C to the structure array 'isolines'
        %-----------------------------------------------------------
        [isolines,~] = scan_lines(CS);
        
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
=======
        % test if the grid is regular or not
        if min(xx(1,:)==xx(end,:)) && min(yy(:,1)==yy(:,end))
            CS1 = contourc(xx(1,:),yy(:,1),psi,H/2);
        else
            figure('visible','off')
            CS1 = contour(xx,yy,psi,H/2);
        end
    end

    % compute the psi field from ssh
    %----------------------------------------------
    if type_detection==2 || type_detection==3

        psi = squeeze(sshh);

        % test if the grid is regular or not
        if min(xx(1,:)==xx(end,:)) && min(yy(:,1)==yy(:,end))
            CS2 = contourc(xx(1,:),yy(:,1),psi,H/2);
        else
            figure('visible','off')
            CS2 = contour(xx,yy,psi,H/2);
>>>>>>> ameda_v2
        end
    end

    % concantene all streamlines
    %----------------------------------------------
    CS = [CS1,CS2];

    % rearrange all the contours in C to the structure array 'isolines'
    %-----------------------------------------------------------
    [isolines,~] = scan_lines(CS);

    % scan streamlines and validate 'centers0' which are alone as potential centers
    %----------------------------------------------

    % Initialization
    j = 1; % first coordinates of the contour scan

    while j <= length(isolines)

        xdata = isolines(j).x; % vortex x's
        ydata = isolines(j).y; % vortex y's

        % conditions to have determine a ptential center
        % 1) a closed contour
        % 2) detected the right eddy center inside the polygon

        if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
                inpolygon(xy_ci,xy_cj,xdata,ydata) && length(xdata)>=n_min

            % searchs the coordinates of the centers in the contour 
            IN = inpolygon(xy_ctsi,xy_ctsj,xdata,ydata);
            [p,~] = find(IN==1); % index of the coordinates in the contour
            nc = length(p); % number of center in the contour

            % only one center include in the streamline
            if nc==1

                type1(ii) = centers_type(ii);
                x1(ii)    = xy_ci;
                y1(ii)    = xy_cj;
                i1(ii)    = C_I;
                j1(ii)    = C_J;

            % the contour contains more than 1 center
            elseif nc>1

                j = length(isolines); % stop the scan

            end
        end

        % increment the counter
        j = j + 1;
    end

end % centers ii loop

%----------------------------------------------
% export potential centers as 'centers'
%----------------------------------------------

% Initialisation
j = 1; % first coordinates of the contour scan

for k=1:length(type1)

    % remove centers with no close contour
    if ~isnan(type1(k))

        centers.type(j) = type1(k);
        centers.x(j)    = x1(k);
        centers.y(j)    = y1(k);
        centers.i(j)    = i1(k);
        centers.j(j)    = j1(k);
        j = j + 1;

    end
end
<<<<<<< HEAD
=======

disp(['    ',num2str(j-1),' potential centers found (',num2str(k-j+1),' max LNAM removed)'])
disp(' ')

>>>>>>> ameda_v2
