function [centers0,centers] = mod_eddy_centers(source,stp,fields)
%[centers0,centers] = mod_eddy_centers(source,stp,fields)
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
%  - bx: used to define smaller areas around each center and find which
%       ones are included (and alone) in a closed streamline. Must take into
%       account the resolution factor!
%  - DH: delta isocontour of ssh
%
%  For a description of the input parameters see mod_eddy_param.m.

%
%  There are two step to select potential center:
%
%  First the max(|LNAM(LOW<0)>K|) are computed.
%
%  Then the potential centers saved are the max LNAM surrounded by at least
%  two closed contour of ssh (or psi by using compute_psi)
%  and have a minimal and a maximal size defined by nRmin and nR_lim.
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
%   Ver. 3.1 2014 LMD from Nencioli et al. routines
%-------------------------
%
%=========================

plo=0; % debug mode

%---------------------------------------------
% read fields and initialisation
%---------------------------------------------

%----------------------------------------------
% load key_source and parameters (use mod_eddy_params.m first)
load('param_eddy_tracking')

%----------------------------------------
% load 2D velocity fields (m/s) for step stp
eval(['[x,y,mask,u,v,ssh] = load_fields_',source,'(stp,resol,deg);'])

% to calculate psi extrapole u and v to 0 in the land
u(isnan(u)) = 0;
v(isnan(v)) = 0;

%----------------------------------------
% initialise centers as structure
centers0 = struct('step',nan,'type',[],'x',[],'y',[],'i',[],'j',[]);
centers = centers0;

%---------------------------------------------
% Max LNAM 'centers0' for a given step k
%---------------------------------------------

disp([' Find potential centers step ',num2str(stp),' %-------------'])

centers0.step = stp;

%---------------------------------------------
% LNAM n OW criteria define contour including potential centers
OW = double(fields.LOW); % Okubo Weiss
LNAM = double(fields.LNAM); % LNAM
LOW = abs(LNAM);
LOW(OW>=0 | isnan(OW)) = 0;

% test if the grid is regular or not
if grid_reg
    CS = contourc(x(1,:),y(:,1),LOW,[K K]);
else
    HF = figure('visible','off');
    CS = contour(x,y,LOW,[K K]);
    close(HF)
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
        %Lm = LNAM.*maskin;
        Lm = LNAM;
        Lm(~in) = NaN;

        % L maximum value inside the contour(j)  
        if any(mask(:).*maskin(:)>0) && max(abs(Lm(:)))~=0

            LC = Lm(abs(Lm)==max(abs(Lm(:))));

            % save coordinates of the L maximum inside the contour(j)
            if mask(find(Lm==LC(1),1))==1

                centers0.type(k) = sign(LC(1));
                centers0.x(k)    = x(find(Lm==LC(1),1));
                centers0.y(k)    = y(find(Lm==LC(1),1));
                [centers0.j(k),centers0.i(k)] = find(Lm==LC(1),1);
                MaxL(k) = LC(1);
                

                % increment the counter
                k = k + 1; % next contour
            end
        end
    end

    % increment the counter
    j = j + n + 1; % series of coordinates of the next contour 
end

disp(['  -> ',num2str(k-1),' max LNAM found step ',num2str(stp)])

if k==1
    disp(['!!! WARNING !!! No LNAM extrema found - check the LNAM computation step ',num2str(stp)])
end

%----------------------------------------------
% Remove 'centers0' with not enough closed and too small streamlines around
% one or two max LNAM. Remove 
% Remove also 'centers0' below 'lat_min' of latitude (basically 5°)
%----------------------------------------------

disp(['  Remove max LNAM without 2 closed streamlines with proper size step ',num2str(stp)])

% all max LNAM coordinates for the given step stp
centers.step = stp;

type1 = nan(1,length(centers0.type));
x1    = type1;
y1    = type1;
j1    = type1;
i1    = type1;
second = type1; % second center to be tested

%----------------------------------------------
% compute each max LNAM in a smaller area
for ii=1:length(centers0.x)

    % fix the indice of the main max LNAM
    C_I = centers0.i(ii);
    C_J = centers0.j(ii);
    % max LNAM coordinates
    xy_ci = centers0.x(ii);
    xy_cj = centers0.y(ii);

    if ~grid_ll || (grid_ll && abs(xy_cj)>lat_min)
        
        % box size around the main max LNAM
        bx = bxi(C_J,C_I);
        % grid size
        Dx = abs(Dxi(C_J,C_I));
        % first deformation radius
        Rd = abs(Rdi(C_J,C_I));
        % coriolis parameter
        f = abs(f_i(C_J,C_I));

        % resize coordinate and velocity matrix 
        % (making sure not to go outside the domain)
        xx = x(max(C_J-bx,1):min(C_J+bx,size(x,1)), ...
            max(C_I-bx,1):min(C_I+bx,size(x,2)));
        yy = y(max(C_J-bx,1):min(C_J+bx,size(y,1)), ...
            max(C_I-bx,1):min(C_I+bx,size(y,2)));
        mk = mask(max(C_J-bx,1):min(C_J+bx,size(mask,1)), ...
            max(C_I-bx,1):min(C_I+bx,size(mask,2)));
        vv = v(max(C_J-bx,1):min(C_J+bx,size(v,1)), ...
            max(C_I-bx,1):min(C_I+bx,size(v,2)));
        uu = u(max(C_J-bx,1):min(C_J+bx,size(u,1)), ...
            max(C_I-bx,1):min(C_I+bx,size(u,2)));

        if type_detection>=2
            sshh = ssh(max(C_J-bx,1):min(C_J+bx,size(ssh,1)), ...
                max(C_I-bx,1):min(C_I+bx,size(ssh,2)));
        end

        % indice of the max LNAM in the small domain
        [cj,ci] = find(yy==xy_cj & xx==xy_ci);

        % indices of all max LNAM in the smaller area
        % (contains at least c_j and c_i)
        cts_j = zeros(length(centers0.y),1);
        cts_i = cts_j;
        for k=1:length(centers0.y)
            try
                [cts_j(k),cts_i(k)] = find(yy==centers0.y(k) & xx==centers0.x(k));
            catch
            end
        end

        % can be empty
        zero_centers = find(cts_j==0 | cts_i==0);
        cts_i(zero_centers) = [];
        cts_j(zero_centers) = [];

        % max LNAM positions in the smaller area
        if ~isempty(cts_i)
            xy_ctsi = zeros(length(cts_i),1);
            xy_ctsj = xy_ctsi;
            for k=1:length(cts_i)
                xy_ctsi(k) = xx(cts_j(k),cts_i(k));
                xy_ctsj(k) = yy(cts_j(k),cts_i(k));
            end
        else
            disp(['!!! ERROR !!! No LNAM maximum computable in the smaller area ',...
                'centers0 ',num2str(ii),' step ',num2str(stp)]) 
            xy_ctsi = [];
            xy_ctsj = [];
        end

        % initialize streamlines to be scanned
        CS1 = [];
        CS2 = [];

        %----------------------------------------------
        % compute the psi field from velocity fields
        if type_detection==1 || type_detection==3

            psi1 = double(compute_psi(xx,yy,mk,uu*f/g*1e3,vv*f/g*1e3,ci,cj,grid_ll));

            % determine contour to be scan
            H = floor(nanmin(psi1(:))):DH:ceil(nanmax(psi1(:)));    
            % impose a limit
            if length(H)>nH_lim
                H = nH_lim;
            end

            % test if the grid is regular or not
            if grid_reg
                CS1 = contourc(xx(1,:),yy(:,1),psi1,H);
            else
                HF = figure('visible','off');
                CS1 = contour(xx,yy,psi1,H);
                close(HF)
            end
        end

        %----------------------------------------------
        % compute the psi field from ssh
        if type_detection==2 || type_detection==3

            psi2 = double(squeeze(sshh));

            % determine contour to be scan
            Hs = floor(nanmin(psi2(:))):DH:ceil(nanmax(psi2(:)));    
            % impose a limit to the number of streamlines scanned
            if length(Hs)>nH_lim
                Hs = nH_lim;
            end

            % test if the grid is regular or not
            if grid_reg
                CS2 = contourc(xx(1,:),yy(:,1),psi2,Hs);
            else
                HF = figure('visible','off');
                CS2 = contour(xx,yy,psi2,Hs);
                close(HF)
            end
        end

        %----------------------------------------------
        % concantene all streamlines
        CS = [CS1,CS2];

        %-----------------------------------------------------------
        % rearrange all the contours in C to the structure array 'isolines'
        [isolines,~] = scan_lines(CS);

        %----------------------------------------------
        % scan streamlines and validate 'centers0' as potential centers

        % Initialization
        j = 1; % first coordinates of the contour scan
        radius = []; % successif radius of streamlines
        
        while j <= length(isolines)

            xdata = isolines(j).x; % vortex x's
            ydata = isolines(j).y; % vortex y's

            % conditions to have determine a ptential center
            % 1) two closed contour
            % 2) respect the maximum and minimal size
            % 3) keep the higher max LNAM center in case of 2 centers

            if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
                    inpolygon(xy_ci,xy_cj,xdata,ydata) && length(xdata)>=n_min

                % search the coordinates of the centers in the contour 
                IN = inpolygon(xy_ctsi,xy_ctsj,xdata,ydata);
                [p,~] = find(IN==1); % index of the coordinates in the contour
                nc = length(p); % number of center in the contour

                % only one center include in the streamline
                if nc==1

                    % increase number of streamlines with the radius value
                    R = mean_radius([xdata;ydata],grid_ll);
                    radius = [radius R(1)];

                    % record only center with 2 streamlines and proper size
                    if length(radius)>=2 && radius(end)>=nRmin*Dx && radius(end)<=nR_lim*Rd

                        disp(['   Validate max LNAM ',num2str(ii),...
                            ' with 2 streamlines at step ',num2str(stp)])

                        type1(ii) = centers0.type(ii);
                        x1(ii)    = xy_ci;
                        y1(ii)    = xy_cj;
                        i1(ii)    = C_I;
                        j1(ii)    = C_J;
 
                        second(ii) = 0; % no second center

                        j = length(isolines); % stop the scan

                    end

                % the contour contains 2 centers
                elseif nc==2 && ii>1

                    % increase number of streamlines with the radius value
                    R = mean_radius([xdata;ydata],grid_ll);
                    radius = [radius R(1)];

                    % second max LNAM centers indice
                    q = find(xy_ctsi(p)~=xy_ci | xy_ctsj(p)~=xy_cj);

                    % record only center with 2 streamlines and proper size
                    if length(radius)>=2 && abs(xy_ctsj(p(q)))>lat_min &&...
                            radius(end)>=nRmin*Dx && radius(end)<=nR_lim*Rd
                        
                        % index the second max LNAM if already scanned
                        jj = find(centers0.x==xy_ctsi(p(q)) &...
                            centers0.y==xy_ctsj(p(q)));
                        
                        % not validated
                        if isnan(x1(jj))
                                
                            % record second center indice
                            second(ii) = jj;
                            
                        % or validated as single eddy
                        elseif second(jj)==0

                            disp(['   Second max LNAM ',num2str(ii),...
                                ' in the same streamline ',num2str(jj),...
                                ' already validated at step ',num2str(stp)])

                        % or already recorded as double eddy
                        elseif second(jj)==ii

                            display(['   2 max LNAM in the same streamline: ',...
                                num2str(ii),' and ',num2str(jj),' step ',num2str(stp)])

                            % search for the higher LNAM
                            if MaxL(jj) > MaxL(ii)

                                display(['    -> validate center ',num2str(jj),...
                                    ' step ',num2str(stp)])
                                second(ii) = NaN;

                            else

                                display(['    -> validate center ',num2str(ii),...
                                    ' step ',num2str(stp)])
                                second(ii) = jj; % second center indice
                                second(jj) = NaN;

                            end

                        % debug warning
                        elseif second(jj)~=ii

                            disp(['   !!! ERROR !!! Second max LNAM ',num2str(jj),...
                                ' include center ',num2str(second(jj)),...
                                ' instead of ',num2str(ii),' at step ',num2str(stp)])

                        end

                        j = length(isolines); % stop the scan
                        
                    end
                    
                % the contour contains more than 2 centers
                elseif nc>2
                    j = length(isolines); % stop the scan
                end

            end

            % increment the isolines counter
            j = j + 1;
        end

        % record max LNAM or max LNAM above latitude 'lat_min' as eddy center
        if ~isnan(second(ii))
            
            type1(ii) = centers0.type(ii);
            x1(ii)    = centers0.x(ii);
            y1(ii)    = centers0.y(ii);
            i1(ii)    = centers0.i(ii);
            j1(ii)    = centers0.j(ii);
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the potential centers detected in the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plo
    figure
    if type_detection==1 || type_detection==3
        contour(xx,yy,psi1,H,'k')
    end
    hold on
    if type_detection==2 || type_detection==3
        contour(xx,yy,psi2,Hs,'r')
    end
    quiver(xx,yy,uu,vv,'k')

    plot(xy_ctsi,xy_ctsj,'r*')
    plot(xy_ci,xy_cj,'k*')

    if ~isnan(x1(ii))
        plot(x1(ii),y1(ii),'og');
    end

    hold off
    title(['Center ',num2str(ii)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

end % centers ii loop

%----------------------------------------------
% resolve double eddies conflicts and export potential centers as 'centers'
%----------------------------------------------

if ~isempty(second==0)

    IND = find(second==0);

    % double eddies with second center in the streamline
    IND1 = second(second>0);

    if ~isempty(IND1)
    
        % double eddies with a weaker center in the streamline
        IND2 = IND1(isnan(second(IND1)));
        for k=1:length(IND2)
            disp (['   Remove max LNAM ',num2str(IND2(k)),...
                ' at step ',num2str(stp)])
            IND = [IND find(second==IND2(k))];
        end

    end

else

    disp (['!!! WARNING or ERROR !!! No potential centers found - ',...
        'check the streamlines scanning process at step ',num2str(stp)])
    stop

end

IND = sort(IND);
centers.type = type1(IND);
centers.x    = x1(IND);
centers.y    = y1(IND);
centers.i    = i1(IND);
centers.j    = j1(IND);

disp([' Potential eddy centers found step ',num2str(stp)])
disp(['  -> ',num2str(length(IND)),' potential centers found'])
disp(['    (',num2str(length(type1)-length(IND)),' max LNAM removed)'])

disp(' ')
