function [centers2,shapes1,shapes2,profil2,warn_shapes,warn_shapes2] = ...
    mod_eddy_shapes(source,stp,fields,centers)
%[centers2,shapes1,shapes2,profil2,warn_shapes,warn_shapes2] =
%       mod_eddy_shapes(source,stp,fields,centers)
%
%  Computes the shapes (if any) of eddies identified with their potential
%  centers by mod_eddy_centers.m and stored as {centers(t)} in
%  [path_out,'eddy_centers',runname]
%
% - 'source' is the type of netcdf data (AVISO, NEMO, ROMS,...) that
%   determine which load_field_'source'is  used in the routine
% - 'fields' is the step 'stp' of the detection_fields computed with
%   mod_eddy_fields.m
%
%  For a description of the input parameters see mod_eddy_param.m

%
%  Then [path_out,'eddy_centers',runname] is updated with the structure
%  {centers2(t)} including coordinates of validated centers as eddies.
%  Fields of the main centers(indice 1) and double centers(indice 2),
%  if any, are:
%  - x1,y1(n): coordinates of the main centers
%  - x2,y2(n): coordinates of the double center if any (NaN otherwise)
%  - ind2(n): indice in centers2 of the second center in interaction
%
%  Also features of validated eddy are saved/updated in
%  [path_out,'eddy_shapes',runname] as structure arrays:
%  - {shapes1(t)} for the main center and
%  - {shapes2(t)} for double eddies (if any).
%
%  Features inside these structure include: 
%  - xy(n): xy is a [2xm double] cell array, containing
%    the coordinates of the m vertices that define the
%    boundaries of the n-th eddy of the t-th step
%  - velmax(n): velmax is the maximum mean velocity along the contour xy
%  - deta(n): deta is the difference between the ssh on the
%    limiting contour and the ssh in the extremum in the middle (max or min)
%  - taumin(n): tau is the minimal eddy turnover time inside
%    the eddy bourndaries
%  - nrho(n): part of the local curvature along the coutour n
%  - rmax(n): rmax is the averaged radius from the eddy centre to
%    the boundary vertices (see 'mean_radius' for more details)
%  - aire(n): aire is the area delimited by the eddy boundaries
%    (see 'mean_radius.m' for more details)
%  - xbary,ybary: barycenter of the ellipse
%  - ellip,theta: the ellipse features fitted on the velmax contour
%  - xy_end, vel_end, deta_end, nrho_end, r_end, aire_end: features for the
%    last contour with one center (shapes1 only)
%
%  Also:
%
%  if extended_diags =1, you save in {shapes1 or 2}:
%  - ke,vort,VortM,OW,LNAM: fields inside the eddy
%
%  if streamlines =1, you save all streamlines scanned during the n-th
%  shapes computation at t=daystreamfunction in {profil2(t)}, with:
%  - step,nc{n}: time step and number of centers at t
%  - eta{n},rmoy{n},vel{n},tau{n},nrhoi{n}: profil scanned
%  - myfit: fitting result with: alpha(n), rsquare(n), rmse(n) which are
%    fitting result from compute_best_fit which are also saved in shapes1
%
%  Another file which contains information on the process to  compute eddy 
%  shapes is saved as [path_out,'warnings_shapes',runname] in
%  {warn_shapes2(t)}:
%  - no_curve(n): 1 if no closed contour of PSI was found
%    around the eddy center;
%  - f(n): Coriolis parameter value at the center
%  - Rd(n): first baroclinic deformation radius of Rossby at the center
%  - gama(n): resolution factor at the center (Rd/DX/deg*res)
%  - bx(n): area where the shape is computed (bx*fac)
%  - calcul_curve(n): 1 if the closed contour as been obtain
%    through the computation of psi field thanks to 'compute_psi'.
%    When both (ssh and psi) are computed, psi is choosen in the
%    case ssh do not allowed to close contour around the center
%  - large_curve1(n): 1 if the velocity of contours outside the velmax
%    contour stay higher than 97% (see k_vel_decay in mod_eddy_param)
%    of the velmax. Then contour is defined as a 'gyre'
%  - large_curve2(n): 1 if the velocity of contours outside the double
%    contour (shapes2) contour stay higher than 97% (see k_vel_decay in
%    mod_eddy_param) of the velmax. Then contour is defined as a 'gyre';
%  - too_weak2(n): 1 if double eddy shapes is removed from the detection
%    due to the presence of 2 double contour or one of the eddies is too
%    small (no velmax recorded).
%
%  (t is the time step index; n is the indice of eddy detected at t)
%
%
%  The function returns also a log file
%  [path_out,'log_eddy_shapes',runname,'.txt'] which contains additional
%  information on the process to compute eddy shapes.
%  (NOTE: if the file already exist, the new log file will be append to it)
%  
%  'eddy_dim' is used to compute eddy shapes. Check the documentation in 
%  eddy_dim.m for further details.
%
%-------------------------
%   Ver. 3.2 Apr.2015 Briac Le Vu
%   Ver. 3.1 2014 LMD from Nencioli et al. routines
%-------------------------
%
%=========================

%---------------------------------------------
% read fields and initialisation
%---------------------------------------------

disp(['Start step ',num2str(stp),' %-------------'])
disp(' ')

%----------------------------------------------
% load key_source and parameters
load('param_eddy_tracking')

%----------------------------------------
% load 2D velocity fields (m/s) for step stp
eval(['[x,y,mask,uu,vv,sshh] = load_fields_',source,'(stp,resol,deg);'])

%----------------------------------------
% initialise centers as structure
centers2 = struct('step',[],'type',[],'x1',[],'y1',[],'x2',[],'y2',[],'dc',[],'ind2',[]);

shapes1 = struct('step',[],'xy',[],'velmax',[],'deta',[],'taumin',[],'nrho',[],...
                'rmax',[],'aire',[],'xbary',[],'ybary',[],'ellip',[],'theta',[],...
                'xy_end',[],'vel_end',[],'deta_end',[],'nrho_end',[],'r_end',[],'aire_end',[]);
shapes2 = struct('step',[],'xy',[],'velmax',[],'deta',[],'taumin',[],'nrho',[],...
                'rmax',[],'aire',[],'xbary',[],'ybary',[],'ellip',[],'theta',[]);

if streamlines
    profil2 = struct('step',[],'nc',[],'eta',[],'rmoy',[],'vel',[],'tau',[],'nrhoi',[],'myfit',[]);

    struct1 = struct('alpha',[],'rsquare',[],'rmse',[]);
    names1 = [fieldnames(shapes1); fieldnames(struct1)];
    shapes1 = cell2struct([struct2cell(shapes1); struct2cell(struct1)], names1, 1);
end

if extended_diags==1
    struct1 = struct('ke',[],'vort',[],'vortM',[],'OW',[],'LNAM',[]);
    names1 = [fieldnames(shapes1); fieldnames(struct1)];
    shapes1 = cell2struct([struct2cell(shapes1); struct2cell(struct1)], names1, 1);
end

% shapes struct for the possible second shape with 2 centers
warn_shapes = struct('no_curve',[],'f',[],'Rd',[],'gama',[],'bx',[],'calcul_curve',[],...
                'large_curve1',[],'large_curve2',[],'too_weak2',[]);
warn_shapes2 = warn_shapes;

%----------------------------------------------
% Compute eddy shape

if streamlines
    profil2.step = centers.step;
else
    profil2 =   [];
end

centers2.step = centers.step;
centers2.ind2 = nan(1,length(centers.type));
shapes1.step = centers.step;
shapes2.step = centers.step;
shapes1.xy = {};
shapes1.xy_end = {};
shapes1.velmax = nan(1,length(centers.type));
shapes2.xy = {};
shapes2.velmax = nan(1,length(centers.type));
warn_shapes.no_curve = ones(1,length(centers.type));

% structures names 
%----------------------------------------------
shapes_name = fieldnames(shapes1);
n1 = find(strcmp(shapes_name,'xy_end'));

%----------------------------------------------
% loop through all centers detected for a given step
for ii=1:length(centers.type)

    disp([' === Center ',num2str(ii),' ==='])

    % center indice
    c_j = centers.j(ii);
    c_i = centers.i(ii);

    % initialization
    bound = 1; % flag that indicates that permit to increase the small area
    fac = 0; % increase factor for the area where PSI is computed
    tmp_large = [1 NaN]; % flag on: No maximum found for single and double eddy
    tmp_CD = []; % centers coordinates in case of double eddy 
    tmp_xy = cell(1,3); % contour for single (max and final) and double eddy
    tmp_allines = []; % streamlines features to be tested every eddy_dim computation
    tmp_rmax = nan(1,3); % radius to be tested every eddy_dim computation
    tmp_velmax = zeros(1,3); % velocity to be tested every eddy_dim computation
    tmp_tau = nan(1,2); % turnover time
    tmp_deta = nan(1,3); % delta ssh 
    tmp_nrho = nan(1,3); % local curvature along the single and double eddy 
    
    while bound
        % factor which increases the area where psi is computed
        % if eddy dimensions are too big (bound=1)
        fac = fac + 1; % this determine larger area

        %----------------------------------------------
        % xy is the computed shape;
        % the others are flags output by eddy_dim, and saved in 'warnings_shapes';
        [CD,xy,allines,rmax,velmax,tau,deta,nrho,large,warn,calcul] =...
            eddy_dim(uu,vv,sshh,mask,x,y,centers,ii,f_i(c_j,c_i),...
            Rdi(c_j,c_i),fac*bxi(c_j,c_i));

        %----------------------------------------------
        % flags exploitation
        if warn
            disp('    -> No significant streamlines closed around the center')
            bound = 0;
        else
            % temporary save eddy_dim(1) until a true eddy is found
            if velmax(1) > tmp_velmax(1)*(1 + epsil) || large(1) < tmp_large(1)
                if tmp_large(1) == 1 || velmax(1)
                    tmp_large(1)  = large(1);
                    tmp_xy(1)     = xy(1);
                    tmp_rmax(1)   = rmax(1);
                    tmp_velmax(1) = velmax(1);
                    tmp_tau(1)    = tau(1);
                    tmp_deta(1)   = deta(1);
                    tmp_nrho(1)   = nrho(1);
                    tmp_allines   = allines;
                end
            end
            % if last contour velocity is still increasing then temporary save eddy_dim(3)
            if isnan(tmp_deta(3)) || abs(deta(3)) > abs(tmp_deta(3))*(1 + epsil) && rmax(3) >= tmp_rmax(1)
                bound = 1;
                % temporary save flags
                tmp_xy(3)     = xy(3);
                tmp_rmax(3)   = rmax(3);
                tmp_velmax(3) = velmax(3);
                tmp_nrho(3)   = nrho(3);
                tmp_deta(3)   = deta(3);
            else
                % stop if no significant increasing of last contour size
                bound = 0;
            end
            % if double eddy not found yet
            if isnan(large(2)) && fac==1
                % test double eddy in larger area
                bound = 1;
            % if double contour is still increasing then temporary save eddy_dim(2)
            elseif velmax(2) > tmp_velmax(2)*(1 + epsil)
                bound = 1;
                tmp_large(2)  = large(2);
                tmp_CD        = CD;
                tmp_xy(2)     = xy(2);
                tmp_rmax(2)   = rmax(2);
                tmp_velmax(2) = velmax(2);
                tmp_deta(2)   = deta(2);
                tmp_nrho(2)   = nrho(2);
                tmp_tau(2)    = tau(2);
                tmp_allines   = allines;
            elseif ~isnan(large(2))
               % stop after 2 scans if no significant increase of the velocity
               bound = 0;
            end
        end
        % if no closed curve more intense in the larger area then
        % final eddy shape is the one computed in the smaller area
        if bound
            disp(['    Big eddy: going to fac = ',num2str(fac+1)])
        elseif fac > 1
            disp(['    No bigger eddy nor interaction stop dezoom at fac ',num2str(fac)])
            % take the last saved eddy_shape
            large   = tmp_large;
            CD      = tmp_CD;
            xy      = tmp_xy;
            allines = tmp_allines;
            rmax    = tmp_rmax;
            velmax  = tmp_velmax;
            tau     = tmp_tau;
            deta    = tmp_deta;
            nrho    = tmp_nrho;
        end
    end % end while loop

    %----------------------------------------------
    % write out which kind of eddy found
    if isempty(xy{1})
        disp('    -> No Eddy')
        xy = cell(1,3);
    elseif rmax(1) < nRmin*Dxi(c_j,c_i) || rmax(3) < nRmin*Dxi(c_j,c_i)
        disp('    -> Too small Eddy')
        xy = cell(1,3);
    elseif large(2) == 0
        disp('    -> Eddy with 2 centers')
    elseif large(1) == 0
        disp('    -> Eddy with 1 center')
    elseif large(2) == 1
        disp('    -> Largest with 2 centers')
    elseif large(1) == 1
        disp('    -> Largest with 1 center')
    end

    %----------------------------------------------------------
    % save eddy_dim results in a struct array shapes end the single eddies
    if ~isempty(xy{3})
        centers2.type(ii) = centers.type(ii);
        centers2.x1(ii)   = centers.x(ii);
        centers2.y1(ii)   = centers.y(ii);
        if rmax(1)>rmax(3)
            disp(['!!! ERROR !!! Rmax > Rend (+',...
                num2str(round(rmax(1)/rmax(3)*100-100)),...
                '%) center ',num2str(ii),' step ',num2str(stp)])
            xy(3) = xy(1);
            rmax(3) = rmax(1);
            velmax(3) = velmax(1);
            deta(3) = deta(1);
            nrho(3) = nrho(1);
        end
        if abs(deta(1))>abs(deta(3))
            disp(['!!! ERROR !!! |Deta_Rmax| > |Deta_Rend| center ',...
                num2str(ii),' step ',num2str(stp)])
            xy(3) = xy(1);
            rmax(3) = rmax(1);
            velmax(3) = velmax(1);
            deta(3) = deta(1);
            nrho(3) = nrho(1);
        end 
        shapes1.xy_end(ii)   = xy(3);
        shapes1.vel_end(ii)  = velmax(3);
        shapes1.deta_end(ii) = deta(3);
        shapes1.nrho_end(ii) = nrho(3);
        shapes1.r_end(ii)    = rmax(3);
        shapes1.aire_end(ii) = pi*rmax(3)*rmax(3);
    else
        disp(['!!! WARNING !!! No contour recorded center ',...
            num2str(ii),' step ',num2str(stp)])
        names = fieldnames(centers2);
        for n=2:length(names)
            centers2.(names{n})(ii) = NaN;
        end
        shapes1.xy{ii} = NaN;
        names = fieldnames(shapes1);
        for n=3:n1-1
            shapes1.(names{n})(ii) = NaN;
        end
        shapes1.xy_end{ii} = NaN;
        for n=n1+1:length(names)
            shapes1.(names{n})(ii) = NaN;
        end
        if streamlines
%            if find(daystreamfunction==stp)
                names = fieldnames(profil2);
                for n=2:length(names)-1
                    profil2.(names{n}){ii} = NaN;
                end
                profil2.myfit(ii).curve = NaN;
                profil2.myfit(ii).err   = NaN;
%            end
        end
    end
    
    %----------------------------------------------------------
    % save eddy_dim results in a struct array shapes1 the single eddies
    if ~isempty(xy{1})
        shapes1.xy(ii)     = xy(1);
        shapes1.velmax(ii) = velmax(1);
        shapes1.taumin(ii) = tau(1);
        shapes1.deta(ii)   = deta(1);
        shapes1.nrho(ii)   = nrho(1);
        shapes1.rmax(ii)   = rmax(1);
        shapes1.aire(ii)   = pi*rmax(1)*rmax(1);
        [xbary,ybary,~,a,b,theta,~] = compute_ellip(xy{1},grid_ll);
        shapes1.xbary(ii) = xbary;
        shapes1.ybary(ii) = ybary;
        if a>=b && a~=0
            shapes1.ellip(ii) = 1-b/a;
        elseif a<b && b~=0
            shapes1.ellip(ii) = 1-a/b;
            if theta > pi/2
                theta = theta-pi/2;
            else
                theta = theta+pi/2;
            end
        else
            shapes1.ellip(ii) = NaN;
            shapes1.theta(ii) = NaN;
        end
        shapes1.theta(ii) = theta;
        if extended_diags==1
            in_eddy = inpolygon(x,y,xy{1}(1,:),xy{1}(2,:));
            ke   = fields.ke(in_eddy~=0);
            vort = fields.vort(in_eddy~=0);
            OW   = fields.OW(in_eddy~=0);
            LNAM = fields.LNAM(in_eddy~=0);
            shapes1.ke(ii)   = nansum(ke(:));
            shapes1.vort(ii) = nanmean(vort(:));
            if nanmean(vort(:)) > 0
                shapes1.vortM(ii) = nanmax(vort);
            elseif nanmean(vort(:)) < 0
                shapes1.vortM(ii) = nanmin(vort);
            else
                display('!!! WARNING !!! vort==0 - check vorticity computation')
            end
            shapes1.OW(ii)   = nanmean(OW(:));
            shapes1.LNAM(ii) = nanmean(LNAM(:));
        end
        
        %----------------------------------------------------------
        % save the streamlines features at i=daystreamfunction
        if streamlines
            if find(daystreamfunction==stp)
                names = fieldnames(profil2);
                for n=2:length(names)-1
                    profil2.(names{n}){ii} = allines(:,n-1)';
                end
                % compute and save curve fitting stats at step i for the eddy ii
                [curve,err] = compute_best_fit(allines,rmax(1),velmax(1));
                if ~isempty(curve)
                    profil2.myfit(ii).curve = curve;
                    profil2.myfit(ii).err   = err;
                    shapes1.alpha(ii)   = curve.a;
                    shapes1.rsquare(ii) = err.rsquare;
                    shapes1.rmse(ii)    = err.rmse;
                else
                    profil2.myfit(ii).curve = NaN;
                    profil2.myfit(ii).err   = NaN;
                    shapes1.alpha(ii)   = NaN;
                    shapes1.rsquare(ii) = NaN;
                    shapes1.rmse(ii)    = NaN;
                end
            else
                names = fieldnames(profil2);
                for n=2:length(names)-1
                    profil2.(names{n}){ii} = NaN;
                end
                profil2.myfit(ii).curve = NaN;
                profil2.myfit(ii).err   = NaN;
                shapes1.alpha(ii)   = NaN;
                shapes1.rsquare(ii) = NaN;
                shapes1.rmse(ii)    = NaN;
            end
        end
    else
        disp(['No contour with maximum velocity center ',...
            num2str(ii),' step ',num2str(stp)])
        shapes1.xy{ii} = NaN;
        names = fieldnames(shapes1);
        for n=3:n1-1
            shapes1.(names{n})(ii) = NaN;
        end
        for n=n1+6:length(names)
            shapes1.(names{n})(ii) = NaN;
        end

    end
    
    %----------------------------------------------------------
    % save eddy_dim results in a struct array shapes2 double eddies
    if ~isempty(xy{2})
        if CD(1,1)==centers.x(ii) && CD(2,1)==centers.y(ii)
            centers2.x2(ii) = CD(1,2);
            centers2.y2(ii) = CD(2,2);
        elseif CD(1,2)==centers.x(ii) && CD(2,2)==centers.y(ii)
            centers2.x2(ii) = CD(1,1);
            centers2.y2(ii) = CD(2,1);
        end
        if grid_ll
            centers2.dc(ii) = sw_dist2(CD(2,:),CD(1,:),'km');
        else
            centers2.dc(ii) = sqrt(diff(CD(1,:)).^2 + diff(CD(2,:)).^2); % km
        end
        shapes2.xy(ii)     = xy(2);
        shapes2.velmax(ii) = velmax(2);
        shapes2.taumin(ii) = tau(2);
        shapes2.deta(ii)   = deta(2);
        shapes2.nrho(ii)   = nrho(2);
        shapes2.rmax(ii)   = rmax(2);
        shapes2.aire(ii)   = pi*rmax(2)*rmax(2);
        [xbary,ybary,~,a,b,theta,~] = compute_ellip(xy{2},grid_ll);
        shapes2.xbary(ii) = xbary;
        shapes2.ybary(ii) = ybary;
        if a>=b && a~=0
            shapes2.ellip(ii) = 1-b/a;
        elseif a<b && b~=0
            shapes2.ellip(ii) = 1-a/b;
            if theta > pi/2
                theta = theta-pi/2;
            else
                theta = theta+pi/2;
            end
        else
            shapes2.ellip(ii) = NaN;
        end
        shapes2.theta(ii) = theta;
    else
        names = fieldnames(centers2);
        for n=5:length(names)
            centers2.(names{n})(ii) = NaN;
        end
        shapes2.xy{ii} = NaN;
        names = fieldnames(shapes2);
        for n=3:length(names)
            shapes2.(names{n})(ii) = NaN;
        end
    end

    %----------------------------------------------------------
    % warnings from shape computation
    warn_shapes.no_curve(ii)     = warn;
    warn_shapes.f(ii)            = f_i(c_j,c_i);
    warn_shapes.Rd(ii)           = Rdi(c_j,c_i);
    warn_shapes.gama(ii)         = gamai(c_j,c_i);
    warn_shapes.bx(ii)           = bxi(c_j,c_i)*fac;
    warn_shapes.calcul_curve(ii) = calcul;
    warn_shapes.large_curve1(ii) = large(1);
    warn_shapes.large_curve2(ii) = large(2);
    warn_shapes.too_weak2(ii)    = 0;

end % ii=last eddy

if isempty(centers2.type)
    disp(['!!! WARNING or ERROR !!! No centers validated step ',...
        num2str(stp),' - check the contour selection'])
    %----------------------------------------------------------
    % Initialize warn_shape for eddies with shapes1
    warn_shapes.no_curve = [];
    warn_shapes2 = warn_shapes;
else
    %----------------------------------------------------------
    % Initialize warn_shape for eddies with shapes1
    warn_shapes2 = warn_shapes;

    %----------------------------------------------------------
    % remove the shapes and centers where no closed contour was found
    replace = find(isnan(shapes1.vel_end));
    names = fieldnames(centers2);
    for n=2:length(names)
        centers2.(names{n})(replace)=[];
    end
    names = fieldnames(shapes1);
    for n=2:length(names)
        shapes1.(names{n})(replace)=[];
    end
    names = fieldnames(shapes2);
    for n=2:length(names)
        shapes2.(names{n})(replace)=[];
    end
    names = fieldnames(warn_shapes2);
    for n=1:length(names)
        warn_shapes2.(names{n})(replace)=[];
    end
    if streamlines
%        if find(daystreamfunction==stp)
            names = fieldnames(profil2);
            for n=2:length(names)
                profil2.(names{n})(replace)=[];
            end
%        end
    end

    %----------------------------------------------------------
    % Remove shapes2 too weak:
    % if 2 shapes1 exist
    %   if 2 shapes2 exist remove the weakest
    % else
    %   replace shapes1 by shapes2
    % end
    %----------------------------------------------------------

    for ii=1:length(shapes2.velmax)

        %----------------------------------------------------------
        % test shapes2(ii) exists
        if ~isnan(shapes2.velmax(ii))
            
            mv=0;
            
            % find indice of the other center
            ind = find(centers2.x1==centers2.x2(ii) &...
                   centers2.y1==centers2.y2(ii));
               
            %----------------------------------------------------------
            % check that double centers are different type (debug)
            if centers2.type(ii)~=centers2.type(ind)
                
                % print out for debugging
                disp(['!!! ERROR !!! Double eddy ',...
                    num2str(ii),' step ',num2str(stp),...
                    ' mistaken around 2 different type of eddy'])
                mv=1;
                
            %----------------------------------------------------------
            % shapes1(ind) exists
            elseif ~isempty(ind)

                %----------------------------------------------------------
                % shapes1(ind) exists and its speed or shapes1(ii) speed
                % is lower than shapes2(ii) speed
                if shapes1.velmax(ind) < shapes2.velmax(ii) ||...
                    shapes1.velmax(ii) < shapes2.velmax(ii)
                
                    disp(['   Double eddy ',num2str(ii),' interact with eddy ',...
                        num2str(ind),' step ',num2str(stp)])

                    centers2.ind2(ii) = ind; % index the second center
                    centers2.ind2(ind) = ii; % index the second center

                    % test shapes2(ind) exists
                    if ~isnan(shapes2.velmax(ind))

                        disp(['   -> Double eddy ',num2str(ii),' with eddy ',num2str(ind),...
                            ' step ',num2str(stp),' is redundant'])

                        % test shapes2 both calculated or both not calculated
                        if warn_shapes2.calcul_curve(ii) == warn_shapes2.calcul_curve(ind)

                            %----------------------------------------------------------
                            % remove shapes2 if weakest
                            if shapes2.velmax(ii) < shapes2.velmax(ind)
                                mv = 1;
                            elseif shapes2.velmax(ii) == shapes2.velmax(ind)
                                if abs(shapes2.deta(ii)) <= abs(shapes2.deta(ind))
                                    mv = 1;
                                end
                            end
                            if mv
                                disp(['    -> Double eddy ',num2str(ii),...
                                    ' weaker than double eddy ',num2str(ind),...
                                    ' at step ',num2str(stp)])
                            end

                        else

                            %----------------------------------------------------------
                            % remove shapes2(ii) if calculated
                            if warn_shapes2.calcul_curve(ii)==1
                                disp(['    -> Calculated double eddy ',num2str(ii),...
                                    ' removed at step ',num2str(stp)])
                                mv= 1;
                            end

                        end

                    end

                %----------------------------------------------------------
                % shapes1(ind) exists but shapes2(ii) speed is too low
                else

                    % print out for debugging
                    disp(['   Double eddy ',num2str(ii),' with eddy ',num2str(ind),...
                        ' too weak in speed at step ',num2str(stp)])
                    mv=1;

                end

            %----------------------------------------------------------
            % if shapes1(ind) doesn't exist remove x2,y2(ii) then
            % copy shapes2 to shapes1_end and also shapes1_max if relevant
            else

                disp(['   Double eddy ',num2str(ii), ' step ',num2str(stp),...
                    ' contains a second core too small'...
                    ' and get a new last contour!!!'])
                mv= 1;

                % remove centers2 if no second eddy
                names = fieldnames(centers2);
                for n=5:length(names)
                    centers2.(names{n})(ii) = NaN;
                end
                
                % replace shapes1_end by mistaken shapes2
                names1 = fieldnames(shapes1);
                names2 = fieldnames(shapes2);
                for n=2:4
                    shapes1.(names1{n+n1-2})(ii) = shapes2.(names2{n})(ii);
                end
                for n=6:8
                    shapes1.(names1{n+n1-3})(ii) = shapes2.(names2{n})(ii);
                end
                
                % change a gyre to a true eddy of necessary
                if ( shapes2.velmax(ii) < k_vel_decay*shapes1.velmax(ii) ||...
                            warn_shapes2.large_curve2(ii) == 0 ) &&...
                            warn_shapes2.large_curve1(ii) == 1

                    % shapes2 met a speed decrease > 3%
                    disp(['   -> Double eddy ',num2str(ii), ' step ',num2str(stp),...
                        ' becomes true eddy!!!'])
                    warn_shapes2.large_curve1(ii) = 0;
                    
                end
                
                % replace also shapes1_max if shapes2 proper shapes1
                if shapes2.velmax(ii) > shapes1.velmax(ii) &&...
                        shapes2.rmax(ii) < nR_lim*Rdi(c_j,c_i) &&...
                        shapes2.nrho(ii) < nrho_lim &&...
                        warn_shapes2.large_curve2(ii) <= warn_shapes2.large_curve1(ii)

                    % shapes2 is a proper speed radius
                    disp(['   -> Double eddy ',num2str(ii), ' step ',num2str(stp),...
                        ' becomes also new speed radius!!!'])
                    for n=2:12
                        shapes1.(names1{n})(ii) = shapes2.(names2{n})(ii);
                    end

                end

            end

            %----------------------------------------------------------
            % at the end remove shapes2 too weak or mistaken as explain above
            if mv

                disp(['   Double eddy ',num2str(ii),' with eddy ',num2str(ind),...
                    ' step ',num2str(stp),' is removed!!!'])
                shapes2.xy{ii} = NaN;
                names = fieldnames(shapes2);
                for n=3:length(names)
                    shapes2.(names{n})(ii) = NaN;
                end
                warn_shapes2.too_weak2(ii) = 1;
                
            end
        end
    end
    
    disp(' ')
    disp(['-------------% End step ',num2str(stp)])
    disp(' ')

end

