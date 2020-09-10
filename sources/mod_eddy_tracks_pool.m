function mod_eddy_tracks_pool(name,update)
%mod_eddy_tracks(name {,update})
%
%  Computes eddy tracking from the eddies features, and saves them
%  as {eddy_tracks(n)} in [path_out,'eddy_tracks_',name].
%
% This routine will use:
%  - Dt (days) is the tolerance of time steps after eddy disapears
%  - update is a flag allowing to update an existing detection matrice:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the {shapes} time serie
%  - V_eddy is the center maximal speed
%  - D_stp, N_can: see mod_eddy_param for more details
%
%  For a description of the input parameters see param_eddy_tracking.m

%
%  Eddy tracks are connected by comparing the detected eddy fields at
%  sucessive time steps. An eddy in the time 't+dt' is assumed to be the
%  new position of an eddy of the same type detected at 't', if the two 
%  centers are found within an area of 'V_eddy*(1+dt)/2 + rmaxi + rmaxj' km
%  centered around the center position at time 't'. Where rmaxj is, for each
%  eddy, his mean radius averaged for the last 5 time steps tracked and
%  rmaxi the new eddy raidus to be tested.
%
%  If no eddies is connected after 't+Dt' the eddy is considered dissipated and
%  the track closed.
%
%  In case two or more eddies are found within the same area, the track is
%  connect to the centers by minimizing the cost funtion of an assignment matrice
%  considering all the new centers at t. The cost function is a NxM matrix:
%       C = sqrt ( d/D ² + dR/Rmax ² + dRo/Ro ² + dt/Dt/2 ²)
%   with N in row is the number of eddies in the past
%        M in colum is the number of new eddies
%
%  
%  Non filtrated tracked eddies are saved/updated as the structure array {traks(n)}
%  in [path_out,'eddy_tracks_',name'] with fields following {centers}
%  and {shapes}, where 'n' is the total number of track recorded and
%  'm' is the step index of a given track:
%  - step(m): steps when the eddy centers forming the track where
%                       detected;
%  - type(m): type of eddy;
%  - x1(m) and y1(m): x and y of the eddy centers forming the track;
%  - shapes1{m} : shapes of the eddies forming the track; (each 
%           cell is a 2xL array containing the x (first row) and 
%           y (second row) coordinates of the 'L' vertices defining
%           the eddy shape for the 'm'-th step.
%   ...
%{
second center in interaction location:
    x2
    y2
    dc is the distance between the two centers
    ind2 is the indice of the second center in shapes2
features of the Vmax profil:
    velmax1 is the maximum average speed along the contour
    tau1 is the minimal revolution period
    deta1 is the extremum delta ssh of the contour
    rmax1 is the radius of the circle equivalent to aire1
    aire1 is the area of the contour
features of the Vend profil as the Vmax profil definition
    velmax3
    tau3
    deta3
    rmax3
    aire3
if streamlines=1, features of the V-R profil:
    alpha
    rsquare
    rmse
if extended_diags=1, additional diags:
    xbary1 and ybary1 are the centroid center of the contour
    ellip1 is the 1-a/b ellipticity
    ke1
    vort1
    vortM1
    OW1
    LNAM1
features and additional diags of the interacting contour
    shapes2
    velmax2
    tau2
    deta2
    rmax2
    aire2
    xbary2
    ybary2
    ellip2
features of eddy in the grid
    Rd is the first deformation radius
    game is the number of grid point to define Rd
flag of the contour
    calcul is 1 if calculated from uv
    large1 is 1 if no maximum is found on V-R profil
    large2 is 1 if no maximum is found for interacting contour
    weak is 1 if a too weak interacting contour in front of velmax1 has been removed
indice of the interaction
    ind of x1,y1 in shapes1
    interaction indice of the interacting eddy in tracks
    interaction2 indice of extra interacting eddy in tracks
flag about the merging and splitting
    split is 1 if eddy come from a splitting
    merge is 1 if eddy end by a merging
    split2 same for extra eddy
    merge2 same for extra eddy
%}  
%
% Also:
%
%  warn_tracks contains flag during the process of filtration
%  which contains warnings relative to the presence of multiple eddies
%  candidats during the tracking within the searching area to determine the tracks
%  and the potential interaction.
%
%  The same information is saved in [path_out,'warnings_tracks_',name], 
%  an Rx3 array (R is the total number of warnings), where each column
%  represent in the order:
%  - step the warning occurred;
%  - n indice of eddy concerned;
%  - number of eddies detected within the area;
%
%-------------------------
%   Ver. 3.2 Feb.2016 Briac Le Vu
%   Ver. 3.1 2014 LMD
%-------------------------
%
%=========================

% load key_source and parameters (use mod_eddy_params.m first)
%----------------------------------------------

load('param_eddy_tracking')

% explicitly define constant for the parfor loop
% !! parfor is slower above 1 or 2 years or detection !!
extended_diags=extended_diags;grid_ll=grid_ll;
dps=dps;Dt=Dt;D_stp=D_stp;T=T;f=f;g=g;V_eddy=V_eddy;N_can=N_can;

% No update by default
if nargin==1
    update = 0;
end

% begin the log file
diary([path_out,'log_eddy_tracks',name,'.txt']);

%----------------------------------------------------------
% load eddy centers and eddy shapes

disp('Load centers and shapes ...')

load([path_out,'eddy_centers',name],'centers2');
load([path_out,'eddy_shapes',name]);

stepF = length(shapes1);

%----------------------------------------------------------
% Names list from centers and shapes structure

var_name1 = fieldnames(centers2);
var_name2 = fieldnames(shapes1);
var_name3 = fieldnames(shapes2);
var_name4 = {'f','Rd','gama','calcul_curve','large_curve1','large_curve2','too_weak2'};

%----------------------------------------------------------
% intitialize/update tracks structure arrays
% - tracks contains all the open and closed tracks which will be saved;

if update && exist([path_out,'eddy_tracks',name,'.mat'],'file')
    
% Build record tracks and active search
    
    load([path_out,'eddy_tracks',name])
    
    step0 = stepF - update + 1;

    tracks_name = fieldnames(tracks);
    warn_name = fieldnames(warn_tracks);
    
% Remove eddies newer than the updating step

    moved = false(1,length(tracks));

    for i=1:length(tracks)
        
        % find step earlier than update date
        ind = tracks(i).step < centers2(step0).step;
        
        % record earlier step
        for n=1:length(tracks_name)
            tracks(i).(tracks_name{n}) = tracks(i).(tracks_name{n})(ind,:);
        end
        for n=[3 8:9]
            warn_tracks(i).(warn_name{n}) = warn_tracks(i).(warn_name{n})(ind,:);
        end
    
        % find step earlier than update date
        ind = warn_tracks(i).step < step0;
        for n=[1:2 4:7 10:11]
            warn_tracks(i).(warn_name{n}) = warn_tracks(i).(warn_name{n})(ind,:);
        end
        
        moved(i) = isempty(tracks(i).step);
        
    end
    
    tracks(moved) = [];
    warn_tracks(moved) = [];
    
    disp(['  update eddy_tracks',name])
    
else
% preallocate tracks array

    tracks = struct('step',[],'type',[],'x1',[],'y1',[],...
                'x2',[],'y2',[],'dc',[],'ind2',[],...
                'shapes1',[],'velmax1',[],'deta1',[],'tau1',[],'nrho1',[],...
                'rmax1',[],'aire1',[],'xbary1',[],'ybary1',[],'ellip1',[],'theta1',[],...
                'shapes3',[],'velmax3',[],'deta3',[],'nrho3',[],'rmax3',[],'aire3',[]);
            
    % add streamlines profiles if needed
    if streamlines
        struct1 = struct('alpha',[],'rsquare',[],'rmse',[]);
        names1 = [fieldnames(tracks); fieldnames(struct1)];
        tracks = cell2struct([struct2cell(tracks); struct2cell(struct1)], names1, 1);
    end
    
    % add computed fields if needed
    if extended_diags==1
        struct1 = struct('ke1',[],'vort1',[],'vortM1',[],'OW1',[],'LNAM1',[]);
        names2 = [fieldnames(tracks); fieldnames(struct1)];
        tracks = cell2struct([struct2cell(tracks); struct2cell(struct1)], names2, 1);
    end
    
    % add double contour 
    struct1 = struct('shapes2',[],'velmax2',[],'deta2',[],'tau2',[],'nrho2',[],...
                'rmax2',[],'aire2',[],'xbary2',[],'ybary2',[],'ellip2',[],'theta2',[]);

    names2 = [fieldnames(tracks); fieldnames(struct1)];
    tracks = cell2struct([struct2cell(tracks); struct2cell(struct1)], names2, 1);
    
    % add warnings shapes and interaction
    struct1 = struct('f',[],'Rd',[],'gama',[],'calcul',[],'large1',[],'large2',[],...
        'weak',[],'ind',[],'interaction',[],'interaction2',[]);
    names3 = [fieldnames(tracks); fieldnames(struct1)];
    tracks = cell2struct([struct2cell(tracks); struct2cell(struct1)], names3, 1);
    
    % build warn_tracks
    warn_tracks = struct('step',[],'last',[],'ind',[],'can',[],...
                        'd',[],'D',[],'C',[],...
                        'interaction',[],'interaction2',[],'far1',[],'far2',[]);
    step0 = 1;

    tracks_name = fieldnames(tracks);
    warn_name = fieldnames(warn_tracks);

end

n1 = find(strcmp(tracks_name,'x2'));
n2 = find(strcmp(tracks_name,'shapes2'));

disp(' ')

%----------------------------------------------------------
% loop through all the steps in which eddies were detected

disp('Tracks eddies')

% /!\ quick fix for multi years tracking
% works only of daily steps !!!
% if stepF-step0 > 366
%     stepF = step0 + 370;
% end

for i=step0:stepF
    
    % variables containing data of all the eddies detected for the
    % current step
    
    disp(' ')
    
    stp = centers2(i).step;
    disp([' Searching step ',num2str(stp),' %------------- '])
    
    disp(' ')
    
    for n=2:length(var_name1)
        eddy.(tracks_name{n}) = centers2(i).(var_name1{n});
    end
    
    N = length(var_name1);
    for n=1:length(var_name2)-1
        eddy.(tracks_name{n+N}) = shapes1(i).(var_name2{n+1});
    end
    
    N = N + length(var_name2)-1;
    for n=1:length(var_name3)-1
        eddy.(tracks_name{n+N}) = shapes2(i).(var_name3{n+1});
    end
    
    N = N + length(var_name3)-1;
    for n=1:length(var_name4)
        eddy.(tracks_name{n+N}) = warn_shapes2(i).(var_name4{n});
    end
    
    %----------------------------------------------------------
    % allocate interacting eddies indice
    
    % initialisation
    eddy.ind = 1:length(eddy.type);
    eddy.interaction = eddy.ind2;
    eddy.interaction2 = nan(1,length(eddy.type));
    
    % allocate first interacting eddy
    Nind = find(~isnan(eddy.ind2));
    ind = eddy.ind2(Nind);
    eddy.interaction(ind) = Nind;

    % allocate case with 3 eddies interacting each others
    Nind2 = ~isnan(eddy.interaction) & eddy.interaction~=eddy.ind2;
    eddy.interaction2(Nind2) = eddy.ind2(Nind2);

    % allocate case with a second interacting eddy
    [u,Ia,~] = unique(ind,'first');
    if ~isscalar(u)
        n = logical(hist(ind,u)-1);
        eddy.interaction2(u(n)) = Nind(Ia(n));
    elseif ~isscalar(~isnan(ind))
        eddy.interaction2(u) = Nind(Ia);
    end

    % check every eddy gots less than 2 candidats
    Nind3 = find(~isnan(eddy.ind2) &...
        eddy.ind2 ~= eddy.interaction & ...
        eddy.ind2 ~= eddy.interaction2 & ...
        eddy.interaction ~= eddy.interaction2,1);

    if ~isempty(Nind3)
        disp(['problem remaining for eddy ',num2str(eddy.ind(Nind3)),' step ',num2str(i)])
    end
    
    %----------------------------------------------------------
    % Begin updating eddy tracks ------------------------------
    %----------------------------------------------------------
    if ~isempty(eddy.type) && ~isempty(eddy.velmax1)% be sure eddies were present that specific step
        
        %----------------------------------------------------------
        % if first time-step, then all the eddies are open tracks 
        % (added to tracks matrix)
        if i==1
            
            for i2=1:length(eddy.type)
                
                tracks(i2).step = stp;

                for n=2:length(tracks_name)
                    tracks(i2).(tracks_name{n}) = eddy.(tracks_name{n})(i2);
                end
                
                % initialize warn_tracks
                warn_tracks(i2).step = stp;
                warn_tracks(i2).last = 0;
                warn_tracks(i2).ind = eddy.ind(i2);
                warn_tracks(i2).can = zeros(1,N_can);
                warn_tracks(i2).d = inf(1,N_can);
                warn_tracks(i2).D = zeros(1,N_can);
                warn_tracks(i2).C = inf(1,N_can);
                warn_tracks(i2).interaction = eddy.interaction(i2);
                warn_tracks(i2).interaction2 = eddy.interaction(i2);
                warn_tracks(i2).far1 = NaN;
                warn_tracks(i2).far2 = NaN;

            end

            % display number of first eddies
            disp(['  -> ',num2str(length(tracks)),' total eddies'])
    
        %----------------------------------------------------------
        % if not first step, then open tracks from previous steps
        % are updated with eddies detected for the current step
        else
            
            %----------------------------------------------------------
            % eddy features at current time
            ind_new = nan(length(eddy.type),13);
            
            ind_new(:,1) = eddy.type; % type
            ind_new(:,2) = eddy.x1;
            ind_new(:,3) = eddy.y1;
            ind_new(:,4) = eddy.xbary1; % barycenter x
            ind_new(:,5) = eddy.ybary1; % barycenter y
            % default center when no barycenter found (no shapes)
            ind_nil = find(ind_new(:,4)==0 & ind_new(:,5)==0);
            ind_new(ind_nil,4) = eddy.x1(ind_nil);
            ind_new(ind_nil,5) = eddy.y1(ind_nil);
            
            if extended_diags==1
                    ind_new(:,13) = eddy.vortM1 * 1e6; % maximale vorticity
            end
            
            ind_new(:,6) = eddy.rmax1; % rmax
            ind_new(:,7) = eddy.velmax1 * T * 1e-3; % velmax in km/day
            ind_new(:,8) = eddy.rmax3; % rmax
            ind_new(:,9) = eddy.velmax3 * T * 1e-3; % velmax in km/day
            
            % find case no shapes1
            ind_nil = find(eddy.rmax1==0);
            
            % replace by shapes3
            ind_new(ind_nil,4) = eddy.x1(ind_nil);
            ind_new(ind_nil,5) = eddy.y1(ind_nil);
            ind_new(ind_nil,6) = eddy.rmax3(ind_nil);
            ind_new(ind_nil,7) = eddy.velmax3(ind_nil) * T * 1e-3;
            
            % Ro and Bu numbers
            ind_new(:,10) = ind_new(:,7) / T ./ eddy.f' ./ ind_new(:,6); % Ro = velmax/f./rmax*1e-3
            ind_new(:,12) = g * eddy.deta1 ./ (eddy.f .* eddy.rmax1).^2 * 1e-3; % Bu = g*deta./(f^2*rmax.^2)*1e-3
            
            % Cost matrice Ntracks in row X Mnew in column
            C = Inf(length(tracks),length(ind_new(:,1)));

            %----------------------------------------------------------
            % loop all tracks from the previous time steps
            parfor i2=1:length(tracks)
                
                %------------------------------------------------------
                % initialise 
                ind_old = nan(1,10);
                Ctot = C(i2,:);
                % time steps since every detection
                last = stp - tracks(i2).step;
                
                j = length(last);
                
                %------------------------------------------------------
                % 1st: find x and y of the eddy centers in tracks by 
                % scanning among the 'Dt' previous
                % days starting with latest step
                %------------------------------------------------------
                
                if last(j)*dps<=Dt
                    
                    %--------------------------------------------------
                    % period to define typical features of the last period
                    L_stp = max(1,j-D_stp):j;

                    % indice with (Y) and without (N) interaction with other eddy
                    Yint = find(~isnan(tracks(i2).interaction(L_stp)));
                    Nint = find(isnan(tracks(i2).interaction(L_stp)));
                    
                    % Indice of the other eddy
                    IYint = tracks(i2).interaction(L_stp(Yint));
                    SYint = tracks(i2).step(L_stp(Yint));
                    
                    % Record features from double eddy
                    Tint = struct('rmax2',[],'velmax2',[],'deta2',[]);
                    Tint.rmax2 = nan;
                    Tint.velmax2 = nan;
                    Tint.deta2 = nan;
                    Tint.f = nan;
                    
                    for k=1:length(SYint)
                        % is double eddy features in i2
                        if ~isnan(tracks(i2).rmax2(L_stp(Yint(k))))
                            Ik = L_stp(Yint(k));
                            Tint.rmax2(k,1) = tracks(i2).rmax2(Ik);
                            Tint.velmax2(k,1) = tracks(i2).velmax2(Ik);
                            Tint.deta2(k,1) = tracks(i2).deta2(Ik);
                            Tint.f(k,1) = tracks(i2).f(Ik);
                        % or in the other interacting eddy
                        else
                            Ik = find(tracks(IYint(k)).step==SYint(k));
                            Tint.rmax2(k,1) = tracks(IYint(k)).rmax2(Ik);
                            Tint.velmax2(k,1) = tracks(IYint(k)).velmax2(Ik);
                            Tint.deta2(k,1) = tracks(IYint(k)).deta2(Ik);
                            Tint.f(k,1) = tracks(IYint(k)).f(Ik);
                        end
                    end
                    
                    %--------------------------------------------------
                    % features of eddy i2 in previous time steps
                    ind_old(1) = tracks(i2).type(j);
                    ind_old(2) = tracks(i2).x1(j);
                    ind_old(3) = tracks(i2).y1(j);
                    ind_old(4) = tracks(i2).xbary1(j); % barycenter x
                    ind_old(5) = tracks(i2).ybary1(j); % barycenter y
                    if ind_old(4)==0 && ind_old(5)==0
                        ind_old(4) = tracks(i2).x1(j);
                        ind_old(5) = tracks(i2).y1(j);
                    end
                    
                    f = tracks(i2).f(j);
                    
                    if extended_diags==1
                            % maximale vorticity
                            ind_old(13) = nanmean(tracks(i2).vortM1(L_stp([Nint;Yint]))) * 1e6;
                    end

                    % rmax from the D_stp+1 previous detections
                    ind_old(6) = nanmean(tracks(i2).rmax1(L_stp([Nint;Yint])));
                    
                    % velmax from the D_stp+1 previous detections
                    ind_old(7) = nanmean(tracks(i2).velmax1(L_stp([Nint;Yint]))) * T * 1e-3;
                    
                    % rmax3 or rmax2 from the D_stp+1 previous detections
                    ind_old(8) = nanmean([tracks(i2).rmax3(L_stp(Nint));Tint.rmax2]);
                    
                    % velmax3 or rmax2 from the D_stp+1 previous detections
                    ind_old(9) = nanmean([tracks(i2).velmax3(L_stp(Nint));Tint.velmax2]) * T * 1e-3;
                    
                    % Ro number Ro = velmax/f./rmax*1e-3 from rmax1
                    ind_old(10) = nanmean(tracks(i2).velmax1(L_stp([Nint;Yint])) ./ ...
                                tracks(i2).rmax1(L_stp([Nint;Yint])) ./ ...
                                tracks(i2).f(L_stp([Nint;Yint]))) * 1e-3;
                            
                    % Ro number Ro = velmax/f./rmax*1e-3 from rmax2
                    ind_old(11) = nanmean(Tint.velmax2 ./ Tint.rmax2 ./ Tint.f) * 1e-3;
                            
                    % Bu_eta number Bu = g*deta./(f^2*rmax.^2)*1e-3 from rmax1 or rmax2
                    ind_old(12) = g * nanmean([tracks(i2).deta1(L_stp(Nint)) ./ ...
                                ( tracks(i2).rmax1(L_stp(Nint)) .* tracks(i2).f(L_stp(Nint)) ).^2;...
                                Tint.deta2 ./ (Tint.rmax2 .* Tint.f).^2]) * 1e-3;

                    %--------------------------------------------------
                    % 2nd: find centers after 'last' steps no more than
                    % 'Dt' days that are within 'D' from the tracked
                    % center then calculate the cost matrix C which is:
                    % C = sqrt ( dt/2Dt ² + <d1,d2>/Dmax ² + <dR1,dR2/3>/SR ² + <dRo1,dRo2>/SRo ² )
                    % with Dmax, SR and SRo from the last D_stp steps
                    %--------------------------------------------------

                    % - assumed maximum distance (D) is:
                    % last D_stp+1 radius + new radius
                    % + typical center deplacement
                    if V_eddy>0
                        D = V_eddy*dps*(1+last(j))/2 + ind_old(6) + ind_new(:,6);
                    else
                    % or + maximal velocity speed = Vmax/2
                        D = ind_old(7)/2*dps*(1+last(j))/2 + ind_old(6) + ind_new(:,6);
                    end
                    
                    % - overall maximum distance assuming Veddy of 6 km/d
                    Dmax = V_eddy*dps*(1+Dt)/2 + ind_old(6) + ind_new(:,6); 
                    
                    % - Distance between centers (d1) and barycenter (d2)
                    d1 = zeros(size(ind_new,1),1);
                    d2 = d1;
                    
                    for i3=1:size(d1)

                        if grid_ll
                            d1(i3) = sw_dist2([ind_new(i3,3) ind_old(3)],...
                                        [ind_new(i3,2) ind_old(2)],'km');
                            d2(i3) = sw_dist2([ind_new(i3,5) ind_old(5)],...
                                        [ind_new(i3,4) ind_old(4)],'km');
                        else
                            d1(i3) = sqrt((ind_new(i3,2) - ind_old(2)).^2 +...
                                    (ind_new(i3,3) - ind_old(3)).^2);
                            d2(i3) = sqrt((ind_new(i3,4) - ind_old(4)).^2 +...
                                    (ind_new(i3,5) - ind_old(5)).^2);
                        end
                        
                    end
                    
                    % - Cost for eddy of same type at less than 'D' km
                    %   center (and barycenter) to center (and barycenter)
                    ind = find( ind_new(:,1)==ind_old(1) & nanmean([d1 d2],2)<=D );
                    
                    if ~isempty(ind)
                        C1 = (last(j)-1)/ Dt / 2; % time absence cost
                        C2 = d1(ind) ./ Dmax(ind); % center distance cost
                        C3 = d2(ind) ./ Dmax(ind); % barycenter distance cost
                        C4 = (ind_new(ind,6) - ind_old(6)) ./ (ind_new(ind,6) + ind_old(6)); % size cost rmax1 contour
                        C5 = (ind_new(ind,8) - ind_old(8)) ./ (ind_new(ind,8) + ind_old(8)); % size cost rmax3 contour
                        C6 = (ind_new(ind,10) - ind_old(10)) ./ (ind_new(ind,10) + ind_old(10)); % intensity cost rmax1
                        C7 = (ind_new(ind,10) - ind_old(11)) ./ (ind_new(ind,10) + ind_old(11)); % intensity cost rmax2 if any
                        
                        Ctot(ind) = sqrt( C1^2 + min(nanmean([C2 C3],2),1).^2 +...
                            min(nanmean([C4 C5],2),1).^2 + min(nanmean([C6 C7],2),1).^2);
                        C(i2,:) = Ctot;
                    end
                    
                    %--------------------------------------------------
                    % Record case in which more than one center is found
                    % within the 'D' area for searching eddies
                    ind = find(Ctot<inf);
                    mov = length(ind);

                    % display and update warning
                    switch mov
                        case 0
                            disp(['  Tracked eddy ',num2str(i2),' (center ',...
                                    num2str(tracks(i2).ind(end)),')',...
                                    ' at stp-',num2str(last(j)),...
                                    ' temporaly lost'])
                        case 1
                        otherwise
                            disp(['  Tracked eddy ',num2str(i2),' (center ',...
                                    num2str(tracks(i2).ind(end)),')',...
                                    ' at stp-',num2str(last(j)),', '...
                                    num2str(mov),' possibilities: ',...
                                    num2str(ind)])
                    end
                    
                    %--------------------------------------------------
                    % Record tracking choices
                    warn_tracks(i2).step = [warn_tracks(i2).step; stp];
                    warn_tracks(i2).last = [warn_tracks(i2).last; last(j)];
                    warn_tracks(i2).can = [warn_tracks(i2).can; [ind zeros(1,N_can-mov)]];
                    warn_tracks(i2).d = [warn_tracks(i2).d; [d1(ind)' inf(1,N_can-mov)]];
                    warn_tracks(i2).D = [warn_tracks(i2).D; [D(ind)' zeros(1,N_can-mov)]];
                    warn_tracks(i2).C = [warn_tracks(i2).C; [Ctot(ind) inf(1,N_can-mov)]];
                    warn_tracks(i2).far1 = [warn_tracks(i2).far1; NaN];
                    warn_tracks(i2).far2 = [warn_tracks(i2).far2; NaN];
            
                end % record too old
                
            end % loop on tracks
            
            %--------------------------------------------------
            % Record case in which more than one center is found
            % within the 'D' area for new eddies at step stp
            disp(' ')
            
            for i3=1:size(C,2)
                
                ind = find(C(:,i3)<inf)';
                mov = length(ind);
                
                % display and update warning
                switch mov
                    case 0
                        disp(['  New eddy ',num2str(i3),' at t',...
                                ' just appears'])
                    case 1
                    otherwise
                        disp(['  New eddy ',num2str(i3),' at t, ',...
                                num2str(mov),' possibilities: ',...
                                num2str(ind)])
                end
                
            end
            
            disp(' ')
            
            %----------------------------------------------------------
            % 3rd: Assignment cost matrice C resolved by Munkres method
            % (Markus Buehren routine). The row index of C corresponds to
            % the old eddies and the column to the new. The result gives
            % the optimal assigment of new to old eddies in case of
            % multiple assigment possible by minimizing the cost function.
            %----------------------------------------------------------
            
            [assign,cost] = assignmentoptimal(C);
            
            if cost==0 && sum(assign)~=0
                display ('Could not find solution')
                display (['Check the computation of Cost matrix ''C'' ',...
                    'with the ''assignmentoptimal.m'' routine'])
                return
            end
                
            %----------------------------------------------------------
            % 4th: update tracks with centers from current step
            %----------------------------------------------------------
            
            % convert eddy indice to tracks indice
            N = length(tracks);
            Nind = zeros(length(eddy.type),1);
            
            % convert old eddy
            Nass = find(assign~=0);
            Nind(assign(Nass)) = Nass;
            added = Nind==0;

            % convert new eddy
            Sind = cumsum(added);
            Nind(added) = N + Sind(added);
            
            %--------------------------------------------------
            % update tracks if any
            for i2=1:length(tracks)
                
                old = assign(i2);
                
                if old~=0

                    % add info of the new eddy center to the open 
                    % track array
                    tracks(i2).step = [tracks(i2).step; stp];
                    
                    for n=2:length(tracks_name)-2
                        tracks(i2).(tracks_name{n}) =...
                            [tracks(i2).(tracks_name{n});...
                            eddy.(tracks_name{n})(old)];
                    end
                    
                    % allocate tracks indice to interaction
                    if ~isnan(eddy.interaction(old))
                        tracks(i2).interaction = [tracks(i2).interaction; ...
                            Nind(eddy.interaction(old))];
                    else
                        tracks(i2).interaction = [tracks(i2).interaction; NaN];
                    end
                    if ~isnan(eddy.interaction2(old))
                        tracks(i2).interaction2 = [tracks(i2).interaction2; ...
                            Nind(eddy.interaction2(old))];
                    else
                        tracks(i2).interaction2 = [tracks(i2).interaction2; NaN];
                    end
                    
                    % record the assigned and interacting eddies indices
                    warn_tracks(i2).ind = [warn_tracks(i2).ind; eddy.ind(old)];
                    warn_tracks(i2).interaction = [warn_tracks(i2).interaction; ...
                        eddy.interaction(old)];
                    warn_tracks(i2).interaction2 = [warn_tracks(i2).interaction2; ...
                        eddy.interaction2(old)];
                    
                end
            end
            
            %----------------------------------------------------------
            % 5th: remaining eddies from present step are new tracks
            % and added to tracks array
            %----------------------------------------------------------
            
            % extract remaining eddies
            for n=2:length(tracks_name)
                eddy.(tracks_name{n}) =...
                    eddy.(tracks_name{n})(added);
            end
            
            if ~isempty(eddy.type)
                for i3=1:length(eddy.type)
                    
                    % start new tracks
                    tracks(length(tracks)+1).step = stp;  
                    
                    % record new tracks
                    for n=2:length(tracks_name)-2
                        tracks(length(tracks)).(tracks_name{n}) =...
                            eddy.(tracks_name{n})(i3);
                    end
                    
                    % allocate tracks indice to interaction
                    if ~isnan(eddy.interaction(i3))
                        tracks(length(tracks)).interaction = Nind(eddy.interaction(i3));
                    else
                        tracks(length(tracks)).interaction = NaN;
                    end
                    if ~isnan(eddy.interaction2(i3))
                        tracks(length(tracks)).interaction2 = Nind(eddy.interaction2(i3));
                    else
                        tracks(length(tracks)).interaction2 = NaN;
                    end
                    
                    % record the assigned and interacting eddies indices
                    warn_tracks(length(warn_tracks)+1).step = stp;
                    warn_tracks(length(warn_tracks)).last = 0;
                    warn_tracks(length(warn_tracks)).ind = eddy.ind(i3);
                    warn_tracks(length(warn_tracks)).can = zeros(1,N_can);
                    warn_tracks(length(warn_tracks)).d = inf(1,N_can);
                    warn_tracks(length(warn_tracks)).D = zeros(1,N_can);
                    warn_tracks(length(warn_tracks)).C = inf(1,N_can);
                    warn_tracks(length(warn_tracks)).interaction = eddy.interaction(i3);
                    warn_tracks(length(warn_tracks)).interaction2 = eddy.interaction2(i3);
                    warn_tracks(length(warn_tracks)).far1 = NaN;
                    warn_tracks(length(warn_tracks)).far2 = NaN;
                    
                end
            end
            
            % display number of new eddies
            disp([' Step ',num2str(stp)])
            disp(['  -> ',num2str(sum(added)),' new eddies'])
            
            % clear some variables (re-initialized at next loop)
            clear ind_new
            
        end % end if first time step or not
        
        %----------------------------------------------------------
        % 6th: eddies are considered dissipated when tracks are not
        % updated for longer than 'Dt' days.
        %----------------------------------------------------------

        moved = 0;

        for i2=1:length(tracks)

            % find interacting tracks
            ind1 = tracks(i2).interaction(end);
            ind2 = tracks(i2).interaction2(end);

            if tracks(i2).step(end) <= stp - Dt/dps &&...
                    tracks(i2).step(end) > stp - Dt/dps - 1

                % create array to change search indice to tracks indice
                moved = moved + 1;

            %----------------------------------------------------------
            % 7th: Remove interaction between eddies too far from each other
            % if 2 shapes1 exist
            %   remove interaction with dc > 3.5 * (rmax1 + rmax2)/2
            % else
            %   remove interaction if dc > 3.5 * rmax
            % end
            %----------------------------------------------------------
            elseif tracks(i2).step(end) == stp

                % for every interaction
                if ~isnan(ind1)
                    
                    j = length(tracks(i2).rmax1);
                    j1= length(tracks(ind1).rmax1);

                    % mmean radius
                    r = nanmean(tracks(i2).rmax1(max(1,j-D_stp):end));
                    r1 = nanmean(tracks(ind1).rmax1(max(1,j1-D_stp):end));

                    % disance between 2 centers
                    if grid_ll
                        d1 = sw_dist2([tracks(i2).y1(end) tracks(ind1).y1(end)],...
                                    [tracks(i2).x1(end) tracks(ind1).x1(end)],'km');
                    else
                        d1 = sqrt((tracks(i2).x1(end) - tracks(ind1).x1(end)).^2 +...
                                (tracks(i2).y1(end) - tracks(ind1).y1(end)).^2);
                    end

                    % remove too far interaction
                    if ( ~isnan(r) && ~isnan(r1) && d1 > (dc_max+1) /2 * (r+r1) ) ||...
                        ( ~isnan(r) && isnan(r1) && d1 > (dc_max+1) * (r) ) ||...
                        ( isnan(r) && ~isnan(r1) && d1 > (dc_max+1) * (r1) )
                        
                        tracks(i2).interaction(end) = NaN;
                        warn_tracks(i2).far1(end) = 1;

                        if isnan(tracks(i2).interaction2(end))
                            tracks(i2).shapes2{end} = NaN;
                            for n=[n1:n1+3 n2+1:n2+10]
                                tracks(i2).(tracks_name{n})(end) = NaN;
                            end
                        end

                        if tracks(ind1).interaction(end)~=i2 && ...
                                tracks(ind1).interaction2(end)~=i2
                            disp(['       remove interaction between ',...
                                num2str(i2),' and ',num2str(ind1)])
                        end

                    else
                        warn_tracks(i2).far1(end) = 0;
                        if tracks(ind1).interaction(end) == i2
                            warn_tracks(ind1).far1(end) = 0;
                        else
                            warn_tracks(ind1).far2(end) = 0;
                        end
                    end
                end

                % for every interaction2
                if ~isnan(ind2)

                    % mean radius
                    r = nanmean(tracks(i2).rmax1(max(1,end-D_stp):end));
                    r2 = nanmean(tracks(ind2).rmax1(max(1,end-D_stp):end));

                    % disance between 2 centers
                    if grid_ll
                        d2 = sw_dist2([tracks(i2).y1(end) tracks(ind2).y1(end)],...
                                    [tracks(i2).x1(end) tracks(ind2).x1(end)],'km');
                    else
                        d2 = sqrt((tracks(i2).x1(end) - tracks(ind2).x1(end)).^2 +...
                                (tracks(i2).y1(end) - tracks(ind2).y1(end)).^2);
                    end

                    % remove too far interaction2
                    if ( ~isnan(r) && ~isnan(r2) && d2 > (dc_max+1) /2 * (r+r2) ) ||...
                        ( ~isnan(r) && isnan(r2) && d2 > (dc_max+1) * (r) ) ||...
                        ( isnan(r) && ~isnan(r2) && d2 > (dc_max+1) * (r2) )
                        
                        tracks(i2).interaction2(end) = NaN;
                        warn_tracks(i2).far2(end) = 1;

                        if isnan(tracks(i2).interaction(end))
                            tracks(i2).shapes2{end} = NaN;
                            for n=[n1:n1+3 n2+1:n2+10]
                                tracks(i2).(tracks_name{n})(end) = NaN;
                            end
                        end

                        if tracks(ind2).interaction(end)~=i2 && ...
                                tracks(ind2).interaction2(end)~=i2
                            disp(['       remove interaction between ',...
                                num2str(i2),' and ',num2str(ind2)])
                        end

                    else
                        warn_tracks(i2).far2(end) = 0;
                        if tracks(ind2).interaction(end) == i2
                            warn_tracks(ind2).far1(end) = 0;
                        else
                            warn_tracks(ind2).far2(end) = 0;
                        end
                    end
                end

            end

        end

        % display number of old eddies
        disp(['  -> ',num2str(moved),' too old eddies'])

    end % end if eddy present at that specific step
    
end % end updating eddy tracks for all step

disp(' ')

%----------------------------------------------------------
% save tracks and warnings in structure array

disp(' ')

disp(['Save tracks in eddy_tracks',name,'.mat ...'])

save([path_out,'eddy_tracks',name],'tracks','warn_tracks','-v7.3')

% close log file
diary off

disp(' ')
