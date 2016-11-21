function step0 = mod_init(stepF,update)
%step0 = mod_init({stepF,update})
%
%   mod_init preallocate structure or update the structure in mat file
%
%   The routine use:
%   - update is a flag allowing to update an existing tracking:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the time serie
%   - stepF: last time step of the computed series since the initial step
%       (be careful to update only serie from the same initial step.)
%
%   The results give the first time step to be computed (step0).
%       (for update == 0, stepO = 1)
%
%-------------------------
%   July 2016 Briac Le Vu
%-------------------------
%
%=========================

% Load default parameters
%----------------------------------------
if nargin==0
    load('param_eddy_tracking','path_out','streamlines','extended_diags','stepF')
    update = 0;
elseif nargin==1
    load('param_eddy_tracking','path_out','streamlines','extended_diags')
    update = 0;
else
    load('param_eddy_tracking','path_out','streamlines','extended_diags')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Updating structure ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if update
    
    step0 = stepF - update+1;
    
    % Load struct array
    %----------------------------------------
    load([path_out,'fields']);
    detection_fields_ni = detection_fields;
    load([path_out,'fields_inter']);
    load([path_out,'eddy_centers']);
    load([path_out,'eddy_shapes'])

    % Keep steps preexisting update
    %----------------------------------------
    detection_fields_ni = detection_fields_ni(1:step0-1);
    detection_fields = detection_fields(1:step0-1);
    centers0 = centers0(1:step0-1);
    centers  = centers(1:step0-1);
    centers2 = centers2(1:step0-1);
    
    shapes1      = shapes1(1:step0-1);
    shapes2      = shapes2(1:step0-1);
    warn_shapes  = warn_shapes(1:step0-1);
    warn_shapes2 = warn_shapes2(1:step0-1);
    
    if streamlines
        profil2  = profil2(1:step0-1);
    end

    % Preallocate steps next to update
    %----------------------------------------
    names = fieldnames(detection_fields);
    for n=1:length(names)
        detection_fields_ni(stepF).(names{n}) = [];
        detection_fields(stepF).(names{n}) = [];
    end

    names = fieldnames(centers);
    for n=1:length(names)
        centers0(stepF).(names{n}) = [];
        centers(stepF).(names{n}) = [];
    end
    
    names = fieldnames(centers2);
    for n=1:length(names)
        centers2(stepF).(names{n}) = [];
    end

    names = fieldnames(shapes1);
    for n=1:length(names)
        shapes1(stepF).(names{n}) = [];
    end

    names = fieldnames(shapes2);
    for n=1:length(names)
        shapes2(stepF).(names{n}) = [];
    end

    names = fieldnames(warn_shapes);
    for n=1:length(names)
        warn_shapes(stepF).(names{n}) = [];
        warn_shapes2(stepF).(names{n}) = [];
    end

    if streamlines
        names = fieldnames(profil2);
        for n=1:length(names)
            profil2(stepF).(names{n}) = [];
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preallocate structure ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else

    step0 = 1;
    
    % Preallocate fields struct array
    %----------------------------------------
    detection_fields_ni(stepF) = struct('step',[],'ke',[],'div',[],...
        'vort',[],'OW',[],'LOW',[],'LNAM',[]);
    detection_fields = detection_fields_ni;
    
    % Preallocate centers struct array
    %----------------------------------------
    centers0(stepF) = struct('step',[],'type',[],'x',[],'y',[],'i',[],'j',[]);
    centers = centers0;
    centers2(stepF) = struct('step',[],'type',[],'x1',[],'y1',[],...
        'x2',[],'y2',[],'dc',[],'ind2',[]);
    
    % Preallocate shapesstruct array
    %----------------------------------------
    shapes1(stepF) = struct('step',[],'xy',[],'velmax',[],'taumin',[],...
        'deta',[],'nrho',[],'rmax',[],'aire',[],...
        'xy_end',[],'vel_end',[],'deta_end',[],'r_end',[],'aire_end',[]);
    shapes2(stepF) = struct('step',[],'xy',[],'velmax',[],...
        'deta',[],'nrho',[],'rmax',[],'aire',[]);
    
    if streamlines
        profil2(stepF) = struct('step',[],'nc',[],'eta',[],'rmoy',[],'vel',[],'tau',[],'myfit',[]);
        struct1(stepF) = struct('alpha',[],'rsquare',[],'rmse',[]);
        names1 = [fieldnames(shapes1); fieldnames(struct1)];
        shapes1 = cell2struct([struct2cell(shapes1); struct2cell(struct1)], names1, 1);
        
    end
    
    if extended_diags==1
        struct2(stepF) = struct('xbary',[],'ybary',[],'ellip',[],...
                    'ke',[],'vort',[],'vortM',[],'OW',[],'LNAM',[]);
        struct3(stepF) = struct('xbary',[],'ybary',[],'ellip',[]);
        names2 = [fieldnames(shapes1); fieldnames(struct2)];
        names3 = [fieldnames(shapes2); fieldnames(struct3)];
        shapes1 = cell2struct([struct2cell(shapes1); struct2cell(struct2)], names2, 1);
        shapes2 = cell2struct([struct2cell(shapes2); struct2cell(struct3)], names3, 1);
    end
    
    % Preallocate warnings struct array
    %----------------------------------------
    warn_shapes(stepF) = struct('no_curve',[],'Rd',[],'gama',[],'bx',[],'calcul_curve',[],...
                    'large_curve1',[],'large_curve2',[],'too_weak2',[]);
    warn_shapes2 = warn_shapes;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save structure array ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save interpolated fields
%----------------------------------------
save([path_out,'fields_inter'],'detection_fields','-v7.3')

% Save non interpolated fields
%----------------------------------------
detection_fields = detection_fields_ni;
save([path_out,'fields'],'detection_fields','-v7.3')

% Save centers
%----------------------------------------
save([path_out,'eddy_centers'],'centers0','centers','centers2','-v7.3')

% save warnings, shapes
%----------------------------------------
if streamlines
    save([path_out,'eddy_shapes'],'shapes1','shapes2',...
        'warn_shapes','warn_shapes2','profil2','-v7.3')
else
    save([path_out,'eddy_shapes'],'shapes1','shapes2',...
        'warn_shapes','warn_shapes2','-v7.3')
end








