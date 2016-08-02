function MAIN_AMEDA(source,runname,parallel,update,stepE)
% MAIN AMEDA routine to loop on time with PARFOR
%   MAIN_AMEDA is the main function of the eddy detection and
%   tracking package. It returns position of the centers, dimensions and 
%   tracks of the eddies detected from the time series of a 2-D velocity 
%   field.
%   It gives also an history of the splitting and merging events.
%
%   - 'source' allows to specify the type of sources file (AVISO, ROMS, NEMO,...)
%     with their specific parameters and Input/Output.
%   - runname is the prefix name of the output.
%   - parellel to use 'parfor' as time loops (=# of procs)
%   - update is a flag allowing to update an existing tracking:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the time serie
%   - stepE is the last time step computed (=stepF by default)
%
%   The algortihm subroutines:
%
%   - mod_eddy_params sets user defined paths and parameters:
%     nc_u nc_v nc_dim b bx r path_in path_out periodic criteres
%     Users should modify keys_sources.m according to their 
%     settings.
%
%   - mod_eddy_centers returns a structure array with the position of the
%     detected eddy centers.
%
%   - mod_eddy_shapes computes dimensions for the detected eddy centers.
%
%   - mod_eddy_tracks computes eddy tracks using the detected centers.
%
%-------------------------
%   June 2016 Briac Le Vu
%-------------------------
%
%=========================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load paths and keys ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start
clear; clc;

% source of data driving the netcdf format
%----------------------------------------
source = 'AVISO';

% use to submit different job in the same out directory at the same time
%----------------------------------------
runname = '_test'; % ex: _d1_d100 or ''

% Produce default parameters in param_eddy_tracking
%----------------------------------------
mod_eddy_params(['keys_sources_',source])
load('param_eddy_tracking')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialiation ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No update by default and ond possibility to change the last time step
if nargin<4
    update = 0;
elseif nargin>4
    stepF=stepE;
end

% Set parallel computation
if parallel>2
    myCluster = parcluster('local');
    delete(myCluster.Jobs)
    matlabpool open parallel
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute LNAM ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate fields struct array if doesn't exist
%----------------------------------------
if update
    step0 = stepF - update+1;
    load([path_out,'fields',runname])
    detection_fields = detection_fields(1:step0-1);
else
    step0 = 1;
    detection_fields(stepF) = struct('step',[],'ke',[],'div',[],'vort',[],'OW',[],'LOW',[],'LNAM',[]);
end

% Compute non interpolated fields step by step
%----------------------------------------
disp(['Compute non interpolated LNAM from step ',num2str(step0),' to step ',num2str(stepF)])
disp('"enlarge coastal mask" by adding b pixels of ocean to the coast and NO INTERPOLATION')
disp(' ')

parfor stp = step0:stepF
%for stp = step0:stepF
    detection_fields(stp) = mod_fields(source,stp,1);
end

% Save non interpolated fields
%----------------------------------------
save([path_out,'fields',runname],'detection_fields','-v7.3')
clear detection_fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute interpolated LNAM ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolated fields may be needed to accurately find centers
%----------------------------------------
if resol>1
    
    % Preallocate fields struct array if doesn't exist
    %----------------------------------------
    if update
        step0 = stepF - update+1;
        load([path_out,'fields_inter',runname])
        detection_fields = detection_fields(1:step0-1);
    else
        step0 = 1;
        detection_fields(stepF) = struct('step',[],'LOW',[],'LNAM',[]);
    end
    
    % Compute interpolated fields step by step
    %----------------------------------------
    disp(['Compute interpolated LNAM from step ',num2str(step0),' to step ',num2str(stepF)])
    disp(['"change resolution" by computing SPLINE INTERPOLATION res=',num2str(resol)])
    disp(' ')
        
    parfor stp = step0:stepF
    %for stp = step0:stepF
        detection_fields(stp) = mod_fields(source,stp,resol);
    end

    % Save interpolated fields
    %----------------------------------------
    save([path_out,'fields_inter',runname],'detection_fields','-v7.3')
    clear detection_fields
    
else

    % Interpolated and non interpolated field are the same
    %----------------------------------------
    disp('Copy non interpolated LNAM')
    disp(' ')
    
    copyfile([path_out,'fields',runname,'.mat'],[path_out,'fields_inter',runname,'.mat'])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find centers ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate fields struct array if doesn't exist
%----------------------------------------
if update && exist([path_out,'eddy_centers',runname,'.mat'],'file')
    step0 = stepF - update+1;
    load([path_out,'eddy_centers',runname]);
    
    if centers0(step0-1).step ~= step0-1
        display('Gap with the last recorded step!')
        return
    end
    
    centers0 = centers0(1:step0-1);
    centers  = centers(1:step0-1);
else
    step0 = 1;
    centers0(stepF) = struct('step',[],'type',[],'x',[],'y',[],'i',[],'j',[]);
    centers = centers0;
end

% Build io matfile
%----------------------------------------
fields_inter_mat = matfile([path_out,'fields_inter',runname,'.mat']);
fields_mat = matfile([path_out,'fields',runname,'.mat']);

% Detection of eddy centers step by step
%----------------------------------------
disp(['Find potential centers from step ',num2str(step0),' to step ',num2str(stepF)])
if resol==1
    disp('Use non interpolated fields')
else
    disp(['Use inteprolated fields resol=',num2str(resol)])
end
disp(' ')

parfor stp = step0:stepF
%for stp = step0:stepF
    % load inter fields at step stp
    detection_fields = fields_inter_mat.detection_fields(:,stp);
    % find max LNAM and potential centers
    [centers0(stp),centers(stp)] = mod_eddy_centers(source,stp,detection_fields);
end

% Save centers in struct array
%----------------------------------------
if update && exist([path_out,'eddy_centers',runname,'.mat'],'file')
    save([path_out,'eddy_centers',runname],'centers0','centers','-append')
else
    save([path_out,'eddy_centers',runname],'centers0','centers','-v7.3')
end
clear centers0 centers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin the log file
%----------------------------------------
diary([path_out,'log_eddy_shapes',runname,'.txt']);

% preallocate shape and warning array
%----------------------------------------------
if update && exist([path_out,'eddy_shapes',runname,'.mat'],'file')
    step0 = stepF - update + 1;
    load([path_out,'eddy_centers',runname])
    load([path_out,'eddy_shapes',runname])
    load([path_out,'warnings_shapes',runname])
    if centers2(step0-1).step ~= step0-1
        display('Gap with the last recorded step!')
        return
    end
    centers2     = centers2(1:step0-1);
    shapes1      = shapes1(1:step0-1);
    shapes2      = shapes2(1:step0-1);
    warn_shapes  = warn_shapes(1:step0-1);
    warn_shapes2 = warn_shapes2(1:step0-1);
    if streamlines
        profil2  = profil2(1:step0-1);
    end
else
    step0 = 1;
    centers2(stepF) = struct('step',[],'type',[],'x1',[],'y1',[],'x2',[],'y2',[],'dc',[],'ds',[]);
    shapes1(stepF) = struct('step',[],'xy',[],'velmax',[],'taumin',[],'deta',[],'rmax',[],'aire',[],...
                    'xy_end',[],'vel_end',[],'deta_end',[],'r_end',[],'aire_end',[]);
    shapes2(stepF) = struct('step',[],'xy',[],'velmax',[],'deta',[],'rmax',[],'aire',[]);
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
    % shapes struct for the possible second shape with 2 centers
    warn_shapes(stepF) = struct('no_curve',[],'fac',[],'bx',[],'calcul_curve',[],...
                    'large_curve1',[],'large_curve2',[],'too_large2',[]);
    warn_shapes2 = warn_shapes;
end

% Build io matfile
%----------------------------------------
centers_mat = matfile([path_out,'eddy_centers',runname,'.mat']);

% Detection of eddy shapes step by step
%----------------------------------------
disp(['Determine contour shapes on ',num2str(bx),'X',num2str(bx),' grid for ',runname])
disp(['from step ',num2str(step0),' to step ',num2str(stepF)])
disp('Use non interpolated fields')
disp(' ')

parfor stp = step0:stepF
%for stp = 1:stepF
    % load fields at step stp
    detection_fields = fields_mat.detection_fields(:,stp);
    % load potential centers at step stp
    potential_centers = centers_mat.centers(:,stp);
    % find eddy shapes
    [centers2(stp),shapes1(stp),shapes2(stp),profil2(stp),warn_shapes(stp),warn_shapes2(stp)] = ...
        mod_eddy_shapes(source,stp,detection_fields,potential_centers);
end

% save warnings, shapes and their centers in structure array
%----------------------------------------
save([path_out,'eddy_centers',runname],'centers2','-append')
if streamlines
    save([path_out,'eddy_shapes',runname],'shapes1','shapes2','warn_shapes','warn_shapes2','profil2')
else
    save([path_out,'eddy_shapes',runname],'shapes1','shapes2','warn_shapes','warn_shapes2')
end
clear centers2 shapes1 shapes2 profil2 warn_shapes warn_shapes2 struct1 struct2 struct3

% close log file
%----------------------------------------
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Track eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Concatenate multi years (use "runname" to loop on years)
%----------------------------------------
name={'2013','2014','2015','2016_apr'};
concat_eddy(name)

% Process tracking at once
%----------------------------------------
runname= ['_',char(name(1)),'_',char(name(end))];
%Tracking of the eddies 2013 to 2016
cut_off=1;% save minimal duration (0=tau(n))
Dt=5;% searching delay tolerance
mod_eddy_tracks(runname,cut_off,Dt,update)








