%MAIN_AMEDA
%
%   MAIN_AMEDA is the main function of the eddy detection and
%   tracking package. It returns position of the centers, dimensions and 
%   tracks of the eddies detected from the time series of a 2-D velocity 
%   field.
%   It gives also an history of the splitting and merging events.
%
%   - 'source' allows to specify the type of sources file (AVISO, ROMS, NEMO,...)
%     with their specific parameters and Input/Output.
%   - cpus to use 'parfor' as time loops (=# of processors)
%       cpus = 1 (default)
%   - update is a flag allowing to update an existing tracking:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the time serie
%   - stepF is the last time step computed
%       stepF = temporal size of the input data
%
%   The algortihm subroutines:
%
%   - mod_eddy_params sets user defined paths and parameters:
%     nc_u nc_v nc_dim b bx r path_in path_out periodic criteres
%     Users should modify keys_sources.m according to their 
%     settings.
%
%   - mod_init initialise or update mat-file.
%
%   - mod_fields compute LNAM.
%
%   - mod_eddy_centers returns a structure array with the position of the
%     detected eddy centers.
%
%   - mod_eddy_shapes computes dimensions for the detected eddy centers.
%
%   - mod_eddy_tracks computes eddy tracks using the detected centers.
%
%   Find the output files in path_out:
%
%   - fields.mat contains detection_fields with LNAM for each step.
%   - eddy_centers.mat contains for each step:
%       * centers0 as the local max(LNAM)
%       * centers as the potential centers
%       * centers2 as the detected eddies
%   - eddy_shapes.mat contains for each step:
%       * shapes1 the eddy features
%       * shapes2 the common double contour features
%       * profil2 the streamlines features scanned around each eddy
%       * warn_shapes the flag for potential centers
%       * warn_shapes2 the flag for detected eddies
%   - eddy_tracks.mat contains eddy centers, features and flags for each eddy
%
%-------------------------
%   June 2016 Briac Le Vu
%-------------------------
%
%=========================

run(['/home6/datahome/fdumas/MATLAB/start_datarmor'])
clear; clc;
disp(version)
ver

%----------------------------------------
% source of data driving the netcdf format
source = 'CROCO';
%source = 'AVISO';

%----------------------------------------
% domaine
keys = 'test_deg2';

%----------------------------------------
% Update option
update = 0; % the serie from the beginning

%----------------------------------------
% Possibility to shorter the serie
%stepF = 10;

%----------------------------------------
% Set parallel computation
cpus=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Produce default parameters in param_eddy_tracking
if exist('stepF','var')
    mod_eddy_params(['keys_sources_',source,'_',keys],stepF)
else
    mod_eddy_params(['keys_sources_',source,'_',keys])
end
run(['keys_sources_',source,'_',keys])
load('param_eddy_tracking','path_out','streamlines','resol','stepF');

%----------------------------------------
% Preallocate structure array and mat-file or prepare update
% !! replace or reinitialise previous results !!
step0 = mod_init(stepF,update);

%----------------------------------------
% Activate matlab pool
cpus=min([cpus,32]);%maximum of 24 procs

if cpus>1
    disp('Check that you have access to "Parallel Computing Toolbox" to use PARPOOL')
    disp('otherwise use MAIN_AMEDA_nopool')
    disp(' ')

    myCluster = parcluster('local');
    delete(myCluster.Jobs)
    mypool = parpool(cpus);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute LNAM ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' === Compute LNAM ==='])
disp(' ')

load([path_out,'fields'],'detection_fields')
detection_fields_ni = detection_fields;

load([path_out,'fields_inter.mat'],'detection_fields')

parfor stp = step0:stepF
    %----------------------------------------
    % Compute non interpolated fields for step stp
    detection_fields_ni(stp) = mod_fields(source,stp,1);
    if resol>1
        %----------------------------------------
        % Compute interpolated fields for step stp
        detection_fields(stp) = mod_fields(source,stp,resol);
    else
        %----------------------------------------
        % Interpolated and non interpolated field are the same
        %disp(' === Interpolated LNAM is the same ===')
        detection_fields(stp) = detection_fields_ni(stp);
    end
end

%----------------------------------------
% Save fields
save([path_out,'fields_inter'],'detection_fields','-v7.3')

if resol>1
    detection_fields = detection_fields_ni;
    save([path_out,'fields'],'detection_fields','-v7.3')
end

clear detection_fields detection_fields_ni

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find centers ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' === Find potential centers ==='])
disp(' ')

load([path_out,'eddy_centers'])

%----------------------------------------
% Build I/O matfile
fields_mat = matfile([path_out,'fields_inter.mat']);

parfor stp = step0:stepF
    % load inter fields at step stp
    %----------------------------------------
    fields = fields_mat.detection_fields(:,stp);
    %----------------------------------------
    % Detection of eddy centers for step stp
    [centers0(stp),centers(stp)] = mod_eddy_centers(source,stp,fields);
end

%----------------------------------------
% Save centers
save([path_out,'eddy_centers'],'centers0','centers','centers2','-v7.3')
clear centers0 centers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' === Determine eddies shapes ==='])
disp(' ')

load([path_out,'eddy_centers'],'centers2')
load([path_out,'eddy_shapes'])

%----------------------------------------
% Build I/O matfile
fields_mat = matfile([path_out,'fields_inter.mat']);
centers_mat = matfile([path_out,'eddy_centers.mat']);

parfor stp = step0:stepF
    %----------------------------------------
    % load fields at step stp
    fields = fields_mat.detection_fields(:,stp);
    %----------------------------------------
    % load potential centers at step stp
    centers = centers_mat.centers(:,stp);
    %----------------------------------------
    % Determination of eddy features for step stp
    if streamlines
        [centers2(stp),shapes1(stp),shapes2(stp),profil2(stp),...
        warn_shapes(stp),warn_shapes2(stp)] = ...
        mod_eddy_shapes(source,stp,fields,centers);
    else
        [centers2(stp),shapes1(stp),shapes2(stp),~,...
        warn_shapes(stp),warn_shapes2(stp)] = ...
        mod_eddy_shapes(source,stp,fields,centers);
    end
    
end

%----------------------------------------
% save warnings, shapes and their centers
save([path_out,'eddy_centers'],'centers2','-append')
if streamlines
    save([path_out,'eddy_shapes'],'shapes1','shapes2',...
        'warn_shapes','warn_shapes2','profil2','-v7.3')
else
    save([path_out,'eddy_shapes'],'shapes1','shapes2',...
        'warn_shapes','warn_shapes2','-v7.3')
end
clear centers2 shapes1 shapes2 profil2 warn_shapes warn_shapes2 struct1 struct2 struct3

%----------------------------------------
% Free workers
delete(gcp('nocreate'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Tracking eddies and record interacting events
name=[''];
mod_eddy_tracks_nopool(name,update);

%----------------------------------------
% Resolve merging and spltting event and filter eddies shorter than cut_off
% need a series longer than 2 turn over time (> 1 month)
mod_merging_splitting(name);





