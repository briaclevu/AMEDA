function MAIN_AMEDA(source,cpus,update,stepF)
%MAIN_AMEDA(source {,cpus,update,stepF})
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load paths and keys ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% Possibility to shorter the serie
stepF = 2;

%----------------------------------------
% Default arguments
if nargin <2
    cpus = 1; % no parallel
    update = 0; % the entire serie
elseif nargin<3
    update = 0;
end

%----------------------------------------
% Produce default parameters in param_eddy_tracking
if exist('stepF','var')
    mod_eddy_params(['keys_sources_',source],stepF)
else
    mod_eddy_params(['keys_sources_',source])
end
load('param_eddy_tracking','path_out','resol','stepF')

%----------------------------------------
% Set parallel computation

%maximum of 12 procs
cpus=min([cpus,stepF-step0+1,4]);

if cpus>1
    
    disp(['Check you have "Parallel Computing Toolbox" to use ',...
        num2str(cpus),' processors'])
    disp(' ')
    
    delete(gcp)
    parpool('local',cpus)
    
end

%----------------------------------------
% Preallocate structure array and mat-file or prepare update
step0 = mod_init(update);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main routines ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Get parameters
load('param_eddy_tracking','path_out','streamlines','resol','stepF')

%----------------------------------------
% Build I/O matfile
disp('Your MATLAB to support "-v7.3" format to get full performance of')
disp('I/O MAT-file and save memory space')
disp(' ')

%----------------------------------------
% !!! Loop on steps using FOR or PARFOR (in hard) !!!
%----------------------------------------

diary([path_out,'log.txt'])

parfor stp = step0:stepF
    
    %----------------------------------------
    % Build I/O matfile
    fields_inter_mat = matfile([path_out,'fields_inter.mat'],'Writable',true);
    fields_mat = matfile([path_out,'fields.mat'],'Writable',true);
    centers_mat = matfile([path_out,'eddy_centers.mat'],'Writable',true);
    shapes_mat = matfile([path_out,'eddy_shapes.mat'],'Writable',true);

%for stp = step0:stepF
    
    %----------------------------------------
    % begin the log file
    %system(['mkdir ',path_out,'/log']);
    %if stp<10
    %    diary([path_out,'log/log_eddy_stp_000',num2str(stp),'.txt']);
    %elseif stp<100
    %    diary([path_out,'log/log_eddy_stp_00',num2str(stp),'.txt']);
    %elseif stp<1000
    %    diary([path_out,'log/log_eddy_stp_0',num2str(stp),'.txt']);
    %lseif stp<10000
    %    diary([path_out,'log/log_eddy_stp_',num2str(stp),'.txt']);
    %end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute LNAM ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------------------------
    % Compute non interpolated fields for step stp
    fields = mod_fields(source,stp,1);
    
    %----------------------------------------
    % Write in I/O matfile step stp
    fields_mat.detection_fields(:,stp) = fields;
    
    if resol>1
        %----------------------------------------
        % Compute interpolated fields for step stp
        fields = mod_fields(source,stp,resol);
    end
    
    %----------------------------------------
    % Write in I/O matfile
    % Interpolated and non interpolated field can be the same
    fields_inter_mat.detection_fields(:,stp) = fields;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find centers ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------------------------
    % Detection of eddy centers for step stp
    [centers0,centers] = mod_eddy_centers(source,stp,fields);
    
    %----------------------------------------
    % Write in I/O matfile step stp
    centers_mat.centers0(:,stp) = centers0;
    centers_mat.centers(:,stp) = centers;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------------------------
    % Determination of eddy features for step stp
    [centers2,shapes1,shapes2,profil2,warn_shapes,warn_shapes2] = ...
        mod_eddy_shapes(source,stp,fields,centers);
    
    %----------------------------------------
    % Write in I/O matfile step stp
    centers_mat.centers2(:,stp) = centers2;
    shapes_mat.shapes1(:,stp) = shapes1;
    shapes_mat.shapes2(:,stp) = shapes2;
    shapes_mat.warn_shapes(:,stp) = warn_shapes;
    shapes_mat.warn_shapes2(:,stp) = warn_shapes2;
    if streamlines
        shapes_mat.profil2(:,stp) = profil2;
    end
    
%----------------------------------------
% close log file
%diary off

end

diary off

%system(['cat ',path_out,'log/log_eddy_stp*.txt >> ',path_out,'log_eddy.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Track eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Tracking and record interacting events
mod_eddy_tracks(update)








