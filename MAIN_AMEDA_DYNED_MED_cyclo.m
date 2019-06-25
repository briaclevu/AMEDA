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

start
clear; clc;

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% domaine
keys = 'DYNED_MED_cyclo';

%----------------------------------------
% Update option
update = 0; % the serie from the beginning

%----------------------------------------
% Possibility to shorter the serie
%stepF = 10;

%----------------------------------------
% inital and final year
Yi = 2000; Yf = 2017;

%----------------------------------------
% Set parallel computation
cpus=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Produce default parameters in param_eddy_tracking
if exist('stepF','var')
    %mod_eddy_params(['keys_sources_',source,'_',keys],stepF)
else
    %mod_eddy_params(['keys_sources_',source,'_',keys])
end
run(['keys_sources_',source,'_',keys])
load('param_eddy_tracking','path_out','streamlines','resol','stepF');

%----------------------------------------
% list of steps (from 2000 to 2017)
% 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017
list=[1 367 732 1097 1462 1828 2193 2558 2923 3289 3654 4019 4384 4750 5115 5480 5845 6211 stepF+1];

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

%----------------------------------------
% process detection in yearly loops

for i=1:length(list)-1%(Yf-2000)+1%1:length(list)-1
    
    stepFF=list(i+1)-list(i);
    dstp = list(i)-1;

    %----------------------------------------
    % Preallocate structure array and mat-file or prepare update
    % !! replace or reinitialise previous results !!
    step0 = mod_init(stepFF,update);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute LNAM ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if i>18
     
    disp([' === Compute LNAM ',num2str(2000-1+i),' ==='])
    disp(' ')

    load([path_out,'fields'],'detection_fields')
    detection_fields_ni = detection_fields;

    load([path_out,'fields_inter.mat'],'detection_fields')

    parfor stp = step0:stepFF
        %----------------------------------------
        % Compute non interpolated fields for step stp
        detection_fields_ni(stp) = mod_fields(source,stp+dstp,1);
        if resol>1
            %----------------------------------------
            % Compute interpolated fields for step stp
            detection_fields(stp) = mod_fields(source,stp+dstp,resol);
        else
            %----------------------------------------
            % Interpolated and non interpolated field are the same
            disp(' === Interpolated LNAM is the same ===')
            detection_fields(stp) = detection_fields_ni(stp);
        end
    end

    %----------------------------------------
    % Save fields
    save([path_out,'fields_inter_',num2str(2000-1+i)],'detection_fields','-v7.3')

    detection_fields = detection_fields_ni;
    save([path_out,'fields_',num2str(2000-1+i)],'detection_fields','-v7.3')
    clear detection_fields detection_fields_ni

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find centers ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp([' === Find potential centers ',num2str(2000-1+i),' ==='])
    disp(' ')

    load([path_out,'eddy_centers'])

    %----------------------------------------
    % Build I/O matfile
    fields_mat = matfile([path_out,'fields_inter_',num2str(2000-1+i),'.mat']);

    parfor stp = step0:stepFF
        % load inter fields at step stp
        %----------------------------------------
        fields = fields_mat.detection_fields(:,stp);
        %----------------------------------------
        % Detection of eddy centers for step stp
        [centers0(stp),centers(stp)] = mod_eddy_centers(source,stp+dstp,fields);
    end

    %----------------------------------------
    % Save centers
    save([path_out,'eddy_centers_',num2str(2000-1+i)],'centers0','centers','-v7.3')
    clear centers0 centers
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp([' === Determine eddies shapes ',num2str(2000-1+i),' ==='])
    disp(' ')

    load([path_out,'eddy_centers'],'centers2')
    load([path_out,'eddy_shapes'])

    %----------------------------------------
    % Build I/O matfile
    fields_mat = matfile([path_out,'fields_inter_',num2str(2000-1+i),'.mat']);
    centers_mat = matfile([path_out,'eddy_centers_',num2str(2000-1+i),'.mat']);

    parfor stp = step0:stepFF
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
                mod_eddy_shapes(source,stp+dstp,fields,centers);
        else
            [centers2(stp),shapes1(stp),shapes2(stp),~,...
                warn_shapes(stp),warn_shapes2(stp)] = ...
                mod_eddy_shapes(source,stp+dstp,fields,centers);
        end
    end

    %----------------------------------------
    % save warnings, shapes and their centers
    save([path_out,'eddy_centers_',num2str(2000-1+i)],'centers2','-append')
    if streamlines
        save([path_out,'eddy_shapes_',num2str(2000-1+i)],'shapes1','shapes2',...
            'warn_shapes','warn_shapes2','profil2','-v7.3')
    else
        save([path_out,'eddy_shapes_',num2str(2000-1+i)],'shapes1','shapes2',...
            'warn_shapes','warn_shapes2','-v7.3')
    end
    clear centers2 shapes1 shapes2 profil2 warn_shapes warn_shapes2 struct1 struct2 struct3
end

end % end loop for year

%----------------------------------------
% Free workers
delete(gcp('nocreate'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Track eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name='_2012';
% 
% %----------------------------------------
% % Tracking eddies and record interacting events
% mod_eddy_tracks_nopool(name,update);
% 
% %----------------------------------------
% % Resolve merging and spltting event and filter eddies shorter than cut_off
% mod_merging_splitting(name);
% 


%----------------------------------------
% concatenete years
concat_eddy(num2cell(Yi:Yf))

%----------------------------------------
% Tracking eddies and record interacting events
name=['_',num2str(Yi),'_',num2str(Yf)];
% mod_eddy_tracks_nopool(name,update);
for i=1:length(list)-1
    update = list(end)-list(i);
    mod_eddy_tracks_nopool(name,update);
end

%----------------------------------------
% Resolve merging and spltting event and filter eddies shorter than cut_off
mod_merging_splitting(name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Reliability Index ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default set of parameter to compute density from satellite tracks
dT=12; Tk=2*dT;
D=4; Dk=3*D;

%----------------------------------------
%% build density tracks
dens=build_density(source,keys,Dk,D,Tk,dT);
save([path_out,'tracks_density_',config,'_T',num2str(dT),'_D',num2str(D)],'dens','-append','-v7.3')

%----------------------------------------
%% associate normalized density to tracks

% source of data driving the netcdf format
name=['_',num2str(Yi),'_',num2str(Yf)];
load([path_out,'eddy_tracks',name])  
tracksT=tracks;
% path_out_IE='/home/blevu/DATA/IeraPetra/datafig/';
% load([path_out_IE,'IEdt'])
% tracksIEDT(1)=ied;
% load([path_out_IE,'IEnrt'])
% tracksIENRT(1)=ied;

% tracks density adapted from Arhan and Colin de Verdière (1985)
load([path_out,'tracks_density_',config,'_T',num2str(dT),'_D',num2str(D)],'dens')
% i=datenum([2016 01 01])-datenum([2000 01 01])+1;
% dens=dens(:,:,i:end);
% save([path_out_IE,'tracks_density_',config,'_T',num2str(dT),'_D',num2str(D)'],'dens');

% compute density for every eddy detection (tracksT structure)
%dens_tracks_IEDT = associate_tracks_density(tracksIEDT,dens,D);
dens_tracks = associate_tracks_density(source,keys,tracksT,dens,D);
% save normalised density associate to detections
% save([path_out_IE,'tracks_density_',config,'_T',...
%         num2str(dT),'_D',num2str(D)],'dens_tracks_IEDT','-append')
save([path_out,'tracks_density_',config,'_T',...
        num2str(dT),'_D',num2str(D)],'dens_tracks','-append','-v7.3')

% % compute density for every eddy detection (tracksT structure)
% dens_tracks_IENRT = associate_tracks_density(tracksIENRT,dens,D);
% % save normalised density associate to detections
% save([path_out_IE,'tracks_density_',config,'_T',...
%         num2str(dT),'_D',num2str(D)],'dens_tracks_IENRT','-append')

%-------------------------------------------------------------------------
%% Compute RI stats on tracks

% probalistic function
RT_seuil=24;% 30(A) / 25-27(C) from 16 to 34
VT_seuil=0.08;% 0.12 independant pour partGhost(A) / 0.08-0.1(C)
dRT=4;% 5(A) / 3(C) from 0 to 10
dVT=0.04;% 0.04(a) / 0.03(C)

% set paramters
dT=12;
D=4;

% RI for DT
load([path_out_IE,'tracks_density_',config,'_T',...
        num2str(dT),'_D',num2str(D)],'dens_tracks_IEDT')
    
% read features
rmax1=dens_tracks_IEDT.rmax1;
velmax1=dens_tracks_IEDT.velmax1;
DensN=dens_tracks_IEDT.dens_int';

PHIR = 1/2*(tanh((rmax1-RT_seuil)/dRT)+1);
PHIV = 1/2*(tanh((velmax1-VT_seuil)/dVT)+1);

% Reliability index
RI_IEDT = DensN.*PHIR.*PHIV;

save([path_out_IE,'tracks_density_',config,'_T',...
        num2str(dT),'_D',num2str(D)],'RI_IEDT','-append')

% RI for NRT
load([path_out_IE,'tracks_density_',config,'_T',...
        num2str(dT),'_D',num2str(D)],'dens_tracks_IENRT')
    
% read features
rmax1=dens_tracks_IENRT.rmax1;
velmax1=dens_tracks_IENRT.velmax1;
DensN=dens_tracks_IENRT.dens_int;

PHIR = 1/2*(tanh((rmax1-RT_seuil)/dRT)+1);
PHIV = 1/2*(tanh((velmax1-VT_seuil)/dVT)+1);

% Reliability index
RI_IENRT = DensN.*PHIR.*PHIV;

save([path_out_IE,'tracks_density_',config,'_T',...
        num2str(dT),'_D',num2str(D)],'RI_IENRT','-append')

