% create folder of netcdf files from traks structure
start
clear all

% configuration to be used has reference of the v1 of the atlas
source = 'AVISO';
%config = 'DYNED_MED_adt';
%config = 'DYNED_MED_cyclo';
keys = 'DYNED_ARA_adt';

% Definition of the parameters specific for the experiment is needed 
run(['keys_sources_',source,'_',keys])
load('param_eddy_tracking')

% load result from AMEDA tracking and merging/splitting
if strcmp(config,'DYNED_MED_cyclo') || strcmp(config,'DYNED_MED')
    
    load([path_out,'eddy_tracks2_2000_2017'])
    %load([path_out,'eddy_tracks2_2010_2011'])
    tracks=tracks2;
    if streamlines
        load([path_out,'eddy_shapes_2000_2017'],'profil2')
        %load([path_out,'eddy_shapes_2010_2011'],'profil2')
    end

    % date and length from "YYYY-MM-DDThh:mm:ssZ"
    datei=datenum([1999 12 31]);
    duralim=0;

    % create folder for the experiment
    path_out2=[path_out,'nc_tracks_2000_2017/'];
    system(['mkdir ',path_out2]);
    
elseif strcmp(config,'DYNED_ARA') || strcmp(config,'DYNED_ARA_cyclo')
    
    load([path_out,'eddy_tracks2_2000_2015'])
    tracks=tracks2;
    if streamlines
        load([path_out,'eddy_shapes_2000_2015'],'profil2')
    end

    % date and length from "YYYY-MM-DDThh:mm:ssZ"
    datei=datenum([1999 12 31]);
    duralim=0;

    % create folder for the experiment
    path_out2=[path_out,'nc_tracks_2000_2015/'];
    system(['mkdir ',path_out2]);
    
end

% load test to get _Fillvalue
FV = double(ncread(['/home/blevu/Resultats/AVISO/MED/eddy_track_test.nc'],'x_max'));
FV1 = -2147483648;
FV2 = -1*FV(end);
NV1 = [];
NV2 = [];
NV3 = [];

% scan and record eddy longer than turnover time (tindracks2 and duralim==0)
for i=1:length(tracks)
    
    %display(['track ',num2str(i)])
    
    %if length(tracks(i).step)>nanmean(tracks(i).tau1)
    if tracks(i).step(end) - tracks(i).step(1)>=duralim
        
        % ID tracks
        if strcmp(config,'DYNED_MED_cyclo') || strcmp(config,'DYNED_MED')
            if i<10
                ID=['0000',num2str(i)];
            elseif i<100
                ID=['000',num2str(i)];
            elseif i<1000
                ID=['00',num2str(i)];
            elseif i<10000
                ID=['0',num2str(i)];
            else
                ID=num2str(i);
            end
        elseif strcmp(config,'DYNED_ARA') || strcmp(config,'DYNED_ARA_cyclo')
            if i<10
                ID=['2000',num2str(i)];
            elseif i<100
                ID=['200',num2str(i)];
            elseif i<1000
                ID=['20',num2str(i)];
            elseif i<10000
                ID=['2',num2str(i)];
            else
                ID=num2str(i+20000);
            end
        end        
        % steps of tracking
        time=tracks(i).step;
                
        % count total number of vertices for max velocity
        clear NV1 NV2 NV3
        for j=1:length(time)
            if isnan(tracks(i).rmax1(j)) || isnan(tracks(i).shapes1{j}(1))
                NV1(j) = 0;
                disp(['no max ',num2str([i,j])])
            else
                NV1(j) = length(tracks(i).shapes1{j});
            end
            if isnan(tracks(i).rmax3(j)) || isnan(tracks(i).shapes3{j}(1)) 
                NV2(j) = 0;
                disp(['no end ',num2str([i,j])])
            else
                if any(tracks(i).rmax3./tracks(i).rmax1<1)
                    NV2(j) = length(tracks(i).shapes1{j});
                else
                    NV2(j) = length(tracks(i).shapes3{j});
                end
            end
            if streamlines
                ind = tracks(i).ind(j);
                NV3(j) = length(profil2(time(j)).rmoy{ind});
                if NV3(j)<1
                    disp(['PROFIL default value tracks ',num2str(i),' step ',num2str(ind(j))])
                end
            end
        end
        
        % file with eddy ID in the name
        nc_file=[path_out2,'eddy_track_',ID,'.nc'];

        % Get the numeric values corresponding to the NETCDF4 and CLASSIC_MODEL constants defined by the NetCDF library.
        %
        cmode = netcdf.getConstant('NETCDF4');
        cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

        % Create a NetCDF-4 file that uses the classic model by specifying the creation mode value, cmode.
        %
        nc = netcdf.create(nc_file,cmode);

        %
        %  Create dimensions
        %
        dimid(1) = netcdf.defDim(nc,'time',length(time));
        % contour max vertices
        dimid(2) = netcdf.defDim(nc,'nv1',sum(NV1));% as CORIOLIS time-series recommendation 
        % contour end vertices
        dimid(3) = netcdf.defDim(nc,'nv2',sum(NV2));% as CORIOLIS time-series recommandation
        % profils V-R vertices
        if streamlines
            dimid(4) = netcdf.defDim(nc,'nv3',sum(NV3));% as CORIOLIS time-series recommandation
        end
        
        %
        %  Create variables
        %
        varid(1) = netcdf.defVar(nc,'time_step','float',dimid(1));
        varid(2) = netcdf.defVar(nc,'x_cen','double',dimid(1));
        varid(3) = netcdf.defVar(nc,'y_cen','double',dimid(1));
        varid(4) = netcdf.defVar(nc,'x_bar','double',dimid(1));
        varid(5) = netcdf.defVar(nc,'y_bar','double',dimid(1));
        varid(6) = netcdf.defVar(nc,'tau_min','double',dimid(1));
        %
        varid(10) = netcdf.defVar(nc,'n_max','int',dimid(1));
        varid(11) = netcdf.defVar(nc,'x_max','double',dimid(2));
        varid(12) = netcdf.defVar(nc,'y_max','double',dimid(2));
        varid(13) = netcdf.defVar(nc,'v_max','double',dimid(1));
        varid(14) = netcdf.defVar(nc,'r_max','double',dimid(1));
        varid(15) = netcdf.defVar(nc,'dssh_max','double',dimid(1));
        varid(16) = netcdf.defVar(nc,'aire_max','double',dimid(1));
        varid(17) = netcdf.defVar(nc,'ellip_max','double',dimid(1));
        varid(18) = netcdf.defVar(nc,'theta_max','double',dimid(1));
        varid(19) = netcdf.defVar(nc,'contour_class','int',dimid(1));
        varid(20) = netcdf.defVar(nc,'Ro','double',dimid(1));
        varid(21) = netcdf.defVar(nc,'VortM','double',dimid(1));
        %
        varid(30) = netcdf.defVar(nc,'n_end','int',dimid(1));
        varid(31) = netcdf.defVar(nc,'x_end','double',dimid(3));
        varid(32) = netcdf.defVar(nc,'y_end','double',dimid(3));
        varid(33) = netcdf.defVar(nc,'v_end','double',dimid(1));
        varid(34) = netcdf.defVar(nc,'r_end','double',dimid(1));
        varid(35) = netcdf.defVar(nc,'dssh_end','double',dimid(1));
        varid(36) = netcdf.defVar(nc,'aire_end','double',dimid(1));
        varid(37) = netcdf.defVar(nc,'alpha','double',dimid(1));
        %
        if streamlines
            varid(40) = netcdf.defVar(nc,'n_pro','int',dimid(1));
            varid(41) = netcdf.defVar(nc,'v_pro','double',dimid(4));
            varid(42) = netcdf.defVar(nc,'r_pro','double',dimid(4));
        end
        %
        varid(50) = netcdf.defVar(nc,'splitting_ID','int',dimid(1));
        varid(51) = netcdf.defVar(nc,'merging_ID','int',dimid(1));
%        varid(52) = netcdf.defVar(nc,'argo_ID','int',dimid(1));
%         varid(53) = netcdf.defVar(nc,'satellite_ID','int',dimid(1));
%         varid(54) = netcdf.defVar(nc,'interaction_ID','int',dimid(1));
        %
        %  Create attributes
        %
        netcdf.putAtt(nc,varid(1),'standard_name','time');
        netcdf.putAtt(nc,varid(1),'units',['days since ',datestr(datei,'yyyy-mm-ddTHH:MM:SSZ')]);
        netcdf.putAtt(nc,varid(1),'axis','T');
        netcdf.putAtt(nc,varid(1),'long_name','time in AVISO files');
        %
        netcdf.putAtt(nc,varid(2),'standard_name','longitude');
        netcdf.putAtt(nc,varid(2),'units','degrees_east');
        netcdf.putAtt(nc,varid(2),'axis','X');
        netcdf.putAtt(nc,varid(2),'long_name','longitude of the eddy center');
        %
        netcdf.putAtt(nc,varid(3),'standard_name','latitude');
        netcdf.putAtt(nc,varid(3),'units','degrees_north');
        netcdf.putAtt(nc,varid(3),'axis','Y');
        netcdf.putAtt(nc,varid(3),'long_name','latitude of the eddy center');
        %
        netcdf.putAtt(nc,varid(4),'standard_name','longitude');
        netcdf.putAtt(nc,varid(4),'units','degrees_east');
        netcdf.putAtt(nc,varid(4),'axis','X');
        netcdf.putAtt(nc,varid(4),'long_name','longitude of the barycenter');
        %
        netcdf.putAtt(nc,varid(5),'standard_name','latitude');
        netcdf.putAtt(nc,varid(5),'units','degrees_north');
        netcdf.putAtt(nc,varid(5),'axis','Y');
        netcdf.putAtt(nc,varid(5),'long_name','latitude of the barycenter');
        %
        netcdf.putAtt(nc,varid(6),'long_name','minimal eddy turnover time ( min(perimeter/v) )');
        netcdf.putAtt(nc,varid(6),'units','days');
        netcdf.putAtt(nc,varid(6),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(10),'long_name','number of vertices for the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(10),'sample_dimension','nv1');
        %
        netcdf.putAtt(nc,varid(11),'long_name','vertices longitude for the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(11),'units','degrees_east');
        %
        netcdf.putAtt(nc,varid(12),'long_name','vertices latitude for the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(12),'units','degrees_north');
        %
        netcdf.putAtt(nc,varid(13),'long_name','maximal mean velocity along a close streamline');
        netcdf.putAtt(nc,varid(13),'units','m/s');
        netcdf.putAtt(nc,varid(13),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(14),'long_name','equivalent radius for the streamline with the maximum mean velocity ( sqrt(aire_max/pi) )');
        netcdf.putAtt(nc,varid(14),'units','km');
        netcdf.putAtt(nc,varid(14),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(15),'long_name','delta ssh inside the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(15),'units','m');
        netcdf.putAtt(nc,varid(15),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(16),'long_name','surface inside the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(16),'units','km^2');
        netcdf.putAtt(nc,varid(16),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(17),'long_name','ellipticity (1-b/a) fitted on the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(17),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(18),'long_name','angle between the longer axe of the ellipse (a) with the longitude westward axe');
        netcdf.putAtt(nc,varid(18),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(19),'long_name','is v_max a true maximum value with a 3% decrease between the streamline with the maximum mean velocity and the last streamline');
        netcdf.putAtt(nc,varid(19),'convention','1 if true 0 if not');
        netcdf.putAtt(nc,varid(19),'_Fillvalue',FV1);
        %
        netcdf.putAtt(nc,varid(20),'long_name','Rossby number deduce from r_max and v_max');
        netcdf.putAtt(nc,varid(20),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(21),'long_name','maximal vorticity inside the streamline with the maximum mean velocity');
        netcdf.putAtt(nc,varid(21),'units','s-1');
        netcdf.putAtt(nc,varid(21),'_Fillvalue',FV2);
        %
        %
        netcdf.putAtt(nc,varid(30),'long_name','number of vertices for the last streamline with one center');
        netcdf.putAtt(nc,varid(30),'sample_dimension','nv2');
        %
        netcdf.putAtt(nc,varid(31),'long_name','vertices longitude for the last streamline with one center');
        netcdf.putAtt(nc,varid(31),'units','degrees_east');
        %
        netcdf.putAtt(nc,varid(32),'long_name','vertices latitude for the last streamline with one center');
        netcdf.putAtt(nc,varid(32),'units','degrees_north');
        %
        netcdf.putAtt(nc,varid(33),'long_name','mean velocity for the last streamline with one center');
        netcdf.putAtt(nc,varid(33),'units','m/s');
        netcdf.putAtt(nc,varid(33),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(34),'long_name','equivalent radius for the last streamline with one center ( sqrt(aire_end/pi) )');
        netcdf.putAtt(nc,varid(34),'units','km');
        netcdf.putAtt(nc,varid(34),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(35),'long_name','delta ssh inside the last streamline with one center');
        netcdf.putAtt(nc,varid(35),'units','m');
        netcdf.putAtt(nc,varid(35),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(36),'long_name','surface of the last streamline with one center');
        netcdf.putAtt(nc,varid(36),'units','km^2');
        netcdf.putAtt(nc,varid(36),'_Fillvalue',FV2);
        %
        netcdf.putAtt(nc,varid(37),'long_name','shape parameter of the v-r profile at the surface');
        netcdf.putAtt(nc,varid(37),'convention','v/v_max = r/r_max.exp( (1-(r/r_max)^alpha) / alpha )');
        netcdf.putAtt(nc,varid(37),'_Fillvalue',FV2);
        %
        %
        if streamlines
            netcdf.putAtt(nc,varid(40),'long_name','number of values for the mean velocity (V) versus mean radius (R) profil');
            netcdf.putAtt(nc,varid(40),'sample_dimension','nv3');
            %
            netcdf.putAtt(nc,varid(41),'long_name','mean velocity (V) along the V-R profil');
            netcdf.putAtt(nc,varid(41),'units','m/s');
            %
            netcdf.putAtt(nc,varid(42),'long_name','mean radius (R) along the V-R profil');
            netcdf.putAtt(nc,varid(42),'units','km');
        end
        %
        %
        netcdf.putAtt(nc,varid(50),'long_name','interacting eddy ID which split');
        netcdf.putAtt(nc,varid(50),'_Fillvalue',FV1);
        %
        netcdf.putAtt(nc,varid(51),'long_name','interacting eddy ID which merge');
        netcdf.putAtt(nc,varid(51),'_Fillvalue',FV1);
        %
%         netcdf.putAtt(nc,varid(52),'long_name','argo profile ID inside the last streamline');
%         netcdf.putAtt(nc,varid(52),'_Fillvalue',FV1);
        %
%         netcdf.putAtt(nc,varid(53),'long_name','satellite tracks ID crossing the last streamline');
%         netcdf.putAtt(nc,varid(53),'_Fillvalue',FV1);
%         %
%         netcdf.putAtt(nc,varid(54),'long_name','Eddy ID interacting');
%         netcdf.putAtt(nc,varid(54),'_Fillvalue',FV1);
        %
        %  Create global attributes
        %
        if strcmp(config,'DYNED_MED_cyclo') || strcmp(config,'DYNED_MED')
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',[date,' make_netcdf_from_tracks.m routine on tracks2_2000_2017 with eddies longer than 2*taumin']);
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',['eddy tracks ID ',ID,' from AMEDA']);
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'summary','This dataset contains AMEDA analysis of 1/8° SSALTO/DUACS products over Mediterranean Sea (http://www.aviso.altimetry.fr)');
        elseif strcmp(config,'DYNED_ARA') || strcmp(config,'DYNED_ARA_cyclo')
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',[date,' make_netcdf_from_tracks.m routine on tracks2_2000_2015 with eddies longer than 2*taumin']);
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',['eddy tracks ID ',ID,' from AMEDA']);
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'summary','This dataset contains AMEDA analysis of 1/8^o SSALTO/DUACS products over Arabian Sea (http://www.aviso.altimetry.fr)');
        end
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'project','DYNED-Atlas');
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','LMD');
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'contact','briac.le-vu@lmd.polytechnique.fr');
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'date_created',date);
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'track_AMEDA_ID',ID);
        if tracks(i).type(1) == 1
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_type','cyclone');
        else
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_type','anticyclone');
        end
        % set less than 20% of the time without clear contour as gure
%         if mean(tracks(i).large1) >= 0.8
%             netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_class','gyre (less than 20% contour_class=1)');
%         else
%             netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_class','eddy (more than 20% contour_class=1)');
%         end
        % eddy origin to be build
        serie=find(tracks(i).split==1);
        if ~isempty(serie)
            ind = unique(tracks(i).interaction(serie));
            for j=1:length(ind)
                serie2 = find(tracks(ind(j)).interaction==i);
                if any(tracks(ind(j)).split(serie2)==1)
                    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_origin','0')
                else
                    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_origin',num2str(ind(j)))
                    break
                end
            end
        else
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_origin','0')
        end
        % eddy death to be build
        serie=find(tracks(i).merge==1);
        if ~isempty(serie)
            ind = unique(tracks(i).interaction(serie));
            for j=1:length(ind)
                serie2 = find(tracks(ind(j)).interaction==i);
                if any(tracks(ind(j)).split(serie2)==1)
                    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_death','0')
                else
                    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_death',num2str(ind(j)))
                    break
                end
            end
        else
            netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'eddy_death','0')
        end
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'number_of_detections',num2str(length(time)));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'life_time',num2str(time(end)-time(1)+1));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'date_first_detection',datestr(time(1)+datei,'yyyy-mm-dd'));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'date_last_detection',datestr(time(end)+datei,'yyyy-mm-dd'));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lon_min',num2str(min(tracks(i).x1)));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lon_max',num2str(max(tracks(i).x1)));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lat_min',num2str(min(tracks(i).y1)));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lat_max',num2str(max(tracks(i).y1)));
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lon_first_detection',num2str(tracks(i).x1(1)));% for v2
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lon_last_detection',num2str(tracks(i).x1(end)));% for v2
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lat_first_detection',num2str(tracks(i).y1(1)));% for v2
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'center_lat_last_detection',num2str(tracks(i).y1(end)));% for v2
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'r_max_avg_in_km',num2str(nanmean(tracks(i).rmax1)));% for v2
        netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'v_max_avg_in_m_per_s',num2str(nanmean(tracks(i).velmax1)));% for v2
        % list of splited and merged eddies with the track !!!make list of non interacting eddies - à revoir...!!!
        UE = unique(tracks(i).interaction(~isnan(tracks(i).interaction)));
        SE = '';
        ME = '';
        for j=1:length(UE)
            indUSE = find(tracks(UE(j)).split==1);
            indUME = find(tracks(UE(j)).merge==1);
            if ~isempty(indUSE)
                indSE = tracks(i).step==tracks(UE(j)).step(indUSE(1));
                if tracks(i).interaction(indSE)==UE(j)
                    if strcmp(SE,'')
                        SE = num2str(UE(j));
                    else
                        SE = [num2str(SE),',',num2str(UE(j))];
                    end
                end
            end
            if ~isempty(indUME)
                indME = tracks(i).step==tracks(UE(j)).step(indUME(end));
                if tracks(i).interaction(indME)==UE(j)
                    if strcmp(ME,'')
                        ME = num2str(UE(j));
                    else
                        ME = [num2str(ME),',',num2str(UE(j))];
                    end
                end
            end
        end
%         netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'list_of_eddies_which_split',SE);
%         netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'list_of_eddies_which_merge',ME);
        %
        % Leave define mode
        %
        netcdf.endDef(nc);
        %
        % Write variables
        %
        netcdf.putVar(nc,varid(1),time);
        netcdf.putVar(nc,varid(2),round(tracks(i).x1*1000)/1000);
        netcdf.putVar(nc,varid(3),round(tracks(i).y1*1000)/1000);
        netcdf.putVar(nc,varid(4),round(tracks(i).xbary1*1000)/1000);
        netcdf.putVar(nc,varid(5),round(tracks(i).ybary1*1000)/1000);
        %
        if any(isnan(tracks(i).tau1)) || any(tracks(i).tau1<=0)
            ind=find(isnan(tracks(i).tau1) | tracks(i).tau1<=0);
            for j=1:length(ind)
                disp(['TAU1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        netcdf.putVar(nc,varid(6),round(tracks(i).tau1*1000)/1000);
        %
        netcdf.putVar(nc,varid(10),NV1);
        n=0;
        for j=1:length(time)
            S=tracks(i).shapes1{j};
            if NV1(j)>0
                netcdf.putVar(nc,varid(11),n,NV1(j),round(S(1,:)*1000)/1000);
                netcdf.putVar(nc,varid(12),n,NV1(j),round(S(2,:)*1000)/1000);
            else
                disp(['SHAPES1 default tracks ',num2str(i),' step ',num2str(j)])
            end
            n=n+NV1(j);
        end
        %
        if any(isnan(tracks(i).velmax1)) || any(tracks(i).velmax1<=0)
            ind=find(isnan(tracks(i).velmax1) | tracks(i).velmax1<=0);
            for j=1:length(ind)
                disp(['VELMAX1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        netcdf.putVar(nc,varid(13),round(tracks(i).velmax1*10000)/10000);
        %
        if any(isnan(tracks(i).rmax1)) || any(tracks(i).rmax1<=0)
            ind=find(isnan(tracks(i).rmax1) | tracks(i).rmax1<=0);
            for j=1:length(ind)
                disp(['RMAX1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        netcdf.putVar(nc,varid(14),round(tracks(i).rmax1*100)/100);
        %
        DE=ones(length(time),1)*FV2;
        for j=1:length(time)
            if tracks(i).type(1) == 1 %cyclone
                if tracks(i).deta1(j)<0
                    DE(j)=tracks(i).deta1(j);
                else
                    display(['DETA1 default value for cyclone ',num2str(i),' step ',num2str(j)])
                end
            else % anticylone
                if tracks(i).deta1(j)>0
                    DE(j)=tracks(i).deta1(j);
                else
                    display(['DETA1 default value for anticyclone ',num2str(i),' step ',num2str(j)])
                end
            end
        end
        netcdf.putVar(nc,varid(15),round(DE*1000)/1000);
        %
        if any(isnan(tracks(i).aire1)) || any(tracks(i).aire1<=0)
            ind=find(isnan(tracks(i).aire1) | tracks(i).aire1<=0);
            for j=1:length(ind)
                disp(['AIRE1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        netcdf.putVar(nc,varid(16),round(tracks(i).aire1*10)/10);
        %
        if any(isnan(tracks(i).ellip1)) || any(tracks(i).ellip1<0)
            ind=find(isnan(tracks(i).ellip1) | tracks(i).ellip1<0);
            for j=1:length(ind)
                disp(['ELLIP1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        netcdf.putVar(nc,varid(17),round(tracks(i).ellip1*100)/100);
        %
        if any(isnan(tracks(i).theta1))
            ind=find(isnan(tracks(i).theta1));
            for j=1:length(ind)
                disp(['THETA1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
%         % !!! only to resolve one shot bug in mod_eddy_shapes
%         for j=1:length(time)
%             ll=tracks(i).shapes1{j};
%             if ~isempty(ll)
%                 [~,~,~,a,b,theta,~]=compute_ellip(ll);
%                 if a<b && b~=0
%                     if theta > pi/2
%                         tracks(i).theta1(j) = theta-pi/2;
%                     elseif theta < pi/2
%                         tracks(i).theta1(j) = theta+pi/2;
%                     end
%                 end
%             end
%         end
        if length(tracks(i).theta1)~=length(tracks(i).ellip1)
            tracks(i).theta1=tracks(i).theta1(1:length(tracks(i).ellip1));
        end
        netcdf.putVar(nc,varid(18),round(tracks(i).theta1*100)/100);
        %
        netcdf.putVar(nc,varid(19),uint8(abs(tracks(i).large1-1)));
        %
        % rossby radius
        f = 2*7.2921e-5*sind(tracks(i).ybary1);
        Ro1 = tracks(i).velmax1./f./tracks(i).rmax1*1e-3;
        if any(isnan(Ro1)) || any(Ro1<=0)
            ind=find(isnan(Ro1) | Ro1<=0);
            for j=1:length(ind)
                disp(['Ro1 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        netcdf.putVar(nc,varid(20),round(tracks(i).type.*Ro1*1000)/1000);
        %
        % voricity maximal
        VO=ones(length(time),1)*FV2;
        for j=1:length(time)
            if tracks(i).type(1) == 1 %cyclone
                if tracks(i).vortM1(j)>0
                    VO(j)=tracks(i).vortM1(j);
                else
                    display(['VORTM1 default value for cyclone ',num2str(i),' step ',num2str(j)])
                end
            else % anticylone
                if tracks(i).vortM1(j)<0
                    VO(j)=tracks(i).vortM1(j);
                else
                    display(['VORTM1 default value for anticyclone ',num2str(i),' step ',num2str(j)])
                end
            end
        end
        netcdf.putVar(nc,varid(21),round(VO*1e8)/1e8);
        %
        netcdf.putVar(nc,varid(30),NV2);
        n = 0;
        for j=1:length(time)
            if any(tracks(i).rmax3./tracks(i).rmax1<1)
                disp(['RMAX3 default value tracks ',num2str(i),' step ',num2str(ind(j))])
                %S=tracks(i).shapes1{j};
            else
                S=tracks(i).shapes3{j};
            end
            if NV2(j)>0
                netcdf.putVar(nc,varid(31),n,NV2(j),round(S(1,:)*1000)/1000);
                netcdf.putVar(nc,varid(32),n,NV2(j),round(S(2,:)*1000)/1000);
            else
                disp(['SHAPES3 default tracks ',num2str(i),' step ',num2str(j)])
            end
            n=n+NV2(j);
        end
        %
        if any(isnan(tracks(i).velmax3)) || any(tracks(i).velmax3<=0)
            ind=find(isnan(tracks(i).velmax3) | tracks(i).velmax3<=0);
            for j=1:length(ind)
                disp(['VELMAX3 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        if any(isnan(tracks(i).rmax3)) || any(tracks(i).rmax3<=0)
            %netcdf.putVar(nc,varid(33),round(tracks(i).velmax1*10000)/10000);
        else
            netcdf.putVar(nc,varid(33),round(tracks(i).velmax3*10000)/10000);
        end
        %
        if any(isnan(tracks(i).rmax3)) || any(tracks(i).rmax3<=0)
            ind=find(isnan(tracks(i).rmax3) | tracks(i).rmax3<=0);
            for j=1:length(ind)
                disp(['RMAX3 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        if any(tracks(i).rmax3./tracks(i).rmax1<1)
            ind=find(tracks(i).rmax3./tracks(i).rmax1<1);
            for j=1:length(ind)
                disp(['RMAX3 smaller than RMAX1 tracks ',num2str(i),' step ',num2str(ind(j))])
            end
            %netcdf.putVar(nc,varid(34),round(tracks(i).rmax1*100)/100);
        else
            netcdf.putVar(nc,varid(34),round(tracks(i).rmax3*100)/100);
        end
        %
        DE=ones(length(time),1)*FV2;
        for j=1:length(time)
            if tracks(i).type(1) == 1 %cyclone
                if tracks(i).deta3(j)<0
                    if any(tracks(i).rmax3./tracks(i).rmax1<1)
                        disp(['RMAX3 default value tracks ',num2str(i),' step ',num2str(ind(j))])
                        %DE(j)=tracks(i).deta1(j);
                    else
                        DE(j)=tracks(i).deta3(j);
                    end
                else
                    display(['DETA3 default value for cyclone ',num2str(i),' step ',num2str(j)])
                end
            else % anticylone
                if tracks(i).deta3(j)>0
                    if any(tracks(i).rmax3./tracks(i).rmax1<1)
                        disp(['RMAX3 default value tracks ',num2str(i),' step ',num2str(ind(j))])
                        %DE(j)=tracks(i).deta1(j);
                    else
                        DE(j)=tracks(i).deta3(j);
                    end
                else
                    display(['DETA3 default value for anticyclone ',num2str(i),' step ',num2str(j)])
                end
            end
        end
        netcdf.putVar(nc,varid(35),round(DE*10000)/10000);
        %
        if any(isnan(tracks(i).aire3)) || any(tracks(i).aire3<=0)
            ind=find(isnan(tracks(i).aire3) | tracks(i).aire3<=0);
            for j=1:length(ind)
                disp(['AIRE3 default value tracks ',num2str(i),' step ',num2str(ind(j))])
            end
        end
        if any(tracks(i).aire3./tracks(i).aire1<1)
            ind=find(tracks(i).aire3./tracks(i).aire1<1);
            for j=1:length(ind)
                disp(['AIRE3 smaller than AIRE1 tracks ',num2str(i),' step ',num2str(ind(j))])
            end
            %netcdf.putVar(nc,varid(36),round(tracks(i).aire1*10)/10);
        else
            netcdf.putVar(nc,varid(36),round(tracks(i).aire3*10)/10);
        end
        %
        if any(tracks(i).alpha<=0)
            ind=find(tracks(i).alpha<=0);
            for j=1:length(ind)
                disp(['ALPHA <0 tracks ',num2str(i),' step ',num2str(ind(j))])
                tracks(i).alpha(ind(j)) = FV2;
            end
        end
        ind = isnan(tracks(i).alpha);
        tracks(i).alpha(ind) = FV2;
        netcdf.putVar(nc,varid(37),round(tracks(i).alpha*100)/100);
        if length(tracks(i).step) > 200 && nanmean(tracks(i).rmax1) > 30 && length(find(~isnan(tracks(i).alpha))) > 10
            disp(['--> Lifetime ',num2str(length(tracks(i).step)),...
            ' | Meansize ',num2str(nanmean(tracks(i).rmax1)),...
            ' | Alpharecord ',num2str(length(find(~isnan(tracks(i).alpha))))])
        end
        %
        if streamlines
            netcdf.putVar(nc,varid(40),NV3);
            n=0;
            for j=1:length(time)
                ind = tracks(i).ind(j);
                V = profil2(time(j)).vel{ind};
                R = profil2(time(j)).rmoy{ind};
                if any(isnan(V)) || any(isnan(R))
                    disp(['V-R default value tracks ',num2str(i),' step ',num2str(ind(j))])
                end
                if NV3(j)>0
                    netcdf.putVar(nc,varid(41),n,NV3(j),round(V*10000)/10000);
                    netcdf.putVar(nc,varid(42),n,NV3(j),round(R*100)/100);
                end
                n=n+NV3(j);
            end
        end
        %
        UE = unique(tracks(i).interaction(~isnan(tracks(i).interaction)));
        SE = zeros(length(tracks(i).interaction),1);
        ME = SE;
        for j=1:length(UE)
            indUSE = find(tracks(UE(j)).split==1);
            indUME = find(tracks(UE(j)).merge==1);
            if ~isempty(indUSE)
                indSE = tracks(i).step==tracks(UE(j)).step(indUSE(1));
                if tracks(i).interaction(indSE)==UE(j)
                    SE(indSE) = UE(j);
                end
            end
            if ~isempty(indUME)
                indME = tracks(i).step==tracks(UE(j)).step(indUME(end));
                if tracks(i).interaction(indME)==UE(j)
                    ME(indME) = UE(j);
                end
            end
        end
        netcdf.putVar(nc,varid(50),SE);
        netcdf.putVar(nc,varid(51),ME);
%         netcdf.putAtt(nc,varid(52),);% argo_ID fill with Cori
%         netcdf.putAtt(nc,varid(53),'long_name','satellite tracks ID crossing the last streamline');
%         netcdf.putVar(nc,varid(54),tracks(i).interaction);
        %
        netcdf.close(nc)
    end
end

%% export only some tracks
% tracks2=tracks;
% 
% duralim=180;
% 
% for i=1:length(tracks)%:-1:length(tracks)-10
%     if tracks(i).step(end) - tracks(i).step(1)<duralim
%         %remove tracks(i)
%         tracks2(i) = [];
%         %short(i)=[];
%     end
% end
% tracks=tracks2;
% %save([path_out,'eddy_tracks2_2000_2015_6m'],'tracks','short','-v7.3')
% save([path_out,'eddy_tracks2_2010_2011_6m'],'tracks','short','-v7.3')
% 

