function concat_eddy(name)
%concat_eddy(name)
%
% horizonatlly concatene mod_fields and mod_eddy results identified by
% 'name' in a mutliyear computation 
%
%-------------------------
%   Feb 2016 Briac Le Vu
%-------------------------
%
%=========================

load('param_eddy_tracking','path_out','streamlines')

disp(['Concatenate year ',num2str(name{1}),'...'])

% % load first 
% load([path_out,'fields_',num2str(name{1})])
% detection_fields_1=detection_fields;

%load([path_out,'eddy_centers_',num2str(name{1})])
%centers2_1=centers2;
%centers_1=centers;
%centers0_1=centers0;

load([path_out,'eddy_shapes_',num2str(name{1})],'profil2')
%shapes1_1=shapes1;
%shapes2_1=shapes2;
%warn_1=warn_shapes;
%warn2_1=warn_shapes2;
if streamlines
    profil2_1=profil2;
end

for i=2:length(name)
    
    disp(['Concatenate year ',num2str(name{i}),'...'])
    
    % load fields, centers, shapes and warns
%     load([path_out,'fields_',num2str(name{i})])
%     detection_fields_2=detection_fields;
    
%    load([path_out,'eddy_centers_',num2str(name{i})])
%    centers2_2=centers2;
%    centers_2=centers;
%    centers0_2=centers0;
    
    load([path_out,'eddy_shapes_',num2str(name{i})],'profil2')
%    shapes1_2=shapes1;
%    shapes2_2=shapes2;
%    warn_2=warn_shapes;
%    warn2_2=warn_shapes2;
    if streamlines
        profil2_2=profil2;
    end

    % concate fields, centers, shapes and warns
%     detection_fields_1=[detection_fields_1,detection_fields_2];

%    centers0_1=[centers0_1,centers0_2];
%    centers_1=[centers_1,centers_2];
%    centers2_1=[centers2_1,centers2_2];

%    shapes1_1=[shapes1_1,shapes1_2];
%    shapes2_1=[shapes2_1,shapes2_2];

%    warn_1=[warn_1,warn_2];
%    warn2_1=[warn2_1,warn2_2];
    
    if streamlines
        profil2_1=[profil2_1,profil2_2];
    end

end

% straight step value
N=length(profil2);
for i=1:N
    %detection_fields_1(i).step=i;
%    centers0_1(i).step=i;
%    centers_1(i).step=i;
%    centers2_1(i).step=i;
%    shapes1_1(i).step=i;
%    shapes2_1(i).step=i;
    if streamlines
        profil2_1(i).step=i;
    end
end

% save concatenation
% detection_fields=detection_fields_1;
% save([path_out,'fields_',num2str(name{1}),'_',num2str(name{end})],'detection_fields','-v7.3')

%centers2=centers2_1;
%centers=centers_1;
%centers0=centers0_1;
%save([path_out,'eddy_centers_',num2str(name{1}),'_',num2str(name{end})],...
%    'centers0','centers','centers2','-v7.3')

%shapes1=shapes1_1;
%shapes2=shapes2_1;
%warn_shapes=warn_1;
%warn_shapes2=warn2_1;
%save([path_out,'eddy_shapes_',num2str(name{1}),'_',num2str(name{end})],...
%     'shapes1','shapes2','warn_shapes','warn_shapes2','-v7.3')

if streamlines
    profil2=profil2_1;
    save([path_out,'eddy_profil_',num2str(name{1}),'_',num2str(name{end})],...
        'profil2','-v7.3')
end


