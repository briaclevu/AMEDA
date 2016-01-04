function mod_fields(b,update)
%  mod_fields(b,update) creates 2D-t fields of various computation.
%  Use velocities from load_field to include a proper mask
%  with b pixel inside land and nan other land.
%  update is a flag allowing to update an existing tracking
%  (update is also = number of steps backward to consider)
%  Saved fields are divergence, velocity, vorticity, Okubo-Weiss and LNAM
%  (calculated using b parameter which could include the resolution)
%
%  For a description of the input parameters see param_eddy_tracking.m.
%
%  Apr.2015 Briac Le Vu
%
%=========================

% define constants
earth_radius = 6378137; % m
% meter (m) per degree of latitude
R=earth_radius*pi/180; % 111320m instead of 111100m

% preallocate fields struct array if doesn't exist
global path_out
global runname

load([path_out,'fields_',runname])
step = size(u,3);

if update
    step0=step-update+1;
else
    detection_fields = struct('ke',{},'div',{},'vort',{},'OW',{},'LOW',{},'LNAM',{});
    step0=1;
end

%----------------------------------------

disp('Compute fields on interpolated grid')

% cycle through time steps
for j=step0:step
    disp(['  Step ',num2str(j),' ... '])
    
    % load 2D velocity fields for stepj (m/s)
    %----------------------------------------
    uu=squeeze(u(:,:,j));
    vv=squeeze(v(:,:,j));

    % Calculation of eddy kinetic energy
    %----------------------------------------
    ke=(uu.^2+vv.^2)/2;
    
    % Calculation of finite spatial element
    %----------------------------------------
    dx=zeros(size(lon));
    dy=dx;
    dux=dx;
    duy=dx;
    dvx=dx;
    dvy=dx;
    
    % Spatial element in term of lon and lat in degree
    dx(2:end-1,2:end-1)=lon(2:end-1,3:end)-lon(2:end-1,1:end-2); %#ok<*COLND>
    dy(2:end-1,2:end-1)=lat(3:end,2:end-1)-lat(1:end-2,2:end-1);
    
    % Calcul finite spatial element meter (m) per degree
    dy=dy*R;
    dx=dx*R.*cosd(lat);
    
    % Compute speed element in m/s
    dux(2:end-1,2:end-1)=(uu(2:end-1,3:end)-uu(2:end-1,1:end-2));
    duy(2:end-1,2:end-1)=(uu(3:end,2:end-1)-uu(1:end-2,2:end-1));
    dvx(2:end-1,2:end-1)=(vv(2:end-1,3:end)-vv(2:end-1,1:end-2));
    dvy(2:end-1,2:end-1)=(vv(3:end,2:end-1)-vv(1:end-2,2:end-1));

    % Calculation of Okubo-Weiss criteria
    %----------------------------------------
    sn=(dux./dx)-(dvy./dy); % shear "cisaillement"
    ss=(dvx./dx)+(duy./dy); % strain "deformation"
    om=(dvx./dx)-(duy./dy); % vorticity "vorticité"

    okubo=sn.^2+ss.^2-om.^2; % in s-2

    % Calculation of divergence
    %----------------------------------------
    div=(dux./dx)+(dvy./dy);

    % Calculation of vorticity field (typically +/-10^-5 s-1)
    %----------------------------------------
    vorticity=om;

    % border is a parameter which prevents the constraints
    % to be applied to points too close to the domain boundaries
    % which would result in an index error
    borders=b+1;
    
    % Calculation of LNAM criteria (Local Normalized Angular Momentum)
    % and LOW criteria (Local Averaged Okubo Weiss)
    %----------------------------------------
    
    % Initialisation
    L=zeros(size(uu));
    LOW=nan(size(uu));

    % calculate LNAM and LOW in all domain pixels
    for i=borders:length(vv(:,1))-borders+1
        for ii=borders:length(vv(1,:))-borders+1

            if ~isnan(vv(i,ii))
                
                % calculate LOW
                OW=okubo(i-b:i+b,ii-b:ii+b);
                LOW(i,ii)=mean(OW(:));
                
                % calculate LNAM
                lonlocal=lon(i-b:i+b,ii-b:ii+b); % lon square sample of length b
                latlocal=lat(i-b:i+b,ii-b:ii+b); % lat square sample of length b
                
                ulocal=uu(i-b:i+b,ii-b:ii+b); % u square sample of length b
                vlocal=vv(i-b:i+b,ii-b:ii+b); % v square sample of length b

                % Local Normalized Angular Momentum
                % Use the midddle of the square as center coordinate
                coordcentre=size(ulocal,1)-b;
                d_loncentre=(lonlocal-lonlocal(coordcentre,coordcentre))*R...
                    .*cosd(latlocal(:,:));
                d_latcentre=(latlocal-latlocal(coordcentre,coordcentre))*R;

                cross=(d_loncentre.*vlocal)-(d_latcentre.*ulocal);
                dot=(ulocal.*d_loncentre)+(vlocal.*d_latcentre);
                produit=sqrt(ulocal.^2+vlocal.^2)...
                    .*sqrt(d_loncentre.^2+d_latcentre.^2);

                if (sum(dot(:))+sum(produit(:))) ~= 0
                    L(i,ii)=sum(cross(:))/(sum(dot(:))+sum(produit(:)));
                else
                    L(i,ii)=0;
                end
            end
        end
    end
    L(isnan(L))=0;

    % save fields in struct array
    %----------------------------------------
    detection_fields(j).ke = ke.*mask;
    detection_fields(j).div = div.*mask;
    detection_fields(j).vort = vorticity.*mask;
    detection_fields(j).OW = okubo.*mask;
    detection_fields(j).LOW = LOW.*mask;
    detection_fields(j).LNAM = L.*mask; % relative value of the LNAM
end

disp(' ')

% export fields
%----------------------------------------
save([path_out,'fields_',runname],'detection_fields','-append')

