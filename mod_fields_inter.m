function mod_fields(update)
%mod_fields({update})
%
%  Creates 2D-t fields of various computation by finite spatial element
%
%  Use velocities save in [path_out,'fields_',runname] by load_field.m
%  including a proper mask with the first b pixel of land (must include
%  the interpolation factor 'res' if any) and NaN in other land pixel.
%
%  update is a flag allowing to update an existing tracking:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the time serie
%
%  For a description of the input parameters see param_eddy_tracking.m.

%
%  Computed 2-Dfields are saved/update in [path_out,'fields_',runname] as:
%  - detection_fields(t).ke: kinetic energy
%  - detection_fields(t).div: divergence
%  - detection_fields(t).vort: vorticity
%  - detection_fields(t).OW: Okubo-Weiss
%  - detection_fields(t).LOW: Local Okubo-Weiss
%  - detection_fields(t).LNAM: Local Normalised Angular Momentum
%  (last 2 fields are calculated using b parameter which could include
%  the interpolation factor 'res')
%
%-------------------------
%  Apr 2015 Briac Le Vu
%-------------------------
%
%=========================

global path_out
global runname
global grid_ll

% No update by default
if nargin==1
    update = 0;
end

% load velocities field (see load_fields for details)
load([path_out,'fields_inter',runname])

% preallocate fields struct array if doesn't exist
stepF = size(u,3);

if update
    step0 = stepF - update+1;
    detection_fields = detection_fields(1:step0-1);
else
    detection_fields = struct('ke',{},'div',{},'vort',{},'OW',{},'LOW',{},'LNAM',{});
    step0 = 1;
end

%----------------------------------------

disp('Compute fields on interpolated grid')

% cycle through time steps
for j=step0:stepF
    
    disp(['  Step ',num2str(j),' ... '])
    
    % load 2D velocity fields for stepj (m/s)
    %----------------------------------------
    uu = squeeze(u(:,:,j));
    vv = squeeze(v(:,:,j));

    % Calculation of eddy kinetic energy
    %----------------------------------------
    ke = (uu.^2 + vv.^2)/2;
    
    % Calculation of finite spatial element
    %----------------------------------------
    dx  = zeros(size(x));
    dy  = dx;
    dux = dx;
    duy = dx;
    dvx = dx;
    dvy = dx;
    
    % Spatial element in deg if grid_ll==1 or in km otherwise
    %----------------------------------------
    dx(2:end-1,2:end-1) = x(2:end-1,3:end) - x(2:end-1,1:end-2); %#ok<*COLND>
    dy(2:end-1,2:end-1) = y(3:end,2:end-1) - y(1:end-2,2:end-1);
    
    if grid_ll
        % define constants
        earth_radius = 6378.137; % km
        % kilometer (km) per degree of latitude
        R = earth_radius*pi/180; % 111.320m
        % Calcul finite spatial element in km
        dx = dx*R.*cosd(y);
        dy = dy*R;
    end
    
    % in meters
    dx = dx*1000; % m
    dy = dy*1000; % m

    % Compute speed element in m/s
    %----------------------------------------
    dux(2:end-1,2:end-1) = (uu(2:end-1,3:end) - uu(2:end-1,1:end-2));
    duy(2:end-1,2:end-1) = (uu(3:end,2:end-1) - uu(1:end-2,2:end-1));
    dvx(2:end-1,2:end-1) = (vv(2:end-1,3:end) - vv(2:end-1,1:end-2));
    dvy(2:end-1,2:end-1) = (vv(3:end,2:end-1) - vv(1:end-2,2:end-1));

    % Calculation of Okubo-Weiss criteria
    %----------------------------------------
    sn = (dux./dx) - (dvy./dy); % shear "cisaillement"
    ss = (dvx./dx) + (duy./dy); % strain "deformation"
    om = (dvx./dx) - (duy./dy); % vorticity "vorticité"

    okubo = sn.^2 + ss.^2 - om.^2; % in s-2

    % Calculation of divergence
    %----------------------------------------
    div = (dux./dx) + (dvy./dy);

    % Calculation of vorticity field (typically +/-10^-5 s-1)
    %----------------------------------------
    vorticity = om;

    % border is a parameter which prevents the constraints
    % to be applied to points too close to the domain boundaries
    % which would result in an index error
    borders = max(b(:)) + 1;
    
    % Calculation of LNAM criteria (Local Normalized Angular Momentum)
    % and LOW criteria (Local Averaged Okubo Weiss)
    %----------------------------------------
    
    % Initialisation
    L   = zeros(size(uu));
    LOW = nan(size(uu));

    % calculate LNAM and LOW in all domain pixels
    %----------------------------------------
    for i=borders:length(vv(:,1))-borders+1
        for ii=borders:length(vv(1,:))-borders+1

            if ~isnan(vv(i,ii))
                
                % calculate LOW
                OW = okubo(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii));
                LOW(i,ii) = mean(OW(:));
                
                % calculate LNAM
                xlocal = x(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % x square sample of length b
                ylocal = y(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % y square sample of length b
                
                ulocal = uu(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % u square sample of length b
                vlocal = vv(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % v square sample of length b

                % Local Normalized Angular Momentum
                % Use the midddle of the square as center coordinate
                
                coordcentre = size(ulocal,1)-b(i,ii);
                
                if grid_ll
                    d_xcentre = (xlocal - xlocal(coordcentre,coordcentre))*R...
                                .*cosd(ylocal(:,:));
                    d_ycentre = (ylocal - ylocal(coordcentre,coordcentre))*R;
                else
                    d_xcentre = (xlocal - xlocal(coordcentre,coordcentre));
                    d_ycentre = (ylocal - ylocal(coordcentre,coordcentre));
                end
                
                cross   = (d_xcentre.*vlocal) - (d_ycentre.*ulocal);
                dot     = (ulocal.*d_xcentre) + (vlocal.*d_ycentre);
                produit = sqrt(ulocal.^2 + vlocal.^2)...
                        .*sqrt(d_xcentre.^2 + d_ycentre.^2);

                if (sum(dot(:))+sum(produit(:))) ~= 0
                    L(i,ii) = sum(cross(:)) / (sum(dot(:)) + sum(produit(:)));
                else
                    L(i,ii) = 0;
                end
                
            end
        end
    end
    
    L(isnan(L)) = 0;

    % save fields in struct array
    %----------------------------------------
    detection_fields(j).ke   = ke.*mask;
    detection_fields(j).div  = div.*mask;
    detection_fields(j).vort = vorticity.*mask;
    detection_fields(j).OW   = okubo.*mask;
    detection_fields(j).LOW  = LOW.*mask;

    % relative value of the LNAM (anticyclone <0 | cyclone >0)
    detection_fields(j).LNAM = L.*mask; 
end

disp(' ')

% export fields
%----------------------------------------
save([path_out,'fields_inter',runname],'detection_fields','-append')
