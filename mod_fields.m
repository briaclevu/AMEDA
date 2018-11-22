function fields = mod_fields(source,stp,resolution)
%fields = mod_fields(source,stp {,resolution})
%
%  Creates 2D fields of various computation by finite spatial element
%  for the step 'stp' of the data of type 'source' (AVISO, NEMO,...) that
%  determine which load_field_'source' is  used in the routine.
%
%  Use velocities get by load_field_'source'.m including a proper mask
%  with the first b pixel of land (must include the interpolation factor
%  'resolution' if any) and NaN in other land pixel.
%
%  For a description of the input parameters see mod_eddy_param.m.

%
%  Computed 2-Dfields are saved/update in [path_out,'fields_',runname] as:
%  - fields.ke: kinetic energy
%  - fields.div: divergence
%  - fields.vort: vorticity
%  - fields.OW: Okubo-Weiss
%  - fields.LOW: Local Okubo-Weiss
%  - fields.LNAM: Local Normalised Angular Momentum
%  (last 2 fields are calculated using b parameter which could include
%  the interpolation factor 'resol')
%
%-------------------------
%  June 2016 Briac Le Vu
%-------------------------
%
%=========================

disp([' <<< Compute fields step ',num2str(stp)])

%----------------------------------------
% load keys_sources and parameters (use mod_eddy_params.m first)
load('param_eddy_tracking','grid_ll','resol','deg','b','bi','f','f_i')

%----------------------------------------
% replace parameters by arguments
if nargin==3
    resol = resolution;
end

%----------------------------------------
% update b parameter and define I/O matfile
if resol>1
    b=bi;
    f=f_i;
end

%----------------------------------------
% load 2D velocity fields (m/s) for step stp
eval(['[x,y,mask,uu,vv,~] = load_fields_',source,'(stp,resol,deg);'])

%----------------------------------------
% Calculation of eddy kinetic energy
ke = (uu.^2 + vv.^2)/2;

%----------------------------------------
% Calculation of finite spatial element
dx  = zeros(size(x));
dy  = dx;
dux = dx;
duy = dx;
dvx = dx;
dvy = dx;

%----------------------------------------
% Spatial element in deg if grid_ll==1 or in km otherwise
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

%----------------------------------------
% Compute speed element in m/s
dux(2:end-1,2:end-1) = (uu(2:end-1,3:end) - uu(2:end-1,1:end-2));
duy(2:end-1,2:end-1) = (uu(3:end,2:end-1) - uu(1:end-2,2:end-1));
dvx(2:end-1,2:end-1) = (vv(2:end-1,3:end) - vv(2:end-1,1:end-2));
dvy(2:end-1,2:end-1) = (vv(3:end,2:end-1) - vv(1:end-2,2:end-1));

%----------------------------------------
% Calculation of Okubo-Weiss criteria
sn = (dux./dx) - (dvy./dy); % shear "cisaillement"
ss = (dvx./dx) + (duy./dy); % strain "deformation"
om = (dvx./dx) - (duy./dy); % vorticity "vorticité"

okubo = sn.^2 + ss.^2 - om.^2; % in s-2

%----------------------------------------
% Calculation of divergence
div = (dux./dx) + (dvy./dy);

%----------------------------------------
% Calculation of vorticity field (typically +/-10^-5 s-1)
vorticity = om.*sign(f);

%----------------------------------------
% border is a parameter which prevents the constraints
% to be applied to points too close to the domain boundaries
% which would result in an index error
borders = max(b(:)) + 1;

%----------------------------------------
% Calculation of LNAM criteria (Local Normalized Angular Momentum)
% and LOW criteria (Local Averaged Okubo Weiss)
%----------------------------------------

disp('Compute LNAM >>>')

% Initialisation
L   = zeros(size(uu));
LOW = nan(size(uu));

%----------------------------------------
% calculate LNAM and LOW in all domain pixels
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
            sumdp = sum(dot(:))+sum(produit(:));

            if sumdp ~= 0
                L(i,ii) = sum(cross(:)) / sumdp * sign(f(i,ii));
            else
                L(i,ii) = 0;
            end

        end
    end
end
    
L(isnan(L)) = 0;

%----------------------------------------
% export fields in struct array
fields.step = stp;
fields.ke   = ke.*mask;
fields.div  = div.*mask;
fields.vort = vorticity.*mask;
fields.OW   = okubo.*mask;
fields.LOW  = LOW.*mask;

% relative value of the LNAM (anticyclone <0 | cyclone >0)
fields.LNAM = L.*mask;

disp(' ')
