function [C,P] = compute_curve(xy,Np,grid_ll)
% calculate curvature for each point of contour xy using fix number of grid
% point
% x = xy(1,:) and y = xy(2,:)
%

% number of total point
N = size(xy,2) - 1;

% adjust Np
if N < Np*2+1
    Np = fix((N-1)/2);
end

% average mean of coordinates
Mx = mean(xy(1,1:N));
My = mean(xy(2,1:N));

% duplicate the contour xy
if xy(1,1)==xy(1,end) && xy(2,1)==xy(2,end)
    xy2=[xy(:,1:end-1) xy xy(:,2:end)];
else
    xy2=[xy xy xy];
end

% segment length along contour xy
if grid_ll
    P = sw_dist2(xy2(2,:),xy2(1,:),'km');
else
    P = sqrt(diff(xy2(2,:)).^2 + diff(xy2(1,:)).^2); % km
end

% initialize curvature
C = zeros(1,N);

% sliding curvature for Np-point segment
for i = 1:N
    
    % find and extract the ith segment
    x = xy2(1,N+i-Np:N+i+Np);
    y = xy2(2,N+i-Np:N+i+Np);
    
    % coordinate of the segment in km
    if grid_ll

        xs(1) = Mx;
        ys(1) = My;

        % initialise
        coord = zeros(2,length(x));

        for pt=1:length(x)

            xs(2) = x(pt);
            ys(2) = y(pt);

            % distances in km of every point from the barycenter
            coord(1,pt) = sign(diff(xs)) * sw_dist2([My My],xs,'km');% dx
            coord(2,pt) = sign(diff(ys)) * sw_dist2(ys,[Mx Mx],'km');% dy
        end

    else
        coord = [x(1:length(x));y(1:length(x))];
    end
    
    try
        % fit a circle on an arc with Taubin (1991)
        ParIn = TaubinNTN(coord');
        R = ParIn(3);
        
        % geometrical fitting (twice coastly)
        %Par = LM([x',y'],ParIn);
        %R = Par(3);
        
        % In or ~IN the xy contour
        IN = inpolygon(mean(x),mean(y),xy(1,:),xy(2,:));
        if IN
            C(i) = 1/R;
        else
            C(i) = -1/R;
        end
    catch
        C(i) = NaN;
    end
    
end

% trunc P
P = P(1:N);

