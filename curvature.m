function [C,P] = curvature(xy,Rd,grid_ll)
% calculate curvature for each point of contour xy using segment of Rd size
% x = xy(1,:) and y = xy(2,:)

% number of point
N = size(xy,2)-1;

% duplicate the contour xy
if xy(1,1)==xy(1,end) && xy(2,1)==xy(2,end)
    xy2=[xy(:,1:end-1) xy xy(:,2:end)];
else
    xy2=[xy xy xy];
end

% segment length along contour xy
if grid_ll
    P = sw_dist(xy2(2,:),xy2(1,:),'km');
else
    P = sqrt(diff(xy2(2,:)).^2 + diff(xy2(1,:)).^2); % km
end

% sliding curvature for Rd-length segment
C = zeros(1,N);
for i = 1:N
    % cumulative sum back and forward at the ith order
    Pcum = cumsum(P(i:end))+cumsum(P(end-i+1:-1:1));
    % find and extract the ith segment
    ind = min(N,find(Pcum>Rd,1));
    x = xy2(1,N+i-ind:N+i+ind);
    y = xy2(2,N+i-ind:N+i+ind);
    % fit a circle on an arc with Taubin (1991)
    ParIn = TaubinNTN([x',y']);
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
end

% trunc P
P = P(1:N);
    

