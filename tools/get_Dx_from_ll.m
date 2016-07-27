function Dx = get_Dx_from_ll(x,y)
% compute the distance matrice from lon and lat matrice
%
% June 2016 B. LE VU
%
disty = sw_dist( y(1:2,1), x(1:2,1), 'km');
distx = nan*x;
for i=2:size(x,1)-1
    for j=2:size(x,2)-1
        distx(i,j) = sw_dist( y(i,[j-1 j+1]), x(i,[j-1 j+1]), 'km')/2;
    end
end
Dx = ( distx + disty )/2;
Dx(1,:)  = Dx(2,:)-diff(Dx(2:3,:));
Dx(end,:)= Dx(end-1,:)+diff(Dx(end-2:end-1,:));
Dx(:,1)  = Dx(:,2);
Dx(:,end)= Dx(:,end-1);
