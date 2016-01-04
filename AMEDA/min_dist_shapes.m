function ds=min_dist_shapes(lonlat1,lonlat2) 
% calcul the distance between 2 eddy edge lonlat1 & lonlat2
ds=nan(length(lonlat1),length(lonlat2));
for i=1:length(lonlat1)
    for j=1:length(lonlat2)
	latitudes=[lonlat1(2,i),lonlat2(2,j)];
	longitudes=[lonlat1(1,i),lonlat2(1,j)];
        ds(i,j)=sw_dist(latitudes,longitudes,'km');
    end
end
ds=nanmin(ds(:));
