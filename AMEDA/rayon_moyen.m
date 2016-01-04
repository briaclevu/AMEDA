function [R,A,P,ll] = rayon_moyen(lonlat)
% calcul du rayon moyen (R), du périmètre (P), de l'aire (A)
% d'une structure fermée (lonlat) par rapport à son barycentre (ll)

% taille du contour
lim=size(lonlat,2)-1;

% calcul du barycentre
somme_lon=sum(lonlat(1,1:lim));
somme_lat=sum(lonlat(2,1:lim));
ll(1)=somme_lon/lim;
ll(2)=somme_lat/lim;

% initialisation
distance=zeros(1,lim+1);
aire=zeros(1,lim);
param=zeros(1,lim);

% Distance = distance entre le barycentre et chaque point
longitudes(1)=ll(1);
latitudes(1)=ll(2);
for point=1:lim+1
    longitudes(2)=lonlat(1,point);
    latitudes(2)=lonlat(2,point);

    distance(point)=sw_dist(latitudes,longitudes,'km');
end

% Distance2 = distance entre 2 points consecutifs du contour
distance2=sw_dist(lonlat(2,:),lonlat(1,:),'km');

% calcul du périmètre
P = sum(distance2(:));

% Calcul de l'aire du tourbillon
for point=1:lim
    a=distance(point);
    b=distance(point+1);
    c=distance2(point);
    s=(a+b+c)/2;
    aire(point)=sqrt(s*(s-a)*(s-b)*(s-c));
end
A = sum(real(aire(:)));

% calcul du rayon moyen
for point=1:lim
    rayonmoy=(distance(point)+distance(point+1))/2;
    param(point)=aire(point)*rayonmoy;
end

R(1) = sqrt(A/(4*atan(1))); % R d'un disque d'aire équivalente
R(2) = sum(param(:))/real(A); % rayon pondéré par l'aire de chaque portion
R(3) = max(distance(1:lim)); % rayon max
R(4) = mean(distance(1:lim)); % rayon moyen de chaque point

