function [lines,lvl] = scan_lines(C)
%[lines,lvl] = scan_lines(C)
%
% Rearrange all the streamlines in C into the structure array 'lines'
% sorted by maximum y coord. Each element of isolines contains
% all the vertices of a given contour level(lvl) in C
%
%  OUTPUT:
%  - lines is a structure array with x, y and max(y) for each streamlines
%  - lvl is the value of the streamline
%
%-------------------------
%   June 2016 by B. LE VU
%-------------------------
%
%=========================

% fill the structure 'isolines'
lines = struct('x',{},'y',{},'l',{});
lvl=[];

% begin two counters
k = 1;
kk = 1;

while k < size(C,2)
    npoints = C(2,k);
    lvl(kk) = C(1,k);
    lines(kk).x = C(1,k+1:k+npoints); % vertex x's
    lines(kk).y = C(2,k+1:k+npoints); % vertex y's
    lines(kk).l = max(C(2,k+1:k+npoints)); % max y of a curve
    kk = kk + 1;
    k = k + npoints + 1;
end

% sort the contours according to their maximum y coord; this way the first
% closed contour across which velocity increases will also be the largest
% one (it's the one which extend further north).
[~,order] = sort([lines(:).l],'ascend');
lines = lines(order);
lvl = lvl(order);
% ! Debug ! Test the contour value scanned
%display([min(diff(lvl)) mean(diff(lvl)) max(diff(lvl))]) % !Debug!

