function L = getRotMat(orientation, elemNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return rotation matrix, L
%
% Input: 1) orientation - # elems x 9 (a, b, f)
%        2) elemNum - current element number
%
% Output: L - rotation matrix
%
% Renee Miller
% 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get local coordinate system z' axis in current element
idx = find(orientation(:,1)==elemNum);

% Get rotation matrix (changing coordinate system)
localXY_plane = [0 0 0 orientation(idx,2:4) orientation(idx,5:7)]; % Origin + two points that define the x-y plane in the local coordinate system
tmp = createBasisTransform3d('g',localXY_plane); % Function finds transformation matrix global x-y plane and local material x-y plane
L = tmp(1:3,1:3);