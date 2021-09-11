function y = homography_transform(x, v)
% HOMOGRAPHY_TRANSFORM applies homographic transform to vectors
%   Y = HOMOGRAPHY_TRANSFORM(X, V) takes a 2xN matrix, each column of which
%   gives the position of a point in a plane. It returns a 2xN matrix whose
%   columns are the input vectors transformed according to the homography
%   V, represented as a 3x3 homogeneous matrix.

% dummy = ones(1, size(x,2));
% dummy1 = [x;dummy];
q = v * x;
p = q(3,:);
y = [q(1,:)./p; q(2,:)./p];
end