function [inside_points, inside_vals] = pointsInsideTetrahedronScaled(P1, P2, P3, P4, scale, ph)
% Vectorized version – MUCH faster than triple for-loop
% Inputs: 4 vertices (3x1 each), scale factor, phantom volume
% Outputs: points inside tetrahedron + phantom values

    % Scale vertices
    P1s = P1 * scale; P2s = P2 * scale; P3s = P3 * scale; P4s = P4 * scale;

    % Build edge matrix
    u = P2s - P1s;
    v = P3s - P1s;
    w = P4s - P1s;
    M = [u, v, w];

    % Axis-aligned bounding box
    vertices = [P1s, P2s, P3s, P4s];
    min_corner = floor(min(vertices, [], 2));
    max_corner = ceil(max(vertices, [], 2));

    % Generate all voxel centers in bounding box
    [X,Y,Z] = ndgrid(min_corner(1):max_corner(1), ...
                     min_corner(2):max_corner(2), ...
                     min_corner(3):max_corner(3));
    pts = [X(:) Y(:) Z(:)] + (-0.5);  % voxel centers

    % Compute barycentric coords (vectorized)
    rel_pts = bsxfun(@minus, pts, P1s');  % subtract P1 from all points
    params = (M \ rel_pts')';             % Nx3 barycentric coords
    alpha = params(:,1); beta = params(:,2); gamma = params(:,3);

    insideMask = (alpha >= 0) & (beta >= 0) & (gamma >= 0) & ...
                 (alpha + beta + gamma <= 1);

    % Select inside points
    inside_points_scaled = pts(insideMask, :);

    % Phantom values (need integer indices!)
    idx = round(inside_points_scaled);
    validMask = all(idx >= 1, 2) & ...
                idx(:,1) <= size(ph,1) & ...
                idx(:,2) <= size(ph,2) & ...
                idx(:,3) <= size(ph,3);

    idx = idx(validMask, :);

    inside_vals = ph(sub2ind(size(ph), idx(:,1), idx(:,2), idx(:,3)));
    inside_points = inside_points_scaled(validMask, :) / scale;
end

