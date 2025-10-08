function [pointsInside, values] = process_tetrahedron_gpu(nodes, volume)
    % Compute bounding box of tetrahedron
    minCorner = floor(min(nodes));
    maxCorner = ceil(max(nodes));

    % Make sure bounds are within volume
    sz = size(volume);
    minCorner = max(minCorner, [1,1,1]);
    maxCorner = min(maxCorner, sz);

    % Generate grid of voxel centroids in bounding box
    [X, Y, Z] = ndgrid(minCorner(1):maxCorner(1), ...
                       minCorner(2):maxCorner(2), ...
                       minCorner(3):maxCorner(3));
    % Shift by -0.5 to get voxel centroid
    X = gpuArray(X) - 0.5;
    Y = gpuArray(Y) - 0.5;
    Z = gpuArray(Z) - 0.5;

    % Flatten
    points = [X(:), Y(:), Z(:)];

    % Check which points are inside tetrahedron using barycentric method
    insideMask = point_in_tetrahedron(points, nodes);

    % Get coordinates inside
    pointsInside = points(insideMask, :);

    % Map to voxel indices (undo -0.5 shift, round to nearest voxel)
    voxelIdx = round(pointsInside + 0.5);

    % Remove out-of-bound voxels (can happen at borders)
    validMask = all(voxelIdx >= 1, 2) & ...
                voxelIdx(:,1) <= sz(1) & ...
                voxelIdx(:,2) <= sz(2) & ...
                voxelIdx(:,3) <= sz(3);

    if isempty(validMask)

        centroid = round(mean(nodes, 1));
        values = volume(centroid(1),centroid(2),centroid(3));
        pointsInside = [centroid(1) - 0.5,centroid(2) - 0.5,centroid(3) - 0.5];
    else
        voxelIdx = voxelIdx(validMask, :);
        pointsInside = pointsInside(validMask, :);
    
        % Convert subscripts to linear indices
        linIdx = sub2ind(sz, voxelIdx(:,1), voxelIdx(:,2), voxelIdx(:,3));
    
        % Get voxel values
        values = volume(linIdx);
        

    end


end

