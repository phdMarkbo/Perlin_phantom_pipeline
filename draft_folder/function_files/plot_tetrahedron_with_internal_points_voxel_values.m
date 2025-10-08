function plot_tetrahedron_with_internal_points_voxel_values(tetraVoxelMap, tetMesh, nodeCoords, tet_ind)
    figure;
    hold on;
    axis equal;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Tetrahedron and Internal Voxel Centroids Colored by Voxel Values');

    % Get first tetrahedron
    tetNodes = nodeCoords(tetMesh(tet_ind,:), :);
    internalPoints = tetraVoxelMap{tet_ind}.coords;
    voxelValues = tetraVoxelMap{tet_ind}.values;

    % Plot the tetrahedron as a semi-transparent 3D surface
    faces = [
        1 2 3;
        1 2 4;
        1 3 4;
        2 3 4
    ];

    patch('Faces', faces, 'Vertices', tetNodes, ...
        'FaceColor', 'cyan', 'FaceAlpha', 0.2, 'EdgeColor', 'k');

    % Normalize voxel values for coloring (optional)
    voxelValues = double(voxelValues);
    cmin = min(voxelValues);
    cmax = max(voxelValues);


    if cmin == cmax
        % Plot with a fixed color
        scatter3(internalPoints(:,1), internalPoints(:,2), internalPoints(:,3), ...
                 50, 'r', 'filled');
        legend('Tetrahedron', sprintf('Voxel centroids (value = %.3f)', cmin));
    else
        % Normalize and plot with color
        colors = (voxelValues - cmin) / (cmax - cmin);
        scatter3(internalPoints(:,1), internalPoints(:,2), internalPoints(:,3), ...
                 50, colors, 'filled');
        colormap('jet');
        colorbar;
        caxis([cmin cmax]);
        legend('Tetrahedron', 'Voxel centroids (colored by value)');
    end

end

