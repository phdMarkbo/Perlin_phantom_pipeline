function tetraVoxelMap = tetra_to_voxel_mapping(voxelVolume, tetMesh, nodeCoords, tetrahedron_index)

    % Parameters

    voxelVolume = gpuArray(voxelVolume);  % Random voxel values on GPU

    % Generate a simple tetrahedral mesh (example: divide cube into tetrahedra)
    
    numTets = size(tetMesh, 1);
    % Output structure
    tetraVoxelMap = cell(size(tetMesh, 1), 1);

    % Loop through each tetrahedron
    
    Displayer = Progressdisp(numTets, BarOrStr="bar"); % Initialize
    for i = 1:numTets
        Displayer.Display_progress(i); % Update progress dynamically
        nodes = nodeCoords(tetMesh(i, :), :);  % 4x3 matrix
        [pointsInside, valuesInside] = process_tetrahedron_gpu(nodes, voxelVolume);
        
        if isempty(valuesInside)
            avgVal = NaN;
        else
            avgVal = mean(gather(valuesInside));
        end

        tetraVoxelMap{i} = struct( ...
            'coords', gather(pointsInside), ...
            'values', gather(valuesInside), ...
            'avgValue', avgVal ...
        );

    end
    if tetrahedron_index > 0 && tetrahedron_index <= size(tetMesh, 1)
        plot_tetrahedron_with_internal_points_voxel_values(tetraVoxelMap, tetMesh, nodeCoords, tetrahedron_index);
    elseif tetrahedron_index == 0
        disp("Not plotting")
    else
        disp("Not a viable tetrahedron index for plotting. Please only use integers between 1 and " + num2str(size(tetMesh, 1)) + ". Use 0 to omit plot.")
    end
end

