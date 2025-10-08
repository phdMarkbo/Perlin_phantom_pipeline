tet_ind_to_be_plotted = 21680;

phantom_struct = load("PerlinPhantom01.mat");
elems = load("element.dat");
nodes_orig = load("node.dat");
nodes_comp = readtable("cupA_01_all_materials_compressed_nodes.txt");

nodes_comp = table2array(nodes_comp);
finalVolume = phantom_struct.finalVolume;

finalVolume(finalVolume==8) = 7;
finalVolume(finalVolume==9) = 8;
finalVolume(finalVolume==11) = 9;
finalVolume(finalVolume==12) = 10;
finalVolume(finalVolume==13) = 11;

nodes_orig = nodes_orig(:, 2:4) * 0.001;
nodes_comp = nodes_comp(1:29479, 2:4);
elems = elems(:, 2:5);

test = nodes_orig * 10000 + 1;

%%
tic;
mapping_struct = tetra_to_voxel_mapping(finalVolume, elems, test,tet_ind_to_be_plotted);
elapsedTime = toc;

fprintf('Time spent calculating: %.2f seconds (%.2f minutes)\n', ...
        elapsedTime, elapsedTime/60);
%mapping_struct_gpu = tetra_to_voxel_mapping(finalVolume, elems, test,tet_ind_to_be_plotted);
%%
info_array = uniqueCountsPerStruct(mapping_struct, "values", nodes_orig, elems);
%%
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

function inside = point_in_tetrahedron(P, T)
    % Barycentric coordinate test for points inside a tetrahedron
    v0 = T(2,:) - T(1,:);
    v1 = T(3,:) - T(1,:);
    v2 = T(4,:) - T(1,:);
    A = [v0', v1', v2'];  % 3x3

    relP = (P - T(1,:))';  % 3xN
    baryCoords = A \ relP;  % 3xN

    u = baryCoords(1,:)';
    v = baryCoords(2,:)';
    w = baryCoords(3,:)';
    t = 1 - u - v - w;

    % Inside test: all barycentric coords in [0,1]
    inside = (u >= 0) & (v >= 0) & (w >= 0) & (t >= 0);
end

function plot_tetrahedron_with_internal_points(tetraVoxelMap, tetMesh, nodeCoords)
    figure;
    hold on;
    axis equal;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Tetrahedron and Internal Voxel Centroids');

    % Get first tetrahedron
    tetNodes = nodeCoords(tetMesh(1,:), :);
    internalPoints = tetraVoxelMap{1}.coords;

    % Plot the tetrahedron as a patch3
    % Each face of the tetrahedron is a triangle
    faces = [
        1 2 3;
        1 2 4;
        1 3 4;
        2 3 4
    ];

    patch('Faces', faces, 'Vertices', tetNodes, ...
        'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'k');

    % Plot internal points
    scatter3(internalPoints(:,1), internalPoints(:,2), internalPoints(:,3), ...
             36, 'r', 'filled');

    legend('Tetrahedron', 'Internal voxel centroids');
end
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

function [pointsInside, values, integralEstimate, meanEstimate] = ...
         process_tetrahedron_gpu_riemann(nodes, volume, samplingStep, useCoarseGrid)
    % PROCESS_TETRAHEDRON_GPU_RIEMANN
    % nodes: 4x3 tetra vertex coordinates (x,y,z)
    % volume: 3D array with voxel values (assumed unit voxel spacing)
    % samplingStep: scalar integer >=1 (equidistant sampling)
    % useCoarseGrid: (optional) true -> shift = -samplingStep/2 (coarse grid),
    %                         false -> shift = -0.5 (true voxel centers)
    %
    % Returns sampled points, their values, an integral estimate, and mean.

    if nargin < 3 || isempty(samplingStep)
        samplingStep = 1;
    end
    if nargin < 4 || isempty(useCoarseGrid)
        useCoarseGrid = true; % recommended
    end

    % Compute bounding box of tetrahedron (in voxel coordinates)
    minCorner = floor(min(nodes, [], 1));
    maxCorner = ceil(max(nodes, [], 1));

    % Limit to volume bounds
    sz = size(volume);
    minCorner = max(minCorner, [1,1,1]);
    maxCorner = min(maxCorner, sz);

    % Choose shift for centering sampled cells
    if useCoarseGrid
        shift = samplingStep/2;      % centers of coarser cells
    else
        shift = 0.5;                 % true voxel centers are at n-0.5
    end

    % Build sampled grid coordinates (note using :samplingStep:)
    [X, Y, Z] = ndgrid(minCorner(1):samplingStep:maxCorner(1), ...
                       minCorner(2):samplingStep:maxCorner(2), ...
                       minCorner(3):samplingStep:maxCorner(3));

    % Move to centroids: original voxel centroid offset is -0.5,
    % so we subtract shift to place sample points
    X = gpuArray(X) - shift;
    Y = gpuArray(Y) - shift;
    Z = gpuArray(Z) - shift;

    % Flatten and move to CPU for point-in-tetra test if necessary
    points = [X(:), Y(:), Z(:)];   % still gpuArray columns

    % Check which points are inside tetrahedron
    insideMask = point_in_tetrahedron(points, nodes); % should accept gpu arrays
    pointsInside = points(insideMask, :);

    % Map sampled points back to voxel indices to read values
    % Undo shift: if shift == 0.5 this maps to true voxel indices; if shift==samplingStep/2
    % mapping will pick the appropriate coarse cell representative index.
    voxelIdx = round(pointsInside + shift); % inverse of subtracting shift

    % Remove out-of-bounds
    validMask = all(voxelIdx >= 1, 2) & ...
                voxelIdx(:,1) <= sz(1) & ...
                voxelIdx(:,2) <= sz(2) & ...
                voxelIdx(:,3) <= sz(3);

    if isempty(pointsInside) || ~any(validMask)
        % Fallback: evaluate at centroid voxel (weight = 1)
        centroid = round(mean(nodes, 1));
        centroid = max(min(centroid, sz), [1,1,1]);
        values = volume(centroid(1), centroid(2), centroid(3));
        pointsInside = [centroid(1) - 0.5, centroid(2) - 0.5, centroid(3) - 0.5];
        % Integral estimate: treat single voxel with unit volume
        integralEstimate = double(values) * 1.0;
        % tetrahedron geometric volume for normalization
        tetVol = abs(det([nodes(2,:)-nodes(1,:); nodes(3,:)-nodes(1,:); nodes(4,:)-nodes(1,:)]) ) / 6;
        meanEstimate = integralEstimate / tetVol;
        return;
    end

    voxelIdx = voxelIdx(validMask, :);
    pointsInside = pointsInside(validMask, :);

    % Convert to linear indices and fetch values
    linIdx = sub2ind(sz, voxelIdx(:,1), voxelIdx(:,2), voxelIdx(:,3));
    values = double(gather(volume(linIdx))); % bring values to CPU and double for sum

    % Riemann weight: each sampled point represents a cube of side samplingStep
    sampleVolume = (samplingStep)^3;

    % Integral estimate = sum(values * sampleVolume)
    integralEstimate = sum(values) * sampleVolume;

    % Geometric volume of tetrahedron (for mean)
    tetVol = abs(det([nodes(2,:)-nodes(1,:); nodes(3,:)-nodes(1,:); nodes(4,:)-nodes(1,:)]) ) / 6;
    if tetVol == 0
        meanEstimate = NaN;
    else
        meanEstimate = integralEstimate / tetVol;
    end
end

function out = uniqueCountsPerStruct(C, fieldname, nodes, elems)
    % C: Nx1 cell array of structs
    % fieldname: string, the field to extract (e.g. 'vals')
    % nodes: Mx3 array of node coordinates
    % elems: Nx4 array of tetrahedral elements (indices into nodes)
    % out: Nx1 struct array, each element has .unique, .counts, .volume, .totalCount

    N = numel(C);
    out(N) = struct('unique', [], 'counts', [], 'volume', [], 'totalCount', []);  % preallocate

    for i = 1:N
        % --- Unique values and counts ---
        vals = C{i}.(fieldname);       % extract values from struct
        [u,~,idx] = unique(vals);      % unique values + mapping
        counts = accumarray(idx(:), 1);% count occurrences

        % --- Volume of tetrahedron ---
        elemNodes = nodes(elems(i,:), :);  % 4x3 array of node coordinates
        % Formula: V = |dot(a, cross(b,c))| / 6
        a = elemNodes(2,:) - elemNodes(1,:);
        b = elemNodes(3,:) - elemNodes(1,:);
        c = elemNodes(4,:) - elemNodes(1,:);
        vol = abs(dot(a, cross(b,c))) / 6;

        % --- Store results in struct array ---
        out(i).unique = u;
        out(i).counts = counts;
        out(i).volume = vol;
        out(i).totalCount = sum(counts);  % sum of counts
    end
end
