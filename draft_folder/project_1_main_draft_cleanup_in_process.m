addpath("/home/user/code/Perlin_phantom_pipeline/draft_folder/function_files")
addpath("/home/user/code/Perlin_phantom_pipeline/draft_folder/load_files/")
addpath("/home/user/code/Perlin_phantom_pipeline/Progressdisp.m/")

tet_ind_to_be_plotted = 1;

phantom_struct = load("PerlinPhantom01.mat");
elems = load("element.dat");
nodes_orig = load("node.dat");
nodes_comp = readtable("cupA_01_all_materials_compressed_nodes.txt");
nodes_comp = table2array(nodes_comp);

finalVolume = phantom_struct.finalVolume;
clear phantom_struct;

finalVolume = replaceValues(finalVolume, [8, 7; 9, 8; 11, 9; 12, 10; 13,11]);

% finalVolume(finalVolume==8) = 7;
% finalVolume(finalVolume==9) = 8;
% finalVolume(finalVolume==11) = 9;
% finalVolume(finalVolume==12) = 10;
% finalVolume(finalVolume==13) = 11;

nodes_orig = nodes_orig(:, 2:4) * 0.001;
nodes_comp = nodes_comp(1:29479, 2:4);
elems = elems(:, 2:5);

test = nodes_orig * 10000 + 1;

%%
%%%SCRIPT1

mapping_struct = tetra_to_voxel_mapping(finalVolume, elems, test,tet_ind_to_be_plotted);

%%
%SCRIPT2

[affine_coords, inside_vals] = pointsInsideTetrahedraByIndices(nodes_orig, nodes_comp, elems, 10000, false, finalVolume);

%%
%SCRIPT3
% 
% clear;
% disp("loading")
% affine_coords = load("inside_ps_coords.mat");
% inside_vals = load("inside_ps_vals.mat");
% affine_coords = affine_coords.inside_ps;
% inside_vals = inside_vals.inside_val;
tic;
compressed_phantom1 = phantom_from_transformed_points(affine_coords, inside_vals);
elapsed_time = toc;
fprintf('script 3 took %.2f seconds.\n', elapsed_time);

%%
%function filled_phantom = phantom_from_transformed_points(affine_coords, inside_vals)
%    affine_coords_nz = (affine_coords - min(affine_coords))*10000;
%    
%    integer_coords_nz = uint16(round(affine_coords_nz)) + 1;
%    
%    max_x = max(integer_coords_nz(:, 1));
%    max_y = max(integer_coords_nz(:, 2));
%    max_z = max(integer_coords_nz(:, 3));
%    bounding_box = false(max_y, max_x, max_z);
%    
%    for i = 1:length(integer_coords_nz)
%    
%        cur_X = integer_coords_nz(i, 1);
%        cur_Y = integer_coords_nz(i, 2);
%        cur_Z = integer_coords_nz(i, 3);
%     
%        bounding_box(cur_Y, cur_X, cur_Z) = true;
%    end
%    
%    strel = ones(50, 50, 50);
%    
%    filled_bounding_vol = imclose(bounding_box,strel);
%    filled_bounding_vol_uint16 = uint16(filled_bounding_vol);
%    filled_bounding_vol_uint16(filled_bounding_vol_uint16 == 0) = intmax("uint16");
%    filled_bounding_vol_uint16(filled_bounding_vol_uint16 == 1) = 0;
%    
%    for i = 1:size(integer_coords_nz,1)
%        
%        cur_X = integer_coords_nz(i, 1);
%        cur_Y = integer_coords_nz(i, 2);
%        cur_Z = integer_coords_nz(i, 3);
%        filled_bounding_vol_uint16(cur_Y, cur_X, cur_Z) = inside_vals(i);
%    
%    end
%    
%    nonzeromask = filled_bounding_vol_uint16 ~= 0; %All breast voxels with assigned values
%    [~, idxClosest] = bwdist(nonzeromask);
%    filled_phantom = filled_bounding_vol_uint16(idxClosest);
%    filled_phantom(filled_phantom==intmax("uint16")) = 0;
%    filled_phantom = filled_phantom .* uint16(filled_bounding_vol); %%mask it
%    %masked_Perlin_phantom = finalVolume(1:750, 1:1700, 1:374) .* uint16(filled_bounding_vol);
%end
%%%
%function [all_points, all_vals] = pointsInsideTetrahedraByIndices(nodes,nodes_trans, tets, scale, doPlot, ph)
%    % pointsInsideTetrahedraByIndices
%    % Inputs:
%    %   nodes - Mx3 array of node coordinates
%    %   tets  - Nx4 array of tetrahedron node indices (1-based)
%    %   scale - scale factor for finer grid resolution
%    %   doPlot - logical flag to plot each tetrahedron
%    % Output:
%    %   all_points - Px3 matrix of all interior points
%
%    %nbr_of_phantom_voxels = size(ph(:), 1);
%
%    all_points = zeros(300000000, 3); %allocating space %%REPLACE WITH SIZE FROM mapped_Struct
%    all_vals = zeros(300000000, 1); %allocating space
%
%    numTets = size(tets, 1);
%    all_points_idx = 1;
%
%    Displayer1 = Progressdisp(numTets, BarOrStr="bar"); % Initialize progress bar
%
%    for i = 1:numTets
%
%        Displayer1.Display_progress(i); % Show progress bar
%
%        idx = tets(i, :);  % Get node indices for tetrahedron i
%        
%        % Extract coordinates of the 4 nodes (transpose to get 3x1 vectors)
%        P1 = nodes(idx(1), :)';
%        P2 = nodes(idx(2), :)';
%        P3 = nodes(idx(3), :)';
%        P4 = nodes(idx(4), :)';
%
%        P1_trans = nodes_trans(idx(1), :)';
%        P2_trans = nodes_trans(idx(2), :)';
%        P3_trans = nodes_trans(idx(3), :)';
%        P4_trans = nodes_trans(idx(4), :)';
%        
%        tet_ref = [P1'; P2'; P3'; P4'];
%        tet_trans = [P1_trans'; P2_trans'; P3_trans'; P4_trans'];
%
%        % Compute interior points
%        [points, vals] = pointsInsideTetrahedronScaled(P1, P2, P3, P4, scale, ph);
%
%        if isempty(points)
%
%            %points = mean(tet_ref);
%            continue
%        end
%
%        % Compute affine transformation
%        points_trans = general_tet_mapping(tet_ref, tet_trans, points')';
%
%        nbr_of_points_in_tetrahedra = size(points, 1);
%
%        start_idx = all_points_idx;
%        end_idx = all_points_idx + nbr_of_points_in_tetrahedra - 1;
%
%        all_points(start_idx:end_idx, :) = points_trans;
%        all_vals(start_idx:end_idx, :) = vals;
%
%        all_points_idx = all_points_idx + nbr_of_points_in_tetrahedra;
%
%        % Plot if requested
%        if doPlot
%            plotTetrahedronWithPoints(points, P1, P2, P3, P4);
%            title(['Tetrahedron_orig #' num2str(i)]);
%        end
%        if doPlot
%            plotTetrahedronWithPoints(points_trans, P1_trans, P2_trans, P3_trans, P4_trans);
%            title(['Tetrahedron #' num2str(i)]);
%        end
%    end
%
%    all_points = all_points(1:end_idx, :); %crop allocated space at the end of nonzero entries
%    all_vals = all_vals(1:end_idx, :);
%
%end
%function tet2_coords = general_tet_mapping(T1, T2, tet1_coords)
%
%    %tet1 ===================
%
%    a1x = T1(1,1); %vertice 1
%    a1y = T1(1,2);
%    a1z = T1(1,3);
%
%    b1x = T1(2,1); %vertice 2
%    b1y = T1(2,2);
%    b1z = T1(2,3);
%
%    c1x = T1(3,1); %vertice 3
%    c1y = T1(3,2);
%    c1z = T1(3,3);
%
%    d1x = T1(4,1); %vertice 4
%    d1y = T1(4,2);
%    d1z = T1(4,3);
%
%    %tet2 ====================
%
%    a2x = T2(1,1); %vertice 1
%    a2y = T2(1,2);
%    a2z = T2(1,3);
%
%    b2x = T2(2,1); %vertice 2
%    b2y = T2(2,2);
%    b2z = T2(2,3);
%
%    c2x = T2(3,1); %vertice 3
%    c2y = T2(3,2);
%    c2z = T2(3,3);
%
%    d2x = T2(4,1); %vertice 4
%    d2y = T2(4,2);
%    d2z = T2(4,3);
%
%    D1 = [d1x; d1y; d1z];
%    D2 = [d2x; d2y; d2z];
%
%    A1 = [a1x b1x c1x; 
%          a1y b1y c1y; 
%          a1z b1z c1z] - D1;
%
%    A2 = [a2x b2x c2x; 
%          a2y b2y c2y; 
%          a2z b2z c2z] - D2;
%    
%    tet2_coords = A2 * (A1 \ (tet1_coords - D1)) + D2;
%
%end
%
%function [inside_points, inside_vals] = pointsInsideTetrahedronScaled(P1, P2, P3, P4, scale, ph)
%% Vectorized version – MUCH faster than triple for-loop
%% Inputs: 4 vertices (3x1 each), scale factor, phantom volume
%% Outputs: points inside tetrahedron + phantom values
%
%    % Scale vertices
%    P1s = P1 * scale; P2s = P2 * scale; P3s = P3 * scale; P4s = P4 * scale;
%
%    % Build edge matrix
%    u = P2s - P1s;
%    v = P3s - P1s;
%    w = P4s - P1s;
%    M = [u, v, w];
%
%    % Axis-aligned bounding box
%    vertices = [P1s, P2s, P3s, P4s];
%    min_corner = floor(min(vertices, [], 2));
%    max_corner = ceil(max(vertices, [], 2));
%
%    % Generate all voxel centers in bounding box
%    [X,Y,Z] = ndgrid(min_corner(1):max_corner(1), ...
%                     min_corner(2):max_corner(2), ...
%                     min_corner(3):max_corner(3));
%    pts = [X(:) Y(:) Z(:)] + (-0.5);  % voxel centers
%
%    % Compute barycentric coords (vectorized)
%    rel_pts = bsxfun(@minus, pts, P1s');  % subtract P1 from all points
%    params = (M \ rel_pts')';             % Nx3 barycentric coords
%    alpha = params(:,1); beta = params(:,2); gamma = params(:,3);
%
%    insideMask = (alpha >= 0) & (beta >= 0) & (gamma >= 0) & ...
%                 (alpha + beta + gamma <= 1);
%
%    % Select inside points
%    inside_points_scaled = pts(insideMask, :);
%
%    % Phantom values (need integer indices!)
%    idx = round(inside_points_scaled);
%    validMask = all(idx >= 1, 2) & ...
%                idx(:,1) <= size(ph,1) & ...
%                idx(:,2) <= size(ph,2) & ...
%                idx(:,3) <= size(ph,3);
%
%    idx = idx(validMask, :);
%
%    inside_vals = ph(sub2ind(size(ph), idx(:,1), idx(:,2), idx(:,3)));
%    inside_points = inside_points_scaled(validMask, :) / scale;
%end
%
%
%function plotTetrahedronWithPoints(points, P1, P2, P3, P4)
%    % plotTetrahedronWithPoints
%    % Inputs:
%    %   points - Nx3 matrix of points inside tetrahedron
%    %   P1, P2, P3, P4 - 3x1 vectors of tetrahedron vertices
%    
%    vertices = [P1, P2, P3, P4];
%    
%    edges = [
%        1 2;
%        1 3;
%        1 4;
%        2 3;
%        2 4;
%        3 4];
%    
%    figure;
%    hold on;
%    scatter3(points(:,1), points(:,2), points(:,3), 36, 'filled');
%    
%    for i = 1:size(edges,1)
%        p_start = vertices(:, edges(i,1));
%        p_end = vertices(:, edges(i,2));
%        plot3([p_start(1), p_end(1)], ...
%              [p_start(2), p_end(2)], ...
%              [p_start(3), p_end(3)], 'k-', 'LineWidth', 1.5);
%    end
%    
%    xlabel('X');
%    ylabel('Y');
%    zlabel('Z');
%    title('Points inside tetrahedron with edges');
%    grid on;
%    axis equal;
%    view(3);
%    hold off;
%end
%
%
%function K = convex_hull_sampled(points, maxPoints, plotHull)
%%CONVEX_HULL_SAMPLED Compute convex hull on a large point set via sampling
%%
%%   points    : M x 3 array of 3D points
%%   maxPoints : maximum number of points to use for hull computation (default 1e6)
%%   plotHull  : logical, true to plot the convex hull (default true)
%%
%%   Returns:
%%       K : indices of the convex hull faces (trisurf format)
%
%if nargin < 2 || isempty(maxPoints)
%    maxPoints = 1e6;
%end
%if nargin < 3
%    plotHull = true;
%end
%
%numPoints = size(points,1);
%
%% Sample points if too many
%if numPoints > maxPoints
%    idx = randperm(numPoints, maxPoints);
%    sampledPoints = points(idx,:);
%else
%    sampledPoints = points;
%end
%
%% Compute convex hull
%K = convhulln(sampledPoints);
%
%% Plot if requested
%if plotHull
%    figure;
%    trisurf(K, sampledPoints(:,1), sampledPoints(:,2), sampledPoints(:,3), ...
%            'FaceColor','cyan', 'FaceAlpha',0.3, 'EdgeColor','k');
%    axis equal;
%    xlabel('X'); ylabel('Y'); zlabel('Z');
%    title(sprintf('Convex Hull of %d Sampled Points', size(sampledPoints,1)));
%    grid on;
%    view(3);
%end
%end
%
%
%
%%%
%function tetraVoxelMap = tetra_to_voxel_mapping(voxelVolume, tetMesh, nodeCoords, tetrahedron_index)
%
%    % Parameters
%
%    voxelVolume = gpuArray(voxelVolume);  % Random voxel values on GPU
%
%    % Generate a simple tetrahedral mesh (example: divide cube into tetrahedra)
%    
%    numTets = size(tetMesh, 1);
%    % Output structure
%    tetraVoxelMap = cell(size(tetMesh, 1), 1);
%
%    % Loop through each tetrahedron
%    
%    Displayer = Progressdisp(numTets, BarOrStr="bar"); % Initialize
%    for i = 1:numTets
%        Displayer.Display_progress(i); % Update progress dynamically
%        nodes = nodeCoords(tetMesh(i, :), :);  % 4x3 matrix
%        [pointsInside, valuesInside] = process_tetrahedron_gpu(nodes, voxelVolume);
%        
%        if isempty(valuesInside)
%            avgVal = NaN;
%        else
%            avgVal = mean(gather(valuesInside));
%        end
%
%        tetraVoxelMap{i} = struct( ...
%            'coords', gather(pointsInside), ...
%            'values', gather(valuesInside), ...
%            'avgValue', avgVal ...
%        );
%
%    end
%    if tetrahedron_index > 0 && tetrahedron_index <= size(tetMesh, 1)
%        plot_tetrahedron_with_internal_points_voxel_values(tetraVoxelMap, tetMesh, nodeCoords, tetrahedron_index);
%    elseif tetrahedron_index == 0
%        disp("Not plotting")
%    else
%        disp("Not a viable tetrahedron index for plotting. Please only use integers between 1 and " + num2str(size(tetMesh, 1)) + ". Use 0 to omit plot.")
%    end
%end
%
%function [pointsInside, values] = process_tetrahedron_gpu(nodes, volume)
%    % Compute bounding box of tetrahedron
%    minCorner = floor(min(nodes));
%    maxCorner = ceil(max(nodes));
%
%    % Make sure bounds are within volume
%    sz = size(volume);
%    minCorner = max(minCorner, [1,1,1]);
%    maxCorner = min(maxCorner, sz);
%
%    % Generate grid of voxel centroids in bounding box
%    [X, Y, Z] = ndgrid(minCorner(1):maxCorner(1), ...
%                       minCorner(2):maxCorner(2), ...
%                       minCorner(3):maxCorner(3));
%    % Shift by -0.5 to get voxel centroid
%    X = gpuArray(X) - 0.5;
%    Y = gpuArray(Y) - 0.5;
%    Z = gpuArray(Z) - 0.5;
%
%    % Flatten
%    points = [X(:), Y(:), Z(:)];
%
%    % Check which points are inside tetrahedron using barycentric method
%    insideMask = point_in_tetrahedron(points, nodes);
%
%    % Get coordinates inside
%    pointsInside = points(insideMask, :);
%
%    % Map to voxel indices (undo -0.5 shift, round to nearest voxel)
%    voxelIdx = round(pointsInside + 0.5);
%
%    % Remove out-of-bound voxels (can happen at borders)
%    validMask = all(voxelIdx >= 1, 2) & ...
%                voxelIdx(:,1) <= sz(1) & ...
%                voxelIdx(:,2) <= sz(2) & ...
%                voxelIdx(:,3) <= sz(3);
%
%    if isempty(validMask)
%
%        centroid = round(mean(nodes, 1));
%        values = volume(centroid(1),centroid(2),centroid(3));
%        pointsInside = [centroid(1) - 0.5,centroid(2) - 0.5,centroid(3) - 0.5];
%    else
%        voxelIdx = voxelIdx(validMask, :);
%        pointsInside = pointsInside(validMask, :);
%    
%        % Convert subscripts to linear indices
%        linIdx = sub2ind(sz, voxelIdx(:,1), voxelIdx(:,2), voxelIdx(:,3));
%    
%        % Get voxel values
%        values = volume(linIdx);
%        
%
%    end
%
%
%end
%
%function inside = point_in_tetrahedron(P, T)
%    % Barycentric coordinate test for points inside a tetrahedron
%    v0 = T(2,:) - T(1,:);
%    v1 = T(3,:) - T(1,:);
%    v2 = T(4,:) - T(1,:);
%    A = [v0', v1', v2'];  % 3x3
%
%    relP = (P - T(1,:))';  % 3xN
%    baryCoords = A \ relP;  % 3xN
%
%    u = baryCoords(1,:)';
%    v = baryCoords(2,:)';
%    w = baryCoords(3,:)';
%    t = 1 - u - v - w;
%
%    % Inside test: all barycentric coords in [0,1]
%    inside = (u >= 0) & (v >= 0) & (w >= 0) & (t >= 0);
%end
%
%function plot_tetrahedron_with_internal_points_voxel_values(tetraVoxelMap, tetMesh, nodeCoords, tet_ind)
%    figure;
%    hold on;
%    axis equal;
%    grid on;
%    view(3);
%    xlabel('X'); ylabel('Y'); zlabel('Z');
%    title('Tetrahedron and Internal Voxel Centroids Colored by Voxel Values');
%
%    % Get first tetrahedron
%    tetNodes = nodeCoords(tetMesh(tet_ind,:), :);
%    internalPoints = tetraVoxelMap{tet_ind}.coords;
%    voxelValues = tetraVoxelMap{tet_ind}.values;
%
%    % Plot the tetrahedron as a semi-transparent 3D surface
%    faces = [
%        1 2 3;
%        1 2 4;
%        1 3 4;
%        2 3 4
%    ];
%
%    patch('Faces', faces, 'Vertices', tetNodes, ...
%        'FaceColor', 'cyan', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
%
%    % Normalize voxel values for coloring (optional)
%    voxelValues = double(voxelValues);
%    cmin = min(voxelValues);
%    cmax = max(voxelValues);
%
%
%    if cmin == cmax
%        % Plot with a fixed color
%        scatter3(internalPoints(:,1), internalPoints(:,2), internalPoints(:,3), ...
%                 50, 'r', 'filled');
%        legend('Tetrahedron', sprintf('Voxel centroids (value = %.3f)', cmin));
%    else
%        % Normalize and plot with color
%        colors = (voxelValues - cmin) / (cmax - cmin);
%        scatter3(internalPoints(:,1), internalPoints(:,2), internalPoints(:,3), ...
%                 50, colors, 'filled');
%        colormap('jet');
%        colorbar;
%        caxis([cmin cmax]);
%        legend('Tetrahedron', 'Voxel centroids (colored by value)');
%    end
%
%end
