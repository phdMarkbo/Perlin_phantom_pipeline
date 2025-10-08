function [all_points, all_vals] = pointsInsideTetrahedraByIndices(nodes,nodes_trans, tets, scale, doPlot, ph)
    % pointsInsideTetrahedraByIndices
    % Inputs:
    %   nodes - Mx3 array of node coordinates
    %   tets  - Nx4 array of tetrahedron node indices (1-based)
    %   scale - scale factor for finer grid resolution
    %   doPlot - logical flag to plot each tetrahedron
    % Output:
    %   all_points - Px3 matrix of all interior points

    %nbr_of_phantom_voxels = size(ph(:), 1);

    all_points = zeros(300000000, 3); %allocating space %%REPLACE WITH SIZE FROM mapped_Struct
    all_vals = zeros(300000000, 1); %allocating space

    numTets = size(tets, 1);
    all_points_idx = 1;

    Displayer1 = Progressdisp(numTets, BarOrStr="bar"); % Initialize progress bar

    for i = 1:numTets

        Displayer1.Display_progress(i); % Show progress bar

        idx = tets(i, :);  % Get node indices for tetrahedron i
        
        % Extract coordinates of the 4 nodes (transpose to get 3x1 vectors)
        P1 = nodes(idx(1), :)';
        P2 = nodes(idx(2), :)';
        P3 = nodes(idx(3), :)';
        P4 = nodes(idx(4), :)';

        P1_trans = nodes_trans(idx(1), :)';
        P2_trans = nodes_trans(idx(2), :)';
        P3_trans = nodes_trans(idx(3), :)';
        P4_trans = nodes_trans(idx(4), :)';
        
        tet_ref = [P1'; P2'; P3'; P4'];
        tet_trans = [P1_trans'; P2_trans'; P3_trans'; P4_trans'];

        % Compute interior points
        [points, vals] = pointsInsideTetrahedronScaled(P1, P2, P3, P4, scale, ph);

        if isempty(points)

            %points = mean(tet_ref);
            continue
        end

        % Compute affine transformation
        points_trans = general_tet_mapping(tet_ref, tet_trans, points')';

        nbr_of_points_in_tetrahedra = size(points, 1);

        start_idx = all_points_idx;
        end_idx = all_points_idx + nbr_of_points_in_tetrahedra - 1;

        all_points(start_idx:end_idx, :) = points_trans;
        all_vals(start_idx:end_idx, :) = vals;

        all_points_idx = all_points_idx + nbr_of_points_in_tetrahedra;

        % Plot if requested
        if doPlot
            plotTetrahedronWithPoints(points, P1, P2, P3, P4);
            title(['Tetrahedron_orig #' num2str(i)]);
        end
        if doPlot
            plotTetrahedronWithPoints(points_trans, P1_trans, P2_trans, P3_trans, P4_trans);
            title(['Tetrahedron #' num2str(i)]);
        end
    end

    all_points = all_points(1:end_idx, :); %crop allocated space at the end of nonzero entries
    all_vals = all_vals(1:end_idx, :);

end
