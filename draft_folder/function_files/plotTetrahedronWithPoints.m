function plotTetrahedronWithPoints(points, P1, P2, P3, P4)
    % plotTetrahedronWithPoints
    % Inputs:
    %   points - Nx3 matrix of points inside tetrahedron
    %   P1, P2, P3, P4 - 3x1 vectors of tetrahedron vertices
    
    vertices = [P1, P2, P3, P4];
    
    edges = [
        1 2;
        1 3;
        1 4;
        2 3;
        2 4;
        3 4];
    
    figure;
    hold on;
    scatter3(points(:,1), points(:,2), points(:,3), 36, 'filled');
    
    for i = 1:size(edges,1)
        p_start = vertices(:, edges(i,1));
        p_end = vertices(:, edges(i,2));
        plot3([p_start(1), p_end(1)], ...
              [p_start(2), p_end(2)], ...
              [p_start(3), p_end(3)], 'k-', 'LineWidth', 1.5);
    end
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Points inside tetrahedron with edges');
    grid on;
    axis equal;
    view(3);
    hold off;
end

