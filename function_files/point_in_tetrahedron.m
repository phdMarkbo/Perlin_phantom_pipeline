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

