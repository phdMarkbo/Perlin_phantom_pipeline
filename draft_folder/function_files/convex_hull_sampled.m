function K = convex_hull_sampled(points, maxPoints, plotHull)
%CONVEX_HULL_SAMPLED Compute convex hull on a large point set via sampling
%
%   points    : M x 3 array of 3D points
%   maxPoints : maximum number of points to use for hull computation (default 1e6)
%   plotHull  : logical, true to plot the convex hull (default true)
%
%   Returns:
%       K : indices of the convex hull faces (trisurf format)

if nargin < 2 || isempty(maxPoints)
    maxPoints = 1e6;
end
if nargin < 3
    plotHull = true;
end

numPoints = size(points,1);

% Sample points if too many
if numPoints > maxPoints
    idx = randperm(numPoints, maxPoints);
    sampledPoints = points(idx,:);
else
    sampledPoints = points;
end

% Compute convex hull
K = convhulln(sampledPoints);

% Plot if requested
if plotHull
    figure;
    trisurf(K, sampledPoints(:,1), sampledPoints(:,2), sampledPoints(:,3), ...
            'FaceColor','cyan', 'FaceAlpha',0.3, 'EdgeColor','k');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Convex Hull of %d Sampled Points', size(sampledPoints,1)));
    grid on;
    view(3);
end
end

