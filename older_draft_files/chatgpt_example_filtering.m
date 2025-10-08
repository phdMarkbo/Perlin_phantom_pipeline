% %% Parameters
% cubeSize = 30;        % number of points per dimension
% R = 5;                % sphere radius
% 
% % Generate cube coordinates
% [x, y, z] = meshgrid(linspace(0,1,cubeSize), linspace(0,1,cubeSize), linspace(0,1,cubeSize));
% 
% %% Step 1: Generate 3D Gaussian noise
% noiseCube = randn(cubeSize, cubeSize, cubeSize);
% 
% %% Step 2: Apply Sobel filter in x-direction
% [gradX, ~, ~] = imgradientxyz(noiseCube, 'sobel');
% gradX = gradX * 5; % amplify streaks
% 
% %% Step 3: Map cube to spherical coordinates
% % Normalize cube coordinates (0-1)
% x_norm = x; y_norm = y; z_norm = z;
% 
% % Map to spherical coordinates
% r = R * x_norm;               % radius scales with x
% theta = 2*pi*y_norm;          % azimuthal angle 0 to 2pi
% phi = pi*z_norm;              % polar angle 0 to pi
% 
% % Convert to Cartesian coordinates (sphere)
% X = r .* sin(phi) .* cos(theta);
% Y = r .* sin(phi) .* sin(theta);
% Z = r .* cos(phi);
% 
% %% Step 4: Map noise values to new spherical coordinates
% % For simplicity, we use the same noise values at each grid point
% % The color represents noise amplitude
% colorOriginal = noiseCube(:);
% colorSobel = gradX(:);
% 
% %% Step 5: Visualization
% figure('Position',[50 50 1800 450]);
% 
% % Original cube noise (slice)
% subplot(1,4,1)
% sliceIdx = round(cubeSize/2);
% imagesc(noiseCube(:,:,sliceIdx))
% axis equal tight; colormap gray;
% title('Original Noise Cube (slice)')
% 
% % Sobel filtered cube (slice)
% subplot(1,4,2)
% imagesc(gradX(:,:,sliceIdx))
% axis equal tight; colormap gray;
% title('Sobel X-direction (slice)')
% 
% % Noise mapped to sphere
% subplot(1,4,3)
% scatter3(X(:), Y(:), Z(:), 20, colorOriginal, 'filled')
% axis equal
% title('Sphere (Original Noise)')
% xlabel('X'); ylabel('Y'); zlabel('Z'); colormap jet; colorbar
% 
% % Sobel mapped to sphere
% subplot(1,4,4)
% scatter3(X(:), Y(:), Z(:), 20, colorSobel, 'filled')
% axis equal
% title('Sphere (Sobel Applied)')
% xlabel('X'); ylabel('Y'); zlabel('Z'); colormap jet; colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%EXAMPLE 2
% %% Parameters
% cubeSize = 50;        % number of points per dimension
% R = 5;                % sphere radius
% 
% % Generate cube coordinates
% [x, y, z] = meshgrid(linspace(0,1,cubeSize), linspace(0,1,cubeSize), linspace(0,1,cubeSize));
% 
% %% Step 1: Gaussian noise
% noiseCube = randn(cubeSize, cubeSize, cubeSize);
% 
% %% Step 2: Sobel in x-direction
% [gradX, ~, ~] = imgradientxyz(noiseCube, 'sobel');
% gradX = gradX * 5; % amplify streaks
% 
% %% Step 3: Map cube to spherical coordinates
% x_norm = x; y_norm = y; z_norm = z;
% 
% r = R * x_norm;               % radius scales with x
% theta = 2*pi*y_norm;          % azimuth
% phi = pi*z_norm;              % polar
% 
% X = r .* sin(phi) .* cos(theta);
% Y = r .* sin(phi) .* sin(theta);
% Z = r .* cos(phi);
% 
% %% Step 4: Full 3D scatter (semi-transparent)
% figure('Position',[100 100 1200 500]);
% 
% % Original noise mapped to sphere
% subplot(1,2,1)
% scatter3(X(:), Y(:), Z(:), 15, noiseCube(:), 'filled')
% axis equal
% colormap jet; colorbar
% title('Sphere (Original Noise)')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(3)
% alpha(0.7)
% 
% % Sobel mapped to sphere
% subplot(1,2,2)
% scatter3(X(:), Y(:), Z(:), 15, gradX(:), 'filled')
% axis equal
% colormap jet; colorbar
% title('Sphere (Sobel Applied)')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(3)
% alpha(0.7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%EXAMPLE 3
% %% Parameters
% cubeSize = 60;        % finer grid for smoother streaks
% R = 5;                % sphere radius
% 
% % Generate cube coordinates
% [x, y, z] = meshgrid(linspace(0,1,cubeSize), linspace(0,1,cubeSize), linspace(0,1,cubeSize));
% 
% %% Step 1: Gaussian noise
% noiseCube = randn(cubeSize, cubeSize, cubeSize);
% 
% %% Step 2: Sobel in x-direction
% [gradX, ~, ~] = imgradientxyz(noiseCube, 'sobel');
% gradX = gradX * 5; % amplify streaks
% 
% %% Step 3: Map cube to spherical coordinates
% x_norm = x; y_norm = y; z_norm = z;
% 
% r = R * x_norm;               % radius scales with x
% theta = 2*pi*y_norm;          % azimuth
% phi = pi*z_norm;              % polar
% 
% X = r .* sin(phi) .* cos(theta);
% Y = r .* sin(phi) .* sin(theta);
% Z = r .* cos(phi);
% 
% %% Step 4: Optional Gaussian smoothing (makes streaks more continuous)
% gradSmooth = imgaussfilt3(gradX, 1.2);   % tweak sigma for smoothness
% 
% %% Step 5: 3D scatter plot (smooth, continuous)
% figure('Position',[100 100 1200 500]);
% 
% % Original noise
% subplot(1,2,1)
% scatter3(X(:), Y(:), Z(:), 5, noiseCube(:), 'filled')
% axis equal tight
% colormap jet; colorbar
% title('Sphere (Original Noise)')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(35,30)
% grid on
% 
% % Smoothed Sobel streaks
% subplot(1,2,2)
% scatter3(X(:), Y(:), Z(:), 5, gradSmooth(:), 'filled')
% axis equal tight
% colormap jet; colorbar
% title('Sphere (Smoothed Sobel)')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(35,30)
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%EXAMPLE 4
% %% Parameters
% cubeSize = 60;    % number of points per dimension
% R = 5;            % half-sphere radius
% 
% % Generate cube coordinates
% [x, y, z] = meshgrid(linspace(0,1,cubeSize), linspace(0,1,cubeSize), linspace(0,1,cubeSize));
% 
% %% Step 1: Gaussian noise
% noiseCube = randn(cubeSize, cubeSize, cubeSize);
% 
% %% Step 2: Sobel in y-direction
% [~, gradY, ~] = imgradientxyz(noiseCube, 'sobel'); 
% gradY = gradY * 5;  % amplify streaks
% 
% %% Step 3: Map cube to half-sphere
% x_norm = x; y_norm = y; z_norm = z;
% 
% r = R * x_norm;                 % radius scales with x
% theta = pi * y_norm;            % azimuth (0 to pi) for half-sphere
% phi = (pi/2) * z_norm;          % polar angle (0 to pi/2)
% 
% % Convert to Cartesian coordinates
% X = r .* sin(phi) .* cos(theta);
% Y = r .* sin(phi) .* sin(theta);
% Z = r .* cos(phi);
% 
% %% Step 4: Optional smoothing
% gradSmooth = imgaussfilt3(gradY, 1.2);
% 
% %% Step 5: Visualize full half-sphere
% figure('Position',[100 100 1800 500]);
% 
% subplot(1,3,1)
% scatter3(X(:), Y(:), Z(:), 5, noiseCube(:), 'filled')
% axis equal tight
% colormap jet; colorbar
% title('Half-Sphere (Original Noise)')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(35,30)
% grid on
% 
% subplot(1,3,2)
% scatter3(X(:), Y(:), Z(:), 5, gradSmooth(:), 'filled')
% axis equal tight
% colormap jet; colorbar
% title('Half-Sphere (Sobel Y-direction)')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(35,30)
% grid on
% 
% %% Step 6: Longitudinal slice in X-Z plane
% % Select a thin azimuthal slice (theta) along Y axis
% thetaSliceCenter = pi/2;       % center of slice in theta
% sliceThickness = pi/20;        % small angular thickness
% 
% sliceMask = (theta >= thetaSliceCenter - sliceThickness/2) & ...
%             (theta <= thetaSliceCenter + sliceThickness/2);
% 
% Xslice = X(sliceMask);
% Zslice = Z(sliceMask);
% sobelSlice = gradSmooth(sliceMask);
% 
% subplot(1,3,3)
% scatter(Xslice, Zslice, 15, sobelSlice, 'filled')
% axis equal tight
% colormap jet; colorbar
% title('Half-Sphere Longitudinal Slice (X-Z, Sobel)')
% xlabel('X'); ylabel('Z')
% grid on

% %%EXAMPLE 5
%% Parameters
%% Parameters
cubeSize = 60;    % number of points per dimension
R = 5;            % half-sphere radius

% Generate cube coordinates
[x, y, z] = meshgrid(linspace(0,1,cubeSize), linspace(0,1,cubeSize), linspace(0,1,cubeSize));

%% Step 1: Gaussian noise
noiseCube = randn(cubeSize, cubeSize, cubeSize);

%% Step 2: Sobel in y-direction
[~, gradY, ~] = imgradientxyz(noiseCube, 'sobel'); 
gradY = gradY * 5;  % amplify streaks

%% Step 3: Map cube to half-sphere
x_norm = x; y_norm = y; z_norm = z;

r = R * x_norm;                 % radius scales with x
theta = pi * y_norm;            % azimuth (0 to pi) for half-sphere
phi = (pi/2) * z_norm;          % polar angle (0 to pi/2)

% Convert to Cartesian coordinates
X = r .* sin(phi) .* cos(theta);
Y = r .* sin(phi) .* sin(theta);
Z = r .* cos(phi);

%% Step 4: Optional smoothing
gradSmooth = imgaussfilt3(gradY, 1.2);

%% Step 5: Visualize full half-sphere (first figure)
figure('Position',[100 100 1200 500]);

subplot(1,2,1)
scatter3(X(:), Y(:), Z(:), 5, noiseCube(:), 'filled')
axis equal tight
colormap jet; colorbar
title('Half-Sphere (Original Noise)')
xlabel('X'); ylabel('Y'); zlabel('Z')
view(35,30)
grid on

subplot(1,2,2)
scatter3(X(:), Y(:), Z(:), 5, gradSmooth(:), 'filled')
axis equal tight
colormap jet; colorbar
title('Half-Sphere (Sobel Y-direction)')
xlabel('X'); ylabel('Y'); zlabel('Z')
view(35,30)
grid on

%% Step 6: Longitudinal slice in X-Z plane (separate figure)
thetaSliceCenter = pi/2;       % center of slice in theta
sliceThickness = pi/20;        % small angular thickness

sliceMask = (theta >= thetaSliceCenter - sliceThickness/2) & ...
            (theta <= thetaSliceCenter + sliceThickness/2);

Xslice = X(sliceMask);
Zslice = Z(sliceMask);
sobelSlice = gradSmooth(sliceMask);

figure('Position',[200 200 600 600])
scatter(Xslice, Zslice, 15, sobelSlice, 'filled')
axis equal tight
colormap jet; colorbar
title('Half-Sphere Longitudinal Slice (X-Z, Sobel)')
xlabel('X'); ylabel('Z')
grid on

