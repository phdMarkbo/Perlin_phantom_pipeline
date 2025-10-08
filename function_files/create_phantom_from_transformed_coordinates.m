function filled_phantom = create_phantom_from_transformed_coordinates(affine_coords, inside_vals)
% PHANTOM_FROM_TRANSFORMED_POINTS THIS REPLACES THE OLD FILLED PHANTOM FUNCTION
% Builds a filled 3D phantom volume from transformed affine coordinates and values.
%
% Inputs:
%   affine_coords - N×3 array of XYZ coordinates (in affine space)
%   inside_vals   - N×1 vector of values to assign to those coordinates
%
% Output:
%   filled_phantom - 3D uint16 volume with filled interior regions

    % --- Normalize coordinates and scale ---
    affine_coords_nz = affine_coords - min(affine_coords, [], 1);
    affine_coords_nz = affine_coords_nz * 10000;

    % --- Convert to integer voxel indices ---
    integer_coords_nz = uint16(round(affine_coords_nz)) + 1;

    % --- Determine bounding box dimensions ---
    max_coords = max(integer_coords_nz, [], 1);
    bounding_box = false(max_coords(2), max_coords(1), max_coords(3)); % (Y,X,Z)

    % --- Populate voxel mask efficiently ---
    lin_idx = sub2ind(size(bounding_box), ...
                      integer_coords_nz(:,2), ...
                      integer_coords_nz(:,1), ...
                      integer_coords_nz(:,3));
    bounding_box(lin_idx) = true;

    % --- Morphological closing to fill internal gaps ---
    se = ones(50, 50, 50);
    filled_bounding_vol = imclose(bounding_box, se);

    % --- Initialize volume with uint16 max value (for later masking) ---
    filled_bounding_vol_uint16 = uint16(filled_bounding_vol);
    filled_bounding_vol_uint16(~filled_bounding_vol) = intmax("uint16");
    filled_bounding_vol_uint16(filled_bounding_vol)  = 0;

    % --- Assign inside values ---
    filled_bounding_vol_uint16(lin_idx) = inside_vals;

    % --- Distance transform filling ---
    nonzeromask = filled_bounding_vol_uint16 ~= 0;
    [~, idxClosest] = bwdist(nonzeromask);
    filled_phantom = filled_bounding_vol_uint16(idxClosest);

    % --- Cleanup and mask ---
    filled_phantom(filled_phantom == intmax("uint16")) = 0;
    filled_phantom = filled_phantom .* uint16(filled_bounding_vol);
end

