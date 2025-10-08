%%Testing mapping list of tetrahderon encapsulated points to another
%%tetrahedron
%phantom_struct = load("PerlinPhantom01.mat");
% elems = load("element.dat");
% nodes_orig = load("node.dat");
% nodes_comp = readtable("cupA_01_all_materials_compressed_nodes.txt");
% 
% nodes_comp = table2array(nodes_comp);
%finalVolume = phantom_struct.finalVolume;
% 
% nodes_orig = nodes_orig(:, 2:4) * 0.001;
% nodes_comp = nodes_comp(:, 2:4);
% elems = elems(:, 2:5);
affine_coords = load("inside_ps_coords.mat");
affine_vals = load("inside_ps_vals_orig.mat");

affine_coords = affine_coords.inside_ps;
affine_vals = affine_vals.inside_val_orig;

tic;
affine_coords_nz = (affine_coords - min(affine_coords))*10000;

integer_coords_nz = uint16(round(affine_coords_nz)) + 1;

max_x = max(integer_coords_nz(:, 1));
max_y = max(integer_coords_nz(:, 2));
max_z = max(integer_coords_nz(:, 3));
bounding_box = false(max_y, max_x, max_z);

for i = 1:length(integer_coords_nz)

    cur_X = integer_coords_nz(i, 1);
    cur_Y = integer_coords_nz(i, 2);
    cur_Z = integer_coords_nz(i, 3);
 
    bounding_box(cur_Y, cur_X, cur_Z) = true;
end

strel = ones(50, 50, 50);

filled_bounding_vol = imclose(bounding_box,strel);
filled_bounding_vol_uint16 = uint16(filled_bounding_vol);
filled_bounding_vol_uint16(filled_bounding_vol_uint16 == 0) = intmax("uint16");
filled_bounding_vol_uint16(filled_bounding_vol_uint16 == 1) = 0;
for i = 1:size(integer_coords_nz,1)
    
    cur_X = integer_coords_nz(i, 1);
    cur_Y = integer_coords_nz(i, 2);
    cur_Z = integer_coords_nz(i, 3);
    filled_bounding_vol_uint16(cur_Y, cur_X, cur_Z) = affine_vals(i);

end

nonzeromask = filled_bounding_vol_uint16 ~= 0; %All breast voxels with assigned values
[~, idxClosest] = bwdist(nonzeromask);
filled_phantom = filled_bounding_vol_uint16(idxClosest);
filled_phantom(filled_phantom==intmax("uint16")) = 0;
filled_phantom = filled_phantom .* uint16(filled_bounding_vol); %%mask it
%masked_Perlin_phantom = finalVolume(1:750, 1:1700, 1:374) .* uint16(filled_bounding_vol); 

elapsedTime = toc;

fprintf('Time spent calculating: %.2f seconds (%.2f minutes)\n', ...
        elapsedTime, elapsedTime/60);
