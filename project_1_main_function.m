clear;
reset(gpuDevice());

addpath("/home/user/code/Perlin_phantom_pipeline/function_files")
addpath("/home/user/code/Perlin_phantom_pipeline/load_files/")
addpath("/home/user/code/Perlin_phantom_pipeline/Progressdisp.m/")

phantom_struct = load("PerlinPhantom01.mat");

elems = load("element.dat");
elems = elems(:, 2:5);
scale_factor = 10000; %from meters to units of 100um
nodes_orig_mm = load("node.dat"); %in mm
nodes_orig_m = nodes_orig_mm(:, 2:4) * 0.001; %in m
nodes_orig_in_100um_norm = nodes_orig_m * scale_factor + 1; %in units of 100 um, shifted to non-zero values

clear nodes_orig_mm;
clear nodes_orig_m;
%nodes_comp = readtable("cupA_01_all_materials_compressed_nodes.txt");
%nodes_comp = table2array(nodes_comp);

%nodes_orig_m = nodes_orig * 0.001;
%nodes_comp = nodes_comp(1:29479, 2:4);
%elems = elems(:, 2:5);

%nodes_orig_scaled = nodes_orig * 10000 + 1;


tet_ind_to_be_plotted = 1;
nbr_of_phantom_nodes = 29479;
last_compression_step_nbr = 9;
materia_ratio_diff_val = 0.05;
E_mod_mat1 = 10000;
E_mod_mat2 = 25000;

replacement_list_FROM_VCT_format = [8, 7; 9, 8; 11, 9; 12, 10; 13,11];
replacement_list_TO_VCT_format = [11, 13; 10, 12; 9, 11; 8,9 ; 7, 8];
finalVolume = phantom_struct.finalVolume;
finalVolume = replaceValues(finalVolume, replacement_list_FROM_VCT_format);
phantom_size = size(finalVolume);
clear phantom_struct;
%%
%%%SCRIPT1

mapping_struct = tetra_to_voxel_mapping(finalVolume, elems, nodes_orig_in_100um_norm,tet_ind_to_be_plotted);
clear finalVolume;
%% FEBio
    
fid = fopen('test.txt', 'w');  % opens and clears the file
if fid == -1
    error('Cannot open model.feb for clearing.');
end
fclose(fid);

test_list_of_ratios = mapValueWithDiff(mapping_struct, materia_ratio_diff_val);
E_module_values = weightedAverage(E_mod_mat1, E_mod_mat2, test_list_of_ratios);
appendValuesToFile_chunked("test.txt", E_module_values);

fid = fopen('test.feb', 'w');  % opens and clears the file
if fid == -1
    error('Cannot open model.feb for clearing.');
end
fclose(fid);

AppendFileToFeb("test.feb", "/home/user/code/Perlin_phantom_pipeline/cupA_template_top.txt");
AppendFileToFeb("test.feb", "test.txt");
AppendFileToFeb("test.feb", "/home/user/code/Perlin_phantom_pipeline/cupA_template_bottom.txt");

%%PERFORM FEBio COMPRESSION TO GET LOG FILE
comp_nodes_array = extractNodalDataBlock("test.log", last_compression_step_nbr, "compressed_nodes2.txt");
nodes_comp = comp_nodes_array(1:nbr_of_phantom_nodes, 2:4); %in mm
clear comp_nodes_array;
clear E_module_values;
clear test_list_of_ratios;
%%
%SCRIPT2

nodes_orig_mm = load("node.dat"); %in mm
nodes_orig_m = nodes_orig_mm(:, 2:4) * 0.001; %in m
clear nodes_orig_mm;
nbr_of_coords = countAllValues(mapping_struct, "values");
coord_values = concatenateFieldValues(mapping_struct, "values");
clear mapping_struct;
clear nodes_orig_in_100um_norm;
reset(gpuDevice())
affine_coords = pointsInsideTetrahedraByIndices(nodes_orig_m, nodes_comp, elems, 10000, false, phantom_size, nbr_of_coords * 2); %allocating more than necessary for safetsy sake
clear elems;
clear nodes_orig_m;
clear nodes_comp;
%%
%SCRIPT3

compressed_phantom = phantom_from_transformed_points(affine_coords, coord_values);
compressed_phantom_vct = replaceValues(compressed_phantom, replacement_list_TO_VCT_format); %Back to openVCT format
