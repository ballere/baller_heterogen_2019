%Load the hydra output .mat file
load('/cbica/projects/pncHeterogeneity/ballerDepHeterogen/results/new_hydra_results_path/HYDRA_results.mat');

%Save the ARI values as a .csv file:

dlmwrite('/cbica/projects/pncHeterogeneity/ballerDepHeterogen/results/new_hydra_results_prefix_ARI.csv',ARI);

%Combine the IDs and clusters assignments (CIDX) into a single matrix

combined = [ID CIDX];

%Define header names

headers = {'bblid','Hydra_k1','Hydra_k2','Hydra_k3','Hydra_k4','Hydra_k5','Hydra_k6','Hydra_k7','Hydra_k8','Hydra_k9','Hydra_k10'};

%Save the combined matrix with headers into a .csv file using the function "csvwrite_with_headers"

csvwrite_with_headers('/cbica/projects/pncHeterogeneity/ballerDepHeterogen/results/new_hydra_results_prefix_HydraSubtypes.csv',combined,headers);

exit;

