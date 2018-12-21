Basedir = '/storage/maullz/Contour_Inference_2018';

outdir  = fullfile(Basedir,'Figures','Supplemental_Figures');

if ~isdir(outdir)
    mkdir(outdir)
end

sim_results_mat = {fullfile(Basedir,'Sim_55_results','60_subjects','Sim_55_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_55_results','120_subjects','Sim_55_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_55_results','240_subjects','Sim_55_240_subjects.mat'); ...
                   fullfile(Basedir,'Sim_38_results','60_subjects','Sim_38_60_subjects.mat'), ...
                   fullfile(Basedir,'Sim_38_results','120_subjects','Sim_38_120_subjects.mat'), ...
                   fullfile(Basedir,'Sim_38_results','240_subjects','Sim_38_240_subjects.mat')};
               
sim_results_mat_dim = size(sim_results_mat);

nominal_levels      = [80, 90, 95];
nominal_prct_levels = nominal_levels/100;
nominal_levels_dim  = size(nominal_levels);

result = load(sim_results_mat{1});
nRlz          = result.nRlz;

std_error_vector = zeros(nominal_levels_dim);
sim_covering_levels_and_subs = zeros([5 sim_results_mat_dim(2) sim_results_mat_dim(1) nominal_levels_dim(2)]);

for i = 1:nominal_levels_dim(2)
    std_error = sqrt((nominal_prct_levels(i)*(1 - nominal_prct_levels(i)))/nRlz);
    std_error_vector(i) = std_error; 
 
    for j = 1:sim_results_mat_dim(1)
        for k = 1:sim_results_mat_dim(2)
            temp = load(char(sim_results_mat(j,k)));
            raw_result                          = sprintf('percentage_success_vector_raw_%d_alternate', nominal_levels(i));
            raw_result_old                      = sprintf('percentage_success_vector_raw_%d', nominal_levels(i));
            observed_result                     = sprintf('percentage_success_vector_observed_%d_alternate', nominal_levels(i));
            observed_result_old                  = sprintf('percentage_success_vector_observed_%d', nominal_levels(i));
            nSubj                               = sprintf('nSubj');
            sim_covering_levels_and_subs(1,k,j,i) = temp.(raw_result);
            sim_covering_levels_and_subs(2,k,j,i) = temp.(observed_result);
            sim_covering_levels_and_subs(3,k,j,i) = temp.(raw_result_old);
            sim_covering_levels_and_subs(4,k,j,i) = temp.(observed_result_old);
            sim_covering_levels_and_subs(5,k,j,i) = temp.(nSubj);
        end
    end
end

results_params = struct('nom_80_results', struct('nominal_level', nominal_levels(1), ...
                                                 'nominal_prct_level', nominal_prct_levels(1), ...
                                                 'std_error', std_error_vector(1), ...
                                                 'sig_2_t_bootstrap_results', sim_covering_levels_and_subs(:,:,1,1),  ...
                                                 'sig_2_old_bootstrap_results', sim_covering_levels_and_subs(:,:,2,1)),  ...
                        'nom_90_results', struct('nominal_level', nominal_levels(2), ...
                                                 'nominal_prct_level', nominal_prct_levels(2), ...
                                                 'std_error', std_error_vector(2), ...
                                                 'sig_2_t_bootstrap_results', sim_covering_levels_and_subs(:,:,1,2), ...
                                                 'sig_2_old_bootstrap_results', sim_covering_levels_and_subs(:,:,2,2)), ...
                        'nom_95_results', struct('nominal_level', nominal_levels(3), ...
                                                 'nominal_prct_level', nominal_prct_levels(3), ...
                                                 'std_error', std_error_vector(3), ...
                                                 'sig_2_t_boostrap_results', sim_covering_levels_and_subs(:,:,1,3), ...
                                                 'sig_2_old_bootstrap_results', sim_covering_levels_and_subs(:,:,2,3)) ...                                                   
                        );
                    
filename = fullfile(outdir, 'SupplementalResults.mat');
save(filename,'results_params');

