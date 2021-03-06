function ConcatenateSims
Sim = 'Sim_20'; 
Nfiles = 100; 
Base = sprintf(['/storage/maullz/Contour_Inference_2018/', Sim, '_results/120_subjects']);

cd(Base)

a = dir([Sim,'_*.mat']);

x = load(a(1).name);
filename = sprintf([Sim,'_120_subjects.mat']);
save(filename,'-struct','x');

for j=2:Nfiles
	x = load(filename);
	y = load(a(j).name);

	vrs = fieldnames(x);

	x.(vrs{2}) = x.(vrs{2}) + y.(vrs{2}); % Adding number of realizations

	for k=9:56
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
	end

	for k=57:60
    		x.(vrs{k}) = [x.(vrs{k}), y.(vrs{k})]; % Concatenates the SupGstore along the 2nd dimension instead of 1st dimension 
	end

	for k=61:108
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
    end

save(filename,'-struct','x');
end

% Recalculating the PercentageSuccessVector now that the SubsetSuccessVector has results for all realizations
x.percentage_success_vector_raw_80 = mean(x.subset_success_vector_raw_80, 1); 
x.percentage_success_vector_raw_90 = mean(x.subset_success_vector_raw_90, 1);
x.percentage_success_vector_raw_95 = mean(x.subset_success_vector_raw_95, 1);
x.percentage_success_vector_raw_80_weighted = mean(x.subset_success_vector_raw_80_weighted, 1); 
x.percentage_success_vector_raw_90_weighted = mean(x.subset_success_vector_raw_90_weighted, 1);
x.percentage_success_vector_raw_95_weighted = mean(x.subset_success_vector_raw_95_weighted, 1);
x.percentage_success_vector_raw_80_linear = mean(x.subset_success_vector_raw_80_linear, 1); 
x.percentage_success_vector_raw_90_linear = mean(x.subset_success_vector_raw_90_linear, 1);
x.percentage_success_vector_raw_95_linear = mean(x.subset_success_vector_raw_95_linear, 1);
x.percentage_success_vector_raw_80_observed_weighted = mean(x.subset_success_vector_raw_80_observed_weighted, 1); 
x.percentage_success_vector_raw_90_observed_weighted = mean(x.subset_success_vector_raw_90_observed_weighted, 1);
x.percentage_success_vector_raw_95_observed_weighted = mean(x.subset_success_vector_raw_95_observed_weighted, 1);

save(filename, '-struct','x');

end
