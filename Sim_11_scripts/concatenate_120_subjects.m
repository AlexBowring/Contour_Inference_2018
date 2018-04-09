function ConcatenateSims
Sim = 'Sim_11'; 
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

	for k=9:54
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
	end

	for k=55:57
    		x.(vrs{k}) = [x.(vrs{k}), y.(vrs{k})]; % Concatenates the SupGstore along the 2nd dimension instead of 1st dimension 
	end

	for k=58:105
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
    	end

	x.(vrs{107}) = [x.(vrs{k}); y.(vrs{k})];

save(filename,'-struct','x');
end

% Recalculating the PercentageSuccessVector now that the SubsetSuccessVector has results for all realizations
x.percentage_success_vector_raw_80_gaussian = mean(x.subset_success_vector_raw_80_gaussian, 1); 
x.percentage_success_vector_raw_90_gaussian = mean(x.subset_success_vector_raw_90_gaussian, 1);
x.percentage_success_vector_raw_95_gaussian = mean(x.subset_success_vector_raw_95_gaussian, 1);
x.percentage_success_vector_raw_80_ero_dil_gaussian = mean(x.subset_success_vector_raw_80_ero_dil_gaussian, 1); 
x.percentage_success_vector_raw_90_ero_dil_gaussian = mean(x.subset_success_vector_raw_90_ero_dil_gaussian, 1);
x.percentage_success_vector_raw_95_ero_dil_gaussian = mean(x.subset_success_vector_raw_95_ero_dil_gaussian, 1);
x.percentage_success_vector_raw_80_signflips = mean(x.subset_success_vector_raw_80_signflips, 1); 
x.percentage_success_vector_raw_90_signflips = mean(x.subset_success_vector_raw_90_signflips, 1);
x.percentage_success_vector_raw_95_signflips = mean(x.subset_success_vector_raw_95_signflips, 1);
x.percentage_success_vector_raw_80_ero_dil_signflips = mean(x.subset_success_vector_raw_80_ero_dil_signflips, 1); 
x.percentage_success_vector_raw_90_ero_dil_signflips = mean(x.subset_success_vector_raw_90_ero_dil_signflips, 1);
x.percentage_success_vector_raw_95_ero_dil_signflips = mean(x.subset_success_vector_raw_95_ero_dil_signflips, 1);

save(filename, '-struct','x');

end
