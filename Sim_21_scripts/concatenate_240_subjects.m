function ConcatenateSims
Sim = 'Sim_21'; 
Nfiles = 100; 
Base = sprintf(['/storage/maullz/Contour_Inference_2018/', Sim, '_results/240_subjects']);

cd(Base)

a = dir([Sim,'_*.mat']);

x = load(a(1).name);
filename = sprintf([Sim,'_240_subjects.mat']);
save(filename,'-struct','x');

for j=2:Nfiles
	x = load(filename);
	y = load(a(j).name);

	vrs = fieldnames(x);

	x.(vrs{2}) = x.(vrs{2}) + y.(vrs{2}); % Adding number of realizations

	for k=9:23
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
	end

    	x.(vrs{24}) = [x.(vrs{24}), y.(vrs{24})]; % Concatenates the SupGstore along the 2nd dimension instead of 1st dimension 

	for k=25:36
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
    	end

save(filename,'-struct','x');
end

% Recalculating the PercentageSuccessVector now that the SubsetSuccessVector has results for all realizations
x.percentage_success_vector_raw_80 = mean(x.subset_success_vector_raw_80, 1); 
x.percentage_success_vector_raw_90 = mean(x.subset_success_vector_raw_90, 1);
x.percentage_success_vector_raw_95 = mean(x.subset_success_vector_raw_95, 1);
x.percentage_success_vector_raw_80_alternate = mean(x.subset_success_vector_raw_80_alternate, 1); 
x.percentage_success_vector_raw_90_alternate = mean(x.subset_success_vector_raw_90_alternate, 1);
x.percentage_success_vector_raw_95_alternate = mean(x.subset_success_vector_raw_95_alternate, 1);

save(filename, '-struct','x');

end
