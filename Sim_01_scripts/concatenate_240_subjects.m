function ConcatenateSims
Sim = 'Sim_01'; 
Nfiles = 30; 
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

	for k=9:12
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
	end
    
	for k = 13:14
		x.(vrs{k}) = [x.(vrs{k}), y.(vrs{k})]; % Concatenates the SupGstore along the 2nd dimension instead of 1st dimension
	end 

	for k=15:24
		x.(vrs{k}) = [x.(vrs{k}); y.(vrs{k})]; % Concatenating the other variables of interest
    end

save(filename,'-struct','x');
end

x.percentage_success_vector_raw = mean(x.subset_success_vector_raw, 1); % Recalculating the PercentageSuccessVector now that the SubsetSuccessVector has results for all realizations
x.percentage_success_vector_cohen = mean(x.subset_success_vector_cohen, 1);
save(filename, '-struct','x');

end
