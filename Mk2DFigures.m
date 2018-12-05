results_mat_dir = '/Users/maullz/Desktop/Contour_Inference_2018/Figures/2D_Figures';

results_mat_file = fullfile(results_mat_dir,'2DResults.mat');

load(results_mat_file)

nominal_vec = ["nom_80_results","nom_90_results","nom_95_results"];
signal_vec  = ["sig_1_std_1_results","sig_1_std_2_results","sig_2_std_1_results","sig_2_std_2_results"];
color_vec   = 'rbrb';

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 1:2
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(5,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.5 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('emp. Covering Rate');
    
    titlename = sprintf('%d%% Nominal Coverage Results', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Sig. 1, Std Dev. 1 (est. boundary)', ...
                     'Sig. 1, Std Dev. 1 (true boundary)', ...
                     'Sig. 1, Std Dev. 2 (est. boundary)', ...
                     'Sig. 1, Std Dev. 2 (true boundary)', ...
                     'Nominal Coverage Level', ...
                     '1.96 * Std Error');
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.25, pos_lgd(3), pos_lgd(4) - 0.2]);
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'Sig_1_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    
    for j = 3:4
        signal = signal_vec(j);
        [true_bdry_cov] = results_params.(results).(signal)(1,:);
        [est_bdry_cov]  = results_params.(results).(signal)(2,:);
        [subs]          = results_params.(results).(signal)(5,:);
        plot(subs,est_bdry_cov,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_cov,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.5 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('emp. Covering Rate');
    
    titlename = sprintf('%d%% Nominal Coverage Results', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Sig. 2, Std Dev. 1 (est. boundary)', ...
                     'Sig. 2, Std Dev. 1 (true boundary)', ...
                     'Sig. 2, Std Dev. 2 (est. boundary)', ...
                     'Sig. 2, Std Dev. 2 (true boundary)', ...
                     'Nominal Coverage Level', ...
                     '1.96 * Std Error');
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.25, pos_lgd(3), pos_lgd(4) - 0.2]);
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'Sig_2_coverage_results.pdf'))