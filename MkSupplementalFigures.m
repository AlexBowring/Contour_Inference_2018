results_mat_dir = '/Users/maullz/Desktop/Contour_Inference_2018/Figures/Supplemental_Figures';

results_mat_file = fullfile(results_mat_dir,'SupplementalResults.mat');

load(results_mat_file)

nominal_vec = ["nom_80_results","nom_90_results","nom_95_results"];
signal_vec  = ["ThreeD_sig_2_t_bootstrap_results","ThreeD_sig_2_old_bootstrap_results","TwoD_sig_2_t_bootstrap_results","TwoD_sig_2_old_bootstrap_results"];
color_vec   = 'rbmg';

% Creating coverage plots for first signal model
figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    for j = 1:2
        signal = signal_vec(j);
        [new_method_bdry_cov] = results_params.(results).(signal)(2,:);
        [old_method_bdry_cov]  = results_params.(results).(signal)(4,:);
        [subs]          = results_params.(results).(signal)(5,:);
        plot(subs,new_method_bdry_cov,[color_vec(2*j - 1) '-' 'x'],'linewidth', 1.5);
        plot(subs,old_method_bdry_cov,[color_vec(2*j) '-' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.2 1]);
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
        lgd = legend('t-bootstrap / BTSN simulation assessment', ...
                     't-bootstrap / SSS simulation assessment', ...
                     'wild bootstrap / BTSN simulation assessment', ...
                     'wild bootstrap / SSS simulation assessment', ...
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
export_fig(fh,fullfile(results_mat_dir,'Supplemental_3D_coverage_results.pdf'))

% Creating coverage plots for second signal model
figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    for j = 1:2
        signal = signal_vec(j+2);
        [new_method_bdry_cov] = results_params.(results).(signal)(2,:);
        [old_method_bdry_cov]  = results_params.(results).(signal)(4,:);
        [subs]          = results_params.(results).(signal)(5,:);
        plot(subs,new_method_bdry_cov,[color_vec(2*j - 1) '-' 'x'],'linewidth', 1.5);
        plot(subs,old_method_bdry_cov,[color_vec(2*j) '-' 'x'],'linewidth', 1.5);
    end 
    % plot the nominal lvl
    nominal_prct_level = results_params.(results).nominal_prct_level;
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level nominal_prct_level], 'k', 'linewidth', 1.5)
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level-1.96*results_params.(results).std_error nominal_prct_level-1.96*results_params.(results).std_error], 'k--')
    plot([subs(1)-5, subs(end)+5], [nominal_prct_level+1.96*results_params.(results).std_error nominal_prct_level+1.96*results_params.(results).std_error], 'k--')
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    ylim([0.2 1]);
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
        lgd = legend('t-Bootstrap / BTSN Simulation Assessment', ...
                     't-Bootstrap / SSS Simulation Assessment', ...
                     'Wild Bootstrap / BTSN Simulation Assessment', ...
                     'Wild Bootstrap / SSS Simulation Assessment', ...
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
export_fig(fh,fullfile(results_mat_dir,'Supplemental_2D_coverage_results.pdf'))