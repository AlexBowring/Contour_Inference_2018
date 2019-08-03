results_mat_dir = '/Users/maullz/Desktop/Contour_Inference_2018/Figures/2D_Figures';

results_mat_file = fullfile(results_mat_dir,'2DResults.mat');

load(results_mat_file)

nominal_vec = ["nom_80_results","nom_90_results","nom_95_results"];
signal_vec  = ["sig_1_std_1_results","sig_1_std_2_results","sig_2_std_1_results","sig_2_std_2_results"];
color_vec   = 'rbrb';

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

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
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Ramp (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
        lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'2D_Sig_1_coverage_results.pdf'))

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
    ylabel('Emp. Covering Rate');
    
    titlename = sprintf('Circle (%d%% Nom.)', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('Homo. Variance (est. bdry)', ...
                     'Homo. Variance (true bdry)', ...
                     'Hetro. Variance (est. bdry)', ...
                     'Hetro. Variance (true bdry)', ...
                     'Nominal Coverage Level', ...
                     '\pm 1.96 \times Simulation Std. Error');
         lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'2D_Sig_2_coverage_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    for j = 1:2
        signal = signal_vec(j);
        [true_bdry_quant] = results_params.(results).(signal)(3,:);
        [est_bdry_quant]  = results_params.(results).(signal)(4,:);
        [subs]          = results_params.(results).(signal)(5,:);
        plot(subs,est_bdry_quant,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_quant,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    % ylim([0.5 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Quantile Value');
    
    titlename = sprintf('%d%% Quantile Results', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('2D Sig. 1, Std Dev. 1 (est. boundary)', ...
                     '2D Sig. 1, Std Dev. 1 (true boundary)', ...
                     '2D Sig. 1, Std Dev. 2 (est. boundary)', ...
                     '2D Sig. 1, Std Dev. 2 (true boundary)');
        lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'2D_Sig_1_quantile_results.pdf'))

figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    for j = 3:4
        signal = signal_vec(j);
        [true_bdry_quant] = results_params.(results).(signal)(3,:);
        [est_bdry_quant]  = results_params.(results).(signal)(4,:);
        [subs]          = results_params.(results).(signal)(5,:);
        plot(subs,est_bdry_quant,[color_vec(j) '-' 'x'],'linewidth', 1.5);
        plot(subs,true_bdry_quant,[color_vec(j) '--' 'x'],'linewidth', 1.5);
    end 
    % set the range for the axis
    xlim([subs(1)-5, subs(end)+5])
    % ylim([0.5 1]);
    % specify tiks
    xticks(subs)
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('Quantile Value');
    
    titlename = sprintf('%d%% Quantile Results', results_params.(results).nominal_level);
    title(titlename);
    
    set(gca, 'fontsize', 18);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('2D Sig. 2, Std Dev. 1 (est. boundary)', ...
                     '2D Sig. 2, Std Dev. 1 (true boundary)', ...
                     '2D Sig. 2, Std Dev. 2 (est. boundary)', ...
                     '2D Sig. 2, Std Dev. 2 (true boundary)');
        lgd.FontSize = 18;
    end   
        
end

lgd_plot = subplot(2,3,5);
axis square;
pos_lgd  = get(lgd_plot,'position');
lgd.FontWeight = 'bold';
set(lgd,'position', [pos_lgd(1), pos_lgd(2) + 0.2, pos_lgd(3), pos_lgd(4) - 0.2]);
set(lgd, 'interpreter', 'tex');
axis(lgd_plot,'off');

set(gcf,'position', [-21,120,1195,682]);
fh = gcf;
set(fh,'color','w');
export_fig(fh,fullfile(results_mat_dir,'2D_Sig_2_quantile_results.pdf'))