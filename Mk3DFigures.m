results_mat_dir = '/Users/maullz/Desktop/Contour_Inference_2018/Figures/3D_Figures';

results_mat_file = fullfile(results_mat_dir,'3DResults.mat');

load(results_mat_file)

nominal_vec = ["nom_80_results","nom_90_results","nom_95_results"];
signal_vec  = ["sig_1_std_1_results","sig_1_std_2_results","sig_2_std_1_results","sig_2_std_2_results","sig_3_std_1_results","sig_3_std_2_results","sig_4_results"];
color_vec   = 'rbrbrbr';

% Creating coverage plots for first signal model
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
    ylim([0.7 1]);
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
        lgd = legend('3D Sig. 1, Std Dev. 1 (est. boundary)', ...
                     '3D Sig. 1, Std Dev. 1 (true boundary)', ...
                     '3D Sig. 1, Std Dev. 2 (est. boundary)', ...
                     '3D Sig. 1, Std Dev. 2 (true boundary)', ...
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_1_coverage_results.pdf'))

% Creating coverage plots for second signal model
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
    ylim([0.7 1]);
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
        lgd = legend('3D Sig. 2, Std Dev. 1 (est. boundary)', ...
                     '3D Sig. 2, Std Dev. 1 (true boundary)', ...
                     '3D Sig. 2, Std Dev. 2 (est. boundary)', ...
                     '3D Sig. 2, Std Dev. 2 (true boundary)', ...
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_2_coverage_results.pdf'))

% Creating coverage plots for third signal model
figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
     
    for j = 5:6
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
    ylim([0.7 1]);
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
        lgd = legend('3D Sig. 3, Std Dev. 1 (est. boundary)', ...
                     '3D Sig. 3, Std Dev. 1 (true boundary)', ...
                     '3D Sig. 3, Std Dev. 2 (est. boundary)', ...
                     '3D Sig. 3, Std Dev. 2 (true boundary)', ...
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_3_coverage_results.pdf'))

% Creating coverage plots for fourth signal model
figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
      
    for j = 7
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
    ylim([0.7 1]);
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
        lgd = legend('3D Sig. 4 (est. boundary)', ...
                     '3D Sig. 4 (true boundary)', ...
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_4_coverage_results.pdf'))

% Creating quantile plots for first signal model
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
    
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('3D Sig. 1, Std Dev. 1 (est. boundary)', ...
                     '3D Sig. 1, Std Dev. 1 (true boundary)', ...
                     '3D Sig. 1, Std Dev. 2 (est. boundary)', ...
                     '3D Sig. 1, Std Dev. 2 (true boundary)');
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_1_quantile_results.pdf'))

% Creating quantile plots for second signal model
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
    
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('3D Sig. 2, Std Dev. 1 (est. boundary)', ...
                     '3D Sig. 2, Std Dev. 1 (true boundary)', ...
                     '3D Sig. 2, Std Dev. 2 (est. boundary)', ...
                     '3D Sig. 2, Std Dev. 2 (true boundary)');
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_2_quantile_results.pdf'))

% Creating quantile plots for third signal model
figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    for j = 5:6
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
    
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('3D Sig. 3, Std Dev. 1 (est. boundary)', ...
                     '3D Sig. 3, Std Dev. 1 (true boundary)', ...
                     '3D Sig. 3, Std Dev. 2 (est. boundary)', ...
                     '3D Sig. 3, Std Dev. 2 (true boundary)');
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_3_quantile_results.pdf'))

% Creating quantile plots for fourth signal model
figure, clf
for i = 1:length(nominal_vec)
    subplot(2,3,i);
    hold on
    results = nominal_vec(i);
    
    for j = 7
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
    
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    
    if i == length(nominal_vec)
        % create legend
        lgd = legend('3D Sig. 4 (est. boundary)', ...
                     '3D Sig. 4 (true boundary)');
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
export_fig(fh,fullfile(results_mat_dir,'3D_Sig_4_quantile_results.pdf'))