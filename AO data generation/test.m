% Compute absolute errors
snr_seperate_all_vector_error = snr_gen_all_vector - snr_seperate_all_vector;
snr_togather_all_vector_error = snr_gen_all_vector - snr_togather_all_vector;
snr_hybrid_all_vector_error = snr_gen_all_vector - snr_hybrid_all_vector;
snr_nonml_all_vector_error = snr_gen_all_vector - snr_nonml_all_vector;

% Compute percentage errors
epsilon = 1e-6;
snr_gen_safe = snr_gen_all_vector;
snr_gen_safe(snr_gen_safe == 0) = epsilon;

snr_seperate_pct_error = (snr_seperate_all_vector_error ./ snr_gen_safe) * 100;
snr_togather_pct_error = (snr_togather_all_vector_error ./ snr_gen_safe) * 100;
snr_hybrid_pct_error = (snr_hybrid_all_vector_error ./ snr_gen_safe) * 100;
snr_nonml_pct_error = (snr_nonml_all_vector_error ./ snr_gen_safe) * 100;

% Visualization settings
fontsize = 16;
linewidth = 3;
boxcolor = [0.2 0.4 0.6]; % dark bluish tone

% Helper function
function plot_styled_boxplot(data, title_text, fontsize, linewidth, boxcolor)
    figure;
    h = boxplot(data, 'Labels', string(1:16), 'Widths', 0.5, 'Colors', boxcolor);
    xlabel('Iteration', 'FontSize', fontsize, 'FontWeight', 'bold');
    ylabel('Weighted Sum Rate Error (%)', 'FontSize', fontsize, 'FontWeight', 'bold');
    title(title_text, 'FontSize', fontsize + 2, 'FontWeight', 'bold');
    set(gca, 'FontSize', fontsize, 'LineWidth', linewidth-1.5);
    ylim([0 65]);
    grid on;

    % Update line width for boxes, whiskers, medians, etc.
    set(h, 'LineWidth', linewidth-1.3);
end

% Generate box plots
plot_styled_boxplot(snr_seperate_pct_error, 'Percentage Error for Approach 1(Seperate Predicton)', fontsize, linewidth, boxcolor);
plot_styled_boxplot(snr_togather_pct_error, 'Percentage Error for Approach 2(Joint Prediction)', fontsize, linewidth, boxcolor);
plot_styled_boxplot(snr_hybrid_pct_error, 'Percentage Error for Approach 3(Beaforming predction followed by phase predciton)', fontsize, linewidth, boxcolor);
plot_styled_boxplot(snr_nonml_pct_error, 'Percentage Error for Random Intialiazation', fontsize, linewidth, boxcolor);
