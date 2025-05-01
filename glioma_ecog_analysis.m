%% ECoG Power Spectral Density Analysis
%
% Author: Koen van der Kuil
% Affiliation: Department of Neuroscience, Erasmus University Medical Center, Rotterdam
%
% Description:
% This script analyzes the power spectral density (PSD) of ECoG recordings during seizure 
% and non-seizure periods across multiple brain regions. The analysis workflow:
%
% 1. Loads seizure metadata and recording information
% 2. Calculates PSDs for each seizure and non-seizure period using Welch's method
% 3. Computes band power in standard frequency bands (delta, theta, alpha, beta, gamma)
% 4. Performs statistical comparisons between seizure and non-seizure states
% 5. Generates visualizations including:
%    - Band power comparisons across channels
%    - Full PSD plots with frequency band highlighting
%    - Statistical significance indicators
%
% The results provide insights into how seizure activity affects spectral power
% across different frequency bands and brain regions.


% Import seizure data and all recording dates
seizures = readtable('.../seizures.csv');
% swds = readtable('.../swds.csv');
all_recordings = readtable('.../all_recording_dates.csv');

% Drop rows with NaN in Racine column and reset index
seizures = seizures(~isnan(seizures.Racine), :);

% Filter for validated seizures only
seizures = seizures(seizures.Validated == 1, :);

% Get actual column names from the table
actual_columns = seizures.Properties.VariableNames;
disp('Available columns in the seizures table:');
disp(actual_columns);

% Define columns to keep (excluding 'chan' and other unnecessary columns)
% Make sure all column names exist in the table
columns_to_keep = {'t_onset', 't_offset', 'date', 'dpi', 'label', ...
                   'duration', 't_onset_rounded', 'duration_rounded', 'RecID', ...
                   'trained', 'MouseID', 'ExperimentType', 'peramp', 'DPI_rounded', ...
                   'Validated', 'Racine', 'Description', 'Generalized_partial'};

% Check if all columns exist in the table
valid_columns = {};
for i = 1:length(columns_to_keep)
    if ismember(columns_to_keep{i}, actual_columns)
        valid_columns{end+1} = columns_to_keep{i};
    else
        warning(['Column "', columns_to_keep{i}, '" not found in the table.']);
    end
end
columns_to_keep = valid_columns;

% Create unique seizures by grouping by all columns
[~, idx] = unique(seizures(:, columns_to_keep), 'rows', 'first');
seizures_unique = seizures(idx, :);

% Display column names
disp(seizures_unique.Properties.VariableNames);

% Merge the dataframes (similar to pandas merge)
merged_df = outerjoin(all_recordings, seizures_unique, 'Keys', {'MouseID', 'dpi', 'ExperimentType', 'peramp'}, 'MergeKeys', true);

% Count seizures per recording
[G, TID] = findgroups(merged_df(:, {'MouseID', 'dpi'}));
seizure_counts = splitapply(@(x) sum(~isnan(x)), merged_df.t_onset, G);
seizure_count_table = [TID table(seizure_counts, 'VariableNames', {'seizure_count'})];

% Join the counts back to the merged dataframe
merged_df = join(merged_df, seizure_count_table, 'Keys', {'MouseID', 'dpi'});

% Fill NaN values with 0 for recordings without seizures
merged_df.seizure_count(isnan(merged_df.seizure_count)) = 0;

% Calculate standard error for error bars
[G, TID] = findgroups(merged_df(:, {'dpi', 'treated'}));
std_seizures = splitapply(@std, merged_df.seizure_count, G);
count_seizures = splitapply(@numel, merged_df.seizure_count, G);
sem_seizures = std_seizures ./ sqrt(count_seizures);
pivot_df_sem = [TID table(sem_seizures, 'VariableNames', {'sem_seizure_count'})];

% Round DPI to 1 day and trim between 20 and 50 DPI
merged_df.dpi_rounded = round(merged_df.dpi);
merged_df = merged_df((merged_df.dpi_rounded >= 20) & (merged_df.dpi_rounded <= 50), :);

% Create a pivot table for visualization
[G, TID] = findgroups(merged_df(:, {'dpi_rounded', 'treated'}));
avg_seizures = splitapply(@mean, merged_df.seizure_count, G);
pivot_df = [TID table(avg_seizures, 'VariableNames', {'avg_seizure_count'})];

% Define colors
colordict = containers.Map({'Treated', 'No treatment'}, {[0.1216, 0.4667, 0.7059], [0, 0, 0]});


%%

% Create a figure
figure('Position', [100, 100, 900, 750], 'Color', 'white');
ax = gca;

% Get unique rounded DPIs
unique_dpis = unique(merged_df.dpi_rounded);
unique_dpis = sort(unique_dpis);

% Create bar positions
bar_width = 0.35;
x_pos = 1:length(unique_dpis);

% Plot bars for untreated group
untreated_means = zeros(size(unique_dpis));
untreated_sems = zeros(size(unique_dpis));
for i = 1:length(unique_dpis)
    dpi = unique_dpis(i);
    untreated_data = merged_df((merged_df.dpi_rounded == dpi) & (merged_df.treated == 0), :);
    untreated_means(i) = mean(untreated_data.seizure_count, 'omitnan');
    untreated_sems(i) = std(untreated_data.seizure_count, 'omitnan') / sqrt(sum(~isnan(untreated_data.seizure_count)));
end
h1 = bar(x_pos - bar_width/2, untreated_means, bar_width, 'FaceColor', 'none', 'EdgeColor', colordict('No treatment'), 'LineWidth', 1);
hold on;

% Plot error bars for untreated
errorbar(x_pos - bar_width/2, untreated_means, untreated_sems, '.', 'Color', colordict('No treatment'), 'LineWidth', 0.75);

% Plot bars for treated group
treated_means = zeros(size(unique_dpis));
treated_sems = zeros(size(unique_dpis));
for i = 1:length(unique_dpis)
    dpi = unique_dpis(i);
    treated_data = merged_df((merged_df.dpi_rounded == dpi) & (merged_df.treated == 1), :);
    treated_means(i) = mean(treated_data.seizure_count, 'omitnan');
    treated_sems(i) = std(treated_data.seizure_count, 'omitnan') / sqrt(sum(~isnan(treated_data.seizure_count)));
end
h2 = bar(x_pos + bar_width/2, treated_means, bar_width, 'FaceColor', 'none', 'EdgeColor', colordict('Treated'), 'LineWidth', 1);

% Plot error bars for treated
errorbar(x_pos + bar_width/2, treated_means, treated_sems, '.', 'Color', colordict('Treated'), 'LineWidth', 0.75);

% Customize the plot
title('Average Number of Seizures Over Time by Treatment Group', 'FontSize', 9, 'FontName', 'Arial');
xlabel('Days Post Injection (DPI)', 'FontSize', 8, 'FontName', 'Arial');
ylabel('Average Number of Seizures', 'FontSize', 8, 'FontName', 'Arial');

% Adjust x-axis labels to show every 5th day
xticks_pos = x_pos(mod(unique_dpis, 5) == 0);
xticks_labels = cellstr(num2str(unique_dpis(mod(unique_dpis, 5) == 0)));
set(ax, 'XTick', xticks_pos, 'XTickLabel', xticks_labels);
xtickangle(45);

% Add significance asterisks
for i = 1:length(unique_dpis)
    dpi = unique_dpis(i);
    untreated = merged_df((merged_df.treated == 0) & (merged_df.dpi_rounded == dpi), :).seizure_count;
    treated = merged_df((merged_df.treated == 1) & (merged_df.dpi_rounded == dpi), :).seizure_count;
    [~, p_value] = ttest2(untreated, treated, 'Vartype', 'unequal');
    if p_value < 0.05
        disp(dpi);
        text(x_pos(i), max([untreated_means(i), treated_means(i)])*1.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
end

% Adjust y-axis to start from 0
ylim([0, max([untreated_means + untreated_sems; treated_means + treated_sems])*1.3]);

% Add a text box with overall statistics
treated_seizures = merged_df.seizure_count(merged_df.treated == 1);
untreated_seizures = merged_df.seizure_count(merged_df.treated == 0);
[~, p_value, ~, stats] = ttest2(treated_seizures, untreated_seizures);
stats_text = sprintf('Overall:\nt-statistic: %.2f\np-value: %.4f', stats.tstat, p_value);
text(0.5, 0.95, stats_text, 'Units', 'normalized', 'FontSize', 6, 'FontName', 'Arial', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

% Modify the legend
legend([h1, h2], {'No treatment', 'Treated'}, 'FontSize', 4, 'FontName', 'Arial', 'Box', 'off');

% Remove top and right spines
box off;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
ax.TickLength = [0.02 0.02];
set(ax, 'FontName', 'Arial', 'FontSize', 8);

% Save the figure
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3 2.5]);
print('F:/ECoG/03_Results/seizures_by_dpi_and_treatment.pdf', '-dpdf', '-r300');

%%

%% Create partial and general seizure dataframes
% Filter for validated seizures with Racine scores
valid_seizures = seizures_unique(~isnan(seizures_unique.Racine) & seizures_unique.Validated == 1, :);

% Separate partial and generalized seizures based on Racine score
partial_seizures = valid_seizures(valid_seizures.Racine <= 3, :);
partial_seizures = partial_seizures(~isnan(partial_seizures.Racine), :);

general_seizures = valid_seizures(valid_seizures.Racine > 3, :);
general_seizures = general_seizures(~isnan(general_seizures.Racine), :);

% Define frequency bands
bands = {'delta', 'alpha', 'theta', 'beta', 'gamma'};
bands_freqs = {[1, 4], [4, 8], [8, 12], [12, 30], [30, 45]};
sfreq = 512;

% Create a table to store power data
general_power_df = table();
% Create a structure to store power spectra for later plotting
power_spectra = struct();
power_spectra.seizure = struct();
power_spectra.nonseizure = struct();
power_spectra.metadata = struct();
power_spectra.freq = [];

% Process each seizure in the general seizures dataframe
fprintf('Processing generalized seizures: 0/%d (0%%)', height(general_seizures));
for i = 1:height(general_seizures)
    % Update progress bar
    fprintf('\rProcessing generalized seizures: %d/%d (%d%%)', i, height(general_seizures), round(i/height(general_seizures)*100));
    
    try
        recording_id = general_seizures.RecID{i};
        mouse_id = strsplit(recording_id, '_');
        mouse_id = mouse_id{1};
        
        % Get seizure timing information
        t_onset = general_seizures.t_onset(i);
        t_offset = general_seizures.t_offset(i);
        
        seizure_start = round(t_onset * sfreq);
        seizure_end = round(t_offset * sfreq);
        seizure_duration = round((t_offset - t_onset) * sfreq);
        
        % Load EDF file
        edf_path = fullfile('D:/Glioma Data Backup/Pten_ECoG/data/', mouse_id, recording_id);
        edf_files = dir(fullfile(edf_path, '*.edf'));
        if isempty(edf_files)
            warning('No EDF file found for recording %s', recording_id);
            continue;
        end
        intan_filename = fullfile(edf_files(1).folder, edf_files(1).name);
        
        % Load and preprocess EDF
        [data, hdr] = edfread(intan_filename);
        sfreq = 512; % hdr.frequency(1);  % Assuming all channels have the same sampling rate
        
        % Get channel names from data variable properties instead of hdr.label
        ch_names = data.Properties.VariableNames;
        
        % Define channel renaming dictionary
        rename_channel_dict = containers.Map({'Ipsi_frontal', 'Ipsi_parietal', 'Contra_frontal', 'Cerebellar'}, ...
                                           {'Ipsi_Frontal', 'Ipsi_Parietal', 'Contra_Frontal', 'Cerebellar'});
        
        % Add buffer around seizure
        %seizure_start = max(1, seizure_start - (30 * sfreq));
        %seizure_end = min(height(data), seizure_end + (30 * sfreq));
        total_duration = height(data);
        
        % Process each channel
        for ch_idx = 1:length(ch_names)
            chan = ch_names{ch_idx};
            
            % Extract full channel data - concatenate all blocks
            chan_data = [];
            for block_idx = 1:length(data.(chan))
                chan_data = [chan_data; data.(chan){block_idx}];
            end
            
            % Get random non-seizure period of same duration
            non_seizure_duration = seizure_end - seizure_start;
            valid_start = [1, seizure_start - non_seizure_duration, seizure_end + 1];
            valid_end = [seizure_start - 1, seizure_end, length(chan_data)];
            
            % Select a random valid segment
            valid_segments = find((valid_end - valid_start) >= non_seizure_duration);
            if isempty(valid_segments)
                warning('No valid non-seizure segment found for recording %s', recording_id);
                continue;
            end
            segment_idx = valid_segments(randi(length(valid_segments)));
            window_start = valid_start(segment_idx) + randi(valid_end(segment_idx) - valid_start(segment_idx) - non_seizure_duration);
            window_stop = window_start + non_seizure_duration;
            
            % Extract seizure and non-seizure data
            window_seizure = chan_data(seizure_start:seizure_end);
            window_nonseizure = chan_data(window_start:window_stop);
            
            % Get standardized channel name for storage
            std_chan = chan;
            if isKey(rename_channel_dict, chan)
                std_chan = rename_channel_dict(chan);
            end
            
            % Calculate band power for seizure and non-seizure periods
            for b = 1:length(bands)
                band_name = bands{b};
                freq_range = bands_freqs{b};
                
                % Calculate power for seizure period - ensure window size is appropriate
                window_size = min(round(sfreq*2), floor(length(window_seizure)/2));
                if window_size < 3  % Skip if window is too small
                    warning('Seizure window too small for spectral analysis in recording %s', recording_id);
                    continue;
                end
                [pxx_seizure, f] = pwelch(window_seizure, hamming(window_size), round(window_size/2), [], sfreq);
                idx_band = f >= freq_range(1) & f <= freq_range(2);
                band_power_seizure = mean(pxx_seizure(idx_band));
                
                % Store full power spectrum for each seizure separately with metadata
                seizure_id = sprintf('seizure_%d', i);
                if ~isfield(power_spectra.seizure, seizure_id)
                    power_spectra.seizure.(seizure_id) = struct();
                    power_spectra.metadata.(seizure_id) = struct();
                    
                    % Store metadata for this seizure
                    power_spectra.metadata.(seizure_id).recording_id = recording_id;
                    power_spectra.metadata.(seizure_id).mouse_id = mouse_id;
                    power_spectra.metadata.(seizure_id).dpi = general_seizures.dpi(i);
                    power_spectra.metadata.(seizure_id).racine = general_seizures.Racine(i);
                    power_spectra.metadata.(seizure_id).t_onset = t_onset;
                    power_spectra.metadata.(seizure_id).t_offset = t_offset;
                    power_spectra.metadata.(seizure_id).duration = t_offset - t_onset;
                end
                
                if ~isfield(power_spectra.seizure.(seizure_id), std_chan)
                    power_spectra.seizure.(seizure_id).(std_chan) = pxx_seizure;
                end
                
                % Store corresponding non-seizure data
                nonseizure_id = sprintf('nonseizure_%d', i);
                if ~isfield(power_spectra.nonseizure, nonseizure_id)
                    power_spectra.nonseizure.(nonseizure_id) = struct();
                end
                
                if ~isfield(power_spectra.nonseizure.(nonseizure_id), std_chan)
                    % Calculate power for non-seizure period
                    window_size = min(round(sfreq*2), floor(length(window_nonseizure)/2));
                    if window_size < 3  % Skip if window is too small
                        warning('Non-seizure window too small for spectral analysis in recording %s', recording_id);
                        continue;
                    end
                    [pxx_nonseizure, ~] = pwelch(window_nonseizure, hamming(window_size), round(window_size/2), [], sfreq);
                    band_power_nonseizure = mean(pxx_nonseizure(idx_band));
                    
                    % Store full power spectrum for non-seizure
                    power_spectra.nonseizure.(nonseizure_id).(std_chan) = pxx_nonseizure;
                else
                    % If already calculated, just get the band power
                    pxx_nonseizure = power_spectra.nonseizure.(nonseizure_id).(std_chan);
                    band_power_nonseizure = mean(pxx_nonseizure(idx_band));
                end
                
                % Store frequency vector once
                if isempty(power_spectra.freq)
                    power_spectra.freq = f;
                end
                
                % Add to table - use cell arrays for string data
                seizure_row = table({recording_id}, {mouse_id}, general_seizures.dpi(i), ...
                                   {chan}, {band_name}, band_power_seizure, {'Generalized'}, ...
                                   'VariableNames', {'rec_id', 'mouse_id', 'dpi', 'chan', ...
                                                    'band', 'power', 'period'});
                                                
                nonseizure_row = table({recording_id}, {mouse_id}, general_seizures.dpi(i), ...
                                      {chan}, {band_name}, band_power_nonseizure, {'No seizure'}, ...
                                      'VariableNames', {'rec_id', 'mouse_id', 'dpi', 'chan', ...
                                                       'band', 'power', 'period'});
                
                general_power_df = [general_power_df; seizure_row; nonseizure_row];
            end
        end
    catch e
        warning('Error processing recording %s: %s', recording_id, e.message);
        disp(['Stack trace: ' getReport(e)]);
        continue;
    end
end
fprintf('\rProcessing generalized seizures: %d/%d (100%%)\n', height(general_seizures), height(general_seizures));

% Map channel names to standard names
fprintf('Mapping channel names: 0/%d (0%%)', height(general_power_df));
for i = 1:height(general_power_df)
    if mod(i, 100) == 0 || i == height(general_power_df)
        fprintf('\rMapping channel names: %d/%d (%d%%)', i, height(general_power_df), round(i/height(general_power_df)*100));
    end
    
    chan = general_power_df.chan{i};
    if isKey(rename_channel_dict, chan)
        general_power_df.chan{i} = rename_channel_dict(chan);
    end
end
fprintf('\rMapping channel names: %d/%d (100%%)\n', height(general_power_df), height(general_power_df));

% Add order column for plotting
fprintf('Adding order column: 0/%d (0%%)', height(general_power_df));
order_map = containers.Map({'Ipsi Frontal', 'Ipsi Parietal', 'Contra Frontal', 'Cerebellar'}, {0, 1, 2, 3});
general_power_df.order = zeros(height(general_power_df), 1);
for i = 1:height(general_power_df)
    if mod(i, 100) == 0 || i == height(general_power_df)
        fprintf('\rAdding order column: %d/%d (%d%%)', i, height(general_power_df), round(i/height(general_power_df)*100));
    end
    
    chan = general_power_df.chan{i};
    if isKey(order_map, chan)
        general_power_df.order(i) = order_map(chan);
    end
end
fprintf('\rAdding order column: %d/%d (100%%)\n', height(general_power_df), height(general_power_df));

% Save power spectra for later use
save('power_spectra_data.mat', 'power_spectra');

%% Plot bandpower data by frequency band with channels on x-axis and seizure state as color
figure('Position', [100, 100, 1200, 800]);

% Get unique bands for plotting
unique_bands = unique(general_power_df.band);
num_bands = length(unique_bands);

% Get unique channels for x-axis
unique_channels = unique(general_power_df.chan);
num_channels = length(unique_channels);

% Define colors for seizure vs non-seizure
seizure_color = [0.8500, 0.3250, 0.0980];  % Orange
nonseizure_color = [0, 0.4470, 0.7410];    % Blue

% Create a subplot for each frequency band
for b = 1:num_bands
    subplot(2, 3, b);
    
    % Filter data for current band
    current_band = unique_bands{b};
    band_data = general_power_df(strcmp(general_power_df.band, current_band), :);
    
    % Initialize arrays to store means and errors
    seizure_means = zeros(1, num_channels);
    seizure_sems = zeros(1, num_channels);
    nonseizure_means = zeros(1, num_channels);
    nonseizure_sems = zeros(1, num_channels);
    
    % Calculate statistics for each channel
    for c = 1:num_channels
        channel = unique_channels{c};
        channel_data = band_data(strcmp(band_data.chan, channel), :);
        
        % Group by period (Generalized vs No seizure)
        seizure_data = channel_data(strcmp(channel_data.period, 'Generalized'), :);
        nonseizure_data = channel_data(strcmp(channel_data.period, 'No seizure'), :);
        
        % Calculate means and standard errors
        if ~isempty(seizure_data)
            seizure_means(c) = mean(seizure_data.power);
            seizure_sems(c) = std(seizure_data.power) / sqrt(height(seizure_data));
        end
        
        if ~isempty(nonseizure_data)
            nonseizure_means(c) = mean(nonseizure_data.power);
            nonseizure_sems(c) = std(nonseizure_data.power) / sqrt(height(nonseizure_data));
        end
    end
    
    % Set up bar positions
    bar_width = 0.35;
    x_pos = 1:num_channels;
    
    % Plot grouped bar chart
    hold on;
    h1 = bar(x_pos - bar_width/2, seizure_means, bar_width, 'FaceColor', seizure_color);
    h2 = bar(x_pos + bar_width/2, nonseizure_means, bar_width, 'FaceColor', nonseizure_color);
    
    % Add error bars
    errorbar(x_pos - bar_width/2, seizure_means, seizure_sems, 'k', 'LineStyle', 'none', 'CapSize', 5);
    errorbar(x_pos + bar_width/2, nonseizure_means, nonseizure_sems, 'k', 'LineStyle', 'none', 'CapSize', 5);
    
    % Set labels and title
    title([upper(current_band(1)) current_band(2:end) ' Band Power'], 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Power (μV²/Hz)', 'FontSize', 10, 'FontName', 'Arial');
    xlabel('Channel', 'FontSize', 10, 'FontName', 'Arial');
    
    % Set x-axis labels
    xticks(x_pos);
    xticklabels(unique_channels);
    xtickangle(45);
    
    % Add legend
    if b == 1
        legend([h1, h2], {'Seizure', 'No Seizure'}, 'Location', 'best', 'FontSize', 8);
    end
    
    % Add statistical comparison
    for c = 1:num_channels
        % Perform t-test between seizure and non-seizure for this channel
        channel = unique_channels{c};
        channel_data = band_data(strcmp(band_data.chan, channel), :);
        
        seizure_values = channel_data.power(strcmp(channel_data.period, 'Generalized'));
        nonseizure_values = channel_data.power(strcmp(channel_data.period, 'No seizure'));
        
        if ~isempty(seizure_values) && ~isempty(nonseizure_values)
            [~, p] = ttest2(seizure_values, nonseizure_values);
            if p < 0.05
                y_pos = max([seizure_means(c) + seizure_sems(c), nonseizure_means(c) + nonseizure_sems(c)]) * 1.1;
                text(c, y_pos, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
            end
        end
    end
    
    % Adjust appearance
    box off;
    grid on;
    hold off;
    
    % Set y-axis limits to be consistent across subplots
    all_values = [seizure_means + seizure_sems, nonseizure_means + nonseizure_sems];
    ylim([0, max(all_values) * 1.2]);
end

% Adjust overall figure appearance
sgtitle('Band Power Comparison by Channel and Seizure State', 'FontSize', 14, 'FontName', 'Arial');

% Save the figure
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 8]);
%print('../bandpower_by_channel.pdf', '-dpdf', '-r300');

%% Plot Power Spectral Density (PSD) for seizure vs non-seizure periods
figure('Position', [100, 100, 1000, 800]);
% Get the structure fields for seizure and non-seizure data
seizure_fields = fieldnames(power_spectra.seizure);
num_seizures = length(seizure_fields);

% Define colors for seizure vs non-seizure
seizure_color = [0.8500, 0.3250, 0.0980];  % Orange
nonseizure_color = [0, 0.4470, 0.7410];    % Blue

% Get frequency vector
freq = power_spectra.freq;

% First, determine the channels from the first seizure
first_seizure = power_spectra.seizure.(seizure_fields{1});
channel_names = fieldnames(first_seizure);
num_channels = length(channel_names);

% Create a subplot for each channel
for ch = 1:num_channels
    channel = channel_names{ch};
    
    % Create subplot
    subplot(num_channels, 1, ch);
    
    % Get the PSD data for this channel directly
    % Instead of averaging across events, plot each event separately
    hold on;
    
    % Plot seizure PSDs
    seizure_legend_added = false;
    for s = 1:num_seizures
        seizure_field = seizure_fields{s};
        if isfield(power_spectra.seizure.(seizure_field), channel) && ...
           isnumeric(power_spectra.seizure.(seizure_field).(channel))
            
            seizure_psd = power_spectra.seizure.(seizure_field).(channel);
            
            % Plot with semi-transparency to see overlapping lines
            if ~seizure_legend_added
                h1 = plot(freq, seizure_psd, 'Color', [seizure_color, 0.7], 'LineWidth', 1);
                seizure_legend_added = true;
            else
                plot(freq, seizure_psd, 'Color', [seizure_color, 0.7], 'LineWidth', 1);
            end
        end
    end
    
    % Plot non-seizure PSDs
    nonseizure_legend_added = false;
    nonseizure_fields = fieldnames(power_spectra.nonseizure);
    for ns = 1:length(nonseizure_fields)
        nonseizure_field = nonseizure_fields{ns};
        if isfield(power_spectra.nonseizure.(nonseizure_field), channel) && ...
           isnumeric(power_spectra.nonseizure.(nonseizure_field).(channel))
            
            nonseizure_psd = power_spectra.nonseizure.(nonseizure_field).(channel);
            
            % Plot with semi-transparency
            if ~nonseizure_legend_added
                h2 = plot(freq, nonseizure_psd, 'Color', [nonseizure_color, 0.7], 'LineWidth', 1);
                nonseizure_legend_added = true;
            else
                plot(freq, nonseizure_psd, 'Color', [nonseizure_color, 0.7], 'LineWidth', 1);
            end
        end
    end
    
    % Add frequency band shading
    band_colors = [0.9 0.9 0.9; 0.85 0.85 0.85; 0.9 0.9 0.9; 0.85 0.85 0.85; 0.9 0.9 0.9];
    for b = 1:length(bands)
        freq_range = bands_freqs{b};
        idx = freq >= freq_range(1) & freq <= freq_range(2);
        if any(idx)
            y_lim = get(gca, 'YLim');
            if isempty(y_lim)
                y_lim = [0, 1]; % Default if no data plotted yet
            end
            
            x_patch = [freq_range(1), freq_range(2), freq_range(2), freq_range(1)];
            y_patch = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
            patch(x_patch, y_patch, band_colors(b,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            
            % Add band name
            text(mean(freq_range), y_lim(2)*0.9, bands{b}, 'HorizontalAlignment', 'center', ...
                 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    
    % Set title and labels
    title(['Channel: ' strrep(channel, '_', ' ')], 'FontSize', 12, 'FontName', 'Arial');
    
    if ch == num_channels  % Only add x-label to bottom subplot
        xlabel('Frequency (Hz)', 'FontSize', 10, 'FontName', 'Arial');
    else
        set(gca, 'XTickLabel', []);  % Remove x-tick labels for all but bottom subplot
    end
    
    ylabel('Power (μV²/Hz)', 'FontSize', 10, 'FontName', 'Arial');
    
    % Set x-axis limits
    xlim([0, 45]);  % Limit to frequencies of interest
    
    % Add legend to first subplot only
    if ch == 1 && exist('h1', 'var') && exist('h2', 'var')
        legend([h1, h2], {'Seizure', 'No Seizure'}, 'Location', 'northeast', 'FontSize', 8);
    end
    
    % Customize appearance
    box off;
    grid on;
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    
    % Use log scale for y-axis to better visualize differences
    set(gca, 'YScale', 'log');
    
    hold off;
end

% Adjust spacing between subplots
tight_spacing = 0.03;
subplot_height = (1 - (num_channels+1)*tight_spacing) / num_channels;
for ch = 1:num_channels
    subplot_position = [(tight_spacing), ...
                        (1 - ch*subplot_height - ch*tight_spacing), ...
                        (1 - 2*tight_spacing), ...
                        subplot_height];
    set(subplot(num_channels, 1, ch), 'Position', subplot_position);
end

% Add overall title
sgtitle('Power Spectral Density: Seizure vs. Non-Seizure Periods', 'FontSize', 14, 'FontName', 'Arial');

% Save the figure
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 8]);
%print(../psd_seizure_vs_nonseizure.pdf', '-dpdf', '-r300');