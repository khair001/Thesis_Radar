clear all; 
close all; 
clc;

%% --- CONFIGURATION AND SETUP ---
BASELINE_DATA_FILE = 'radar_logs/baseline500frame_test_data.mat';
INTERFERENCE_DATA_FILE = 'radar_logs/interference500frame_test_data.mat';
BASELINE_VAR_NAME = 'allBaselineObjects';
INTERFERENCE_VAR_NAME = 'allInterferenceObjects';

% --- DBSCAN CLUSTERING PARAMETERS ---
DBSCAN_EPSILON = 0.4; 
DBSCAN_MINPTS = 6;

% --- CLASSIFICATION PARAMETERS ---
% Any object with an absolute average velocity below this value (in m/s) will be
% classified as 'static'.
STATIC_VELOCITY_THRESHOLD = 0.05; % (m/s)

% --- VISUALIZATION SETTINGS ---
MAX_RANGE_DISPLAY = 5.0; % Consistent axis limits for visualization (meters)
NOISE_COLOR = [0.8 0.8 0.8]; % Gray for noise

%% --- LOAD DATA ---
disp('--- Loading Radar Data for Analysis ---');
baselineObjects = [];
interferenceObjects = [];

% Load Baseline Data
if exist(BASELINE_DATA_FILE, 'file')
    try
        loadedData = load(BASELINE_DATA_FILE);
        if isfield(loadedData, BASELINE_VAR_NAME)
            baselineObjects = loadedData.(BASELINE_VAR_NAME);
            disp(['- Successfully loaded ', num2str(length(baselineObjects)), ' raw baseline detections.']);
        else
            warning(['Could not find variable ''', BASELINE_VAR_NAME, ''' in ', BASELINE_DATA_FILE, '. Attempting fallback...']);
            fieldNames = fieldnames(loadedData);
            for i = 1:length(fieldNames)
                if isstruct(loadedData.(fieldNames{i}))
                    baselineObjects = loadedData.(fieldNames{i});
                    disp(['- Loaded fallback variable ''', fieldNames{i}, ''' with ', num2str(length(baselineObjects)), ' raw detections.']);
                    break;
                end
            end
        end
    catch ME
        disp(['ERROR: Failed to load baseline data: ', ME.message]);
    end
else
    disp(['WARNING: Baseline data file not found: ', BASELINE_DATA_FILE]);
end

% Load Interference Data
if exist(INTERFERENCE_DATA_FILE, 'file')
    try
        loadedData = load(INTERFERENCE_DATA_FILE);
        if isfield(loadedData, INTERFERENCE_VAR_NAME)
            interferenceObjects = loadedData.(INTERFERENCE_VAR_NAME);
            disp(['- Successfully loaded ', num2str(length(interferenceObjects)), ' raw interference detections.']);
        else
            warning(['Could not find variable ''', INTERFERENCE_VAR_NAME, ''' in ', INTERFERENCE_DATA_FILE, '. Attempting fallback...']);
            fieldNames = fieldnames(loadedData);
            for i = 1:length(fieldNames)
                if isstruct(loadedData.(fieldNames{i}))
                    interferenceObjects = loadedData.(fieldNames{i});
                    disp(['- Loaded fallback variable ''', fieldNames{i}, ''' with ', num2str(length(interferenceObjects)), ' raw detections.']);
                    break;
                end
            end
        end
    catch ME
        disp(['ERROR: Failed to load interference data: ', ME.message]);
    end
else
    disp(['WARNING: Interference data file not found: ', INTERFERENCE_DATA_FILE]);
end

if isempty(baselineObjects) && isempty(interferenceObjects)
    disp('No data loaded. Analysis cannot proceed.');
    return;
end
disp('--- Data Loading Complete. Starting Analysis ---');

%% --- PERFORM DBSCAN CLUSTERING, CLASSIFICATION, AND METRICS ---
fprintf('\n------------------------------------------------\n');
disp('--- Analysis: DBSCAN Clustering & Summary ---');

% Cluster and classify baseline detections
[baselineClusters, numStaticBaseline, numMovingBaseline, baselineOutlierIdx] = ...
    clusterAndClassify(baselineObjects, DBSCAN_EPSILON, DBSCAN_MINPTS, STATIC_VELOCITY_THRESHOLD);

% Cluster and classify interference detections
[interferenceClusters, numStaticInterference, numMovingInterference, interferenceOutlierIdx] = ...
    clusterAndClassify(interferenceObjects, DBSCAN_EPSILON, DBSCAN_MINPTS, STATIC_VELOCITY_THRESHOLD);

% Calculate Performance Metrics
fprintf('--- Performance Metrics ---\n');
% Baseline Metrics
num_total_baseline = length(baselineObjects);
num_noise_baseline = sum(baselineOutlierIdx);
num_clusters_baseline = length(baselineClusters);
fpr_baseline = num_noise_baseline / num_total_baseline;

% Interference Metrics
num_total_interference = length(interferenceObjects);
num_noise_interference = sum(interferenceOutlierIdx);
num_clusters_interference = length(interferenceClusters);
fpr_interference = num_noise_interference / num_total_interference;

% Print Summary
fprintf('------------------------------------------------\n');
fprintf('Analysis for BASELINE Data:\n');
fprintf('Total Raw Detections: %d\n', num_total_baseline);
fprintf('Estimated Objects (Clustered): %d\n', num_clusters_baseline);
fprintf('  - Static Objects: %d\n', numStaticBaseline);
fprintf('  - Moving Objects: %d\n', numMovingBaseline);
fprintf('  - Noise Detections: %d\n', num_noise_baseline);
fprintf('Performance Metrics:\n');
fprintf('  - False Positive Rate (FPR): %.2f%%\n', fpr_baseline * 100);

fprintf('\nAnalysis for INTERFERENCE Data:\n');
fprintf('Total Raw Detections: %d\n', num_total_interference);
fprintf('Estimated Objects (Clustered): %d\n', num_clusters_interference);
fprintf('  - Static Objects: %d\n', numStaticInterference);
fprintf('  - Moving Objects: %d\n', numMovingInterference);
fprintf('  - Noise Detections: %d\n', num_noise_interference);
fprintf('Performance Metrics:\n');
fprintf('  - False Positive Rate (FPR): %.2f%%\n', fpr_interference * 100);
fprintf('------------------------------------------------\n');

%% --- VISUALIZATION: 3D SCATTER PLOTS (Combined in Subplots) ---
fprintf('\n--- Visualizing Clustered 3D Spatial Distribution ---\n');
% Create one figure with 2 subplots
figure('Name', '3D Radar Object Visualization (Static vs Moving)', 'NumberTitle', 'off', 'Position', [50 50 1400 700]);

%% ------------------------------------------------------------------------
% Subplot 1: Baseline Objects
ax_baseline = subplot(1,2,1); % left subplot
view(ax_baseline, 3);
grid(ax_baseline, 'on'); hold(ax_baseline, 'on');
h_radar = plot3(ax_baseline, 0, 0, 0, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'blue', 'DisplayName', 'Radar');
h_noise = plot3(ax_baseline, nan, nan, nan, 'o', 'MarkerFaceColor', NOISE_COLOR, 'MarkerEdgeColor', 'k', 'MarkerSize', 3, 'DisplayName', 'Noise Detections');
legend_handles = [h_radar, h_noise];
legend_names = {'Radar', 'Noise Detections'};

if ~isempty(baselineObjects)
    % Plot noise detections
    set(h_noise, 'XData', [baselineObjects(baselineOutlierIdx).x], ...
                 'YData', [baselineObjects(baselineOutlierIdx).y], ...
                 'ZData', [baselineObjects(baselineOutlierIdx).z]);
    
    % Cluster plotting
    numTotalClusters = length(baselineClusters);
    clusterColors = lines(numTotalClusters);
    
    for i = 1:numTotalClusters
        clusterPoints = baselineClusters(i).detections;
        markerStyle = 'o';
        if baselineClusters(i).isMoving
            markerStyle = 's'; % Square for moving
        else % isStatic
            markerStyle = 'o'; % Circle for static
        end
        
        h_obj = plot3(ax_baseline, [clusterPoints.x], [clusterPoints.y], [clusterPoints.z], ...
              'Marker', markerStyle, ...
              'MarkerFaceColor', clusterColors(i,:), ...
              'MarkerEdgeColor', 'k', ...
              'MarkerSize', 6, ...
              'DisplayName', ['Object ', num2str(i)]);
        
        legend_handles(end+1) = h_obj;
        legend_names{end+1} = ['Object ', num2str(i), ' (', ifelse(baselineClusters(i).isStatic, 'Static', 'Moving'), ')'];
    end
end
ax_baseline.XLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
ax_baseline.YLim = [0, MAX_RANGE_DISPLAY];
ax_baseline.ZLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
xlabel(ax_baseline, 'X (m)'); ylabel(ax_baseline, 'Y (m) - Forward'); zlabel(ax_baseline, 'Z (m) - Height');
title(ax_baseline, 'Baseline Objects (Static vs Moving)');
legend(ax_baseline, legend_handles, legend_names, 'Location', 'bestoutside');
hold(ax_baseline, 'off');

%% ------------------------------------------------------------------------
% Subplot 2: Interference Objects
ax_interference = subplot(1,2,2); % right subplot
view(ax_interference, 3);
grid(ax_interference, 'on'); hold(ax_interference, 'on');
h_radar_int = plot3(ax_interference, 0, 0, 0, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'blue', 'DisplayName', 'Radar');
h_noise_int = plot3(ax_interference, nan, nan, nan, 'o', 'MarkerFaceColor', NOISE_COLOR, 'MarkerEdgeColor', 'k', 'MarkerSize', 3, 'DisplayName', 'Noise Detections');
legend_handles_int = [h_radar_int, h_noise_int];
legend_names_int = {'Radar', 'Noise Detections'};

if ~isempty(interferenceObjects)
    % Plot noise detections
    set(h_noise_int, 'XData', [interferenceObjects(interferenceOutlierIdx).x], ...
                     'YData', [interferenceObjects(interferenceOutlierIdx).y], ...
                     'ZData', [interferenceObjects(interferenceOutlierIdx).z]);
    
    numTotalClusters_int = length(interferenceClusters);
    clusterColors_int = lines(numTotalClusters_int);
    for i = 1:numTotalClusters_int
        clusterPoints = interferenceClusters(i).detections;
        markerStyle_int = 'o';
        if interferenceClusters(i).isMoving
            markerStyle_int = 's'; % Square for moving
        else % isStatic
            markerStyle_int = 'o'; % Circle for static
        end
        
        h_obj_int = plot3(ax_interference, [clusterPoints.x], [clusterPoints.y], [clusterPoints.z], ...
              'Marker', markerStyle_int, ...
              'MarkerFaceColor', clusterColors_int(i,:), ...
              'MarkerEdgeColor', 'k', ...
              'MarkerSize', 6, ...
              'DisplayName', ['Object ', num2str(i)]);
        
        legend_handles_int(end+1) = h_obj_int;
        legend_names_int{end+1} = ['Object ', num2str(i), ' (', ifelse(interferenceClusters(i).isStatic, 'Static', 'Moving'), ')'];
    end
end
ax_interference.XLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
ax_interference.YLim = [0, MAX_RANGE_DISPLAY];
ax_interference.ZLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
xlabel(ax_interference, 'X (m)'); ylabel(ax_interference, 'Y (m) - Forward'); zlabel(ax_interference, 'Z (m) - Height');
title(ax_interference, 'Interference Objects (Static vs Moving)');
legend(ax_interference, legend_handles_int, legend_names_int, 'Location', 'bestoutside');
hold(ax_interference, 'off');

%% --- ANALYSIS: RANGE-VELOCITY DISTRIBUTION (Side-by-Side) ---
fprintf('\n------------------------------------------------\n');
disp('--- Analysis: Range-Velocity Scatter Plot ---');
figure('Name', 'Range-Velocity Distribution: Baseline vs Interference', 'NumberTitle', 'off', 'Position', [100 100 1000 500]);
% Subplot for Baseline
subplot(1,2,1);
hold on;
if ~isempty(baselineObjects)
    scatter([baselineObjects.range], [baselineObjects.velocity], 15, 'filled');
    line([0 MAX_RANGE_DISPLAY], [STATIC_VELOCITY_THRESHOLD STATIC_VELOCITY_THRESHOLD], 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Moving Threshold');
    line([0 MAX_RANGE_DISPLAY], [-STATIC_VELOCITY_THRESHOLD -STATIC_VELOCITY_THRESHOLD], 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1.5);
end
title('Baseline Range vs. Velocity');
xlabel('Range (m)');
ylabel('Velocity (m/s)');
xlim([0 MAX_RANGE_DISPLAY]);
grid on;
hold off;
% Subplot for Interference
subplot(1,2,2);
hold on;
if ~isempty(interferenceObjects)
    scatter([interferenceObjects.range], [interferenceObjects.velocity], 15, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
    line([0 MAX_RANGE_DISPLAY], [STATIC_VELOCITY_THRESHOLD STATIC_VELOCITY_THRESHOLD], 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Moving Threshold');
    line([0 MAX_RANGE_DISPLAY], [-STATIC_VELOCITY_THRESHOLD -STATIC_VELOCITY_THRESHOLD], 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1.5);
end
title('Interference Range vs. Velocity');
xlabel('Range (m)');
ylabel('Velocity (m/s)');
xlim([0 MAX_RANGE_DISPLAY]);
grid on;
hold off;

%% --- ANALYSIS: RANGE-AZIMUTH POLAR DISTRIBUTION (Side-by-Side) ---
fprintf('\n------------------------------------------------\n');
disp('--- Analysis: Range-Azimuth Polar Distribution ---');
figure('Name', 'Range-Azimuth Polar Distribution: Baseline vs Interference', 'NumberTitle', 'off', 'Position', [100 100 1000 500]);
% Subplot for Baseline Polar
subplot(1,2,1);
if ~isempty(baselineObjects)
    polarplot([baselineObjects.azimuth_deg]*pi/180, [baselineObjects.range], '.', 'MarkerSize', 10, 'Color', [0 0.4470 0.7410]);
end
title('Baseline Range-Azimuth');
rlim([0 MAX_RANGE_DISPLAY]);
thetalim([-90 90]);
grid on;
% Subplot for Interference Polar
subplot(1,2,2);
if ~isempty(interferenceObjects)
    polarplot([interferenceObjects.azimuth_deg]*pi/180, [interferenceObjects.range], '.', 'MarkerSize', 10, 'Color', [0.8500 0.3250 0.0980]);
end
title('Interference Range-Azimuth');
rlim([0 MAX_RANGE_DISPLAY]);
thetalim([-90 90]);
grid on;

%% --- VISUALIZATION: METRIC DISTRIBUTIONS (HISTOGRAMS) ---
fprintf('\n------------------------------------------------\n');
disp('--- Analysis: Distribution of Key Metrics ---');
figure('Name', 'Metric Distributions: Baseline vs Interference', 'NumberTitle', 'off', 'Position', [100 100 1200 800]);
% --- Range Distribution ---
subplot(2,2,1);
hold on;
if ~isempty(baselineObjects)
    histogram([baselineObjects.range], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0 0.4470 0.7410]);
end
if ~isempty(interferenceObjects)
    histogram([interferenceObjects.range], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.7);
end
title('Range Distribution'); xlabel('Range (m)'); ylabel('Probability Density');
legend('Baseline', 'Interference'); grid on; hold off;
% --- Velocity Distribution ---
subplot(2,2,2);
hold on;
if ~isempty(baselineObjects)
    histogram([baselineObjects.velocity], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0 0.4470 0.7410]);
end
if ~isempty(interferenceObjects)
    histogram([interferenceObjects.velocity], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.7);
end
title('Velocity Distribution'); xlabel('Velocity (m/s)'); ylabel('Probability Density');
legend('Baseline', 'Interference'); grid on; hold off;
% --- SNR Distribution ---
subplot(2,2,3);
hold on;
if ~isempty(baselineObjects)
    histogram([baselineObjects.snr_db], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0 0.4470 0.7410]);
end
if ~isempty(interferenceObjects)
    histogram([interferenceObjects.snr_db], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.7);
end
title('SNR Distribution'); xlabel('SNR (dB)'); ylabel('Probability Density');
legend('Baseline', 'Interference'); grid on; hold off;
% --- Noise Distribution ---
subplot(2,2,4);
hold on;
if ~isempty(baselineObjects)
    histogram([baselineObjects.noise_db], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0 0.4470 0.7410]);
end
if ~isempty(interferenceObjects)
    histogram([interferenceObjects.noise_db], 'BinMethod', 'auto', 'Normalization', 'probability', 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.7);
end
title('Noise Distribution'); xlabel('Noise (dB)'); ylabel('Probability Density');
legend('Baseline', 'Interference'); grid on; hold off;
disp(' ');
disp('--- Analysis Complete ---');
disp('Check the generated figures for visualization of the data distributions.');

%% --- HELPER FUNCTIONS ---
% Helper function for clustering and classification
function [clusteredObjects, numStatic, numMoving, outlierIdx] = clusterAndClassify(allDetections, epsilon, minPts, velocityThreshold)
    if isempty(allDetections)
        clusteredObjects = [];
        numStatic = 0; numMoving = 0; 
        outlierIdx = [];
        return;
    end
    
    points = [[allDetections.x]', [allDetections.y]', [allDetections.z]'];
    clusterIdx = customDBSCAN(points, epsilon, minPts);
    
    outlierIdx = (clusterIdx == -1);
    
    numClusters = max(clusterIdx);
    
    % Structure updated to remove 'isInterference' field
    clusteredObjects = struct('detections', {}, 'avgVelocity', {}, 'isStatic', {}, 'isMoving', {});
    
    numStatic = 0;
    numMoving = 0;
    
    for i = 1:numClusters
        clusterDetections = allDetections(clusterIdx == i);
        
        avgVel = mean(abs([clusterDetections.velocity]));
        
        % Simple Static/Moving Classification Logic
        isMoving = avgVel > velocityThreshold;
        isStatic = ~isMoving;
        
        % Structure fields updated to remove interference fields
        clusteredObjects(i).detections = clusterDetections;
        clusteredObjects(i).avgVelocity = avgVel;
        clusteredObjects(i).isStatic = isStatic;
        clusteredObjects(i).isMoving = isMoving;
        
        if isStatic
            numStatic = numStatic + 1;
        else % isMoving
            numMoving = numMoving + 1;
        end
    end
end

% A simple function to replicate the ternary operator for clarity
function result = ifelse(condition, true_result, false_result)
    if condition
        result = true_result;
    else
        result = false_result;
    end
end

%% ========================= CUSTOM DBSCAN FUNCTION ========================
% This function performs DBSCAN clustering without requiring the
% Statistics and Machine Learning Toolbox.
function clusterIdx = customDBSCAN(points, epsilon, minPts)
    numPoints = size(points, 1);
    clusterIdx = zeros(numPoints, 1); % 0 means unvisited
    clusterCounter = 0;
    
    % Manually calculate pairwise Euclidean distances without toolboxes
    Z = zeros(numPoints, numPoints);
    for i = 1:numPoints
        for j = i+1:numPoints
            dist = sqrt(sum((points(i,:) - points(j,:)).^2));
            Z(i,j) = dist;
            Z(j,i) = dist;
        end
    end
    
    for i = 1:numPoints
        if clusterIdx(i) == 0
            % Find neighborhood of current point
            neighbors = find(Z(i, :) <= epsilon);
            
            % If neighborhood is too small, mark as noise
            if length(neighbors) < minPts
                clusterIdx(i) = -1;
            else
                % Expand cluster
                clusterCounter = clusterCounter + 1;
                clusterIdx(i) = clusterCounter;
                
                % Process neighbors
                k = 1;
                while k <= length(neighbors)
                    neighborIdx = neighbors(k);
                    if clusterIdx(neighborIdx) == 0
                        clusterIdx(neighborIdx) = clusterCounter;
                        
                        % Find neighborhood of neighbor
                        newNeighbors = find(Z(neighborIdx, :) <= epsilon);
                        
                        % If neighbor is a core point, add its neighbors to the list
                        if length(newNeighbors) >= minPts
                            neighbors = [neighbors, setdiff(newNeighbors, neighbors)];
                        end
                    elseif clusterIdx(neighborIdx) == -1
                        clusterIdx(neighborIdx) = clusterCounter;
                    end
                    k = k + 1;
                end
            end
        end
    end
end
