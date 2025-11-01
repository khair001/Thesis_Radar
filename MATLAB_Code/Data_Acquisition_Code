% MATLAB mmWave Radar Real-time 3D Observer
% This script connects to an mmWave radar, sends a configuration,
% parses detected object data (TLV Type 1 and Type 7), and
% visualizes it in a real-time 3D plot.
%
% This is designed to be the 'Observer' radar in an interference experiment.
% It now includes functionality to log all detected and filtered objects,
% with guaranteed saving via onCleanup AND periodic saving.
% To run this:
% 2. In MATLAB Command Window, ensure your current folder contains this .m file
%    and your 'my_profile.cfg' file.


clear all; close all; clc; 
% --- CONFIGURATION ---
% IMPORTANT: Update these COM port names and config file 
cfgPortName = 'COM7';  % Configuration port 
dataPortName = 'COM8'; % Data port 
configFile = 'my_profile.cfg'; % Ensure this file is in the same directory
cfgBaudRate = 115200;
dataBaudRate = 921600;

% --- RADAR FRAME HEADER CONSTANTS ---
% Defined in mmw_output.h (MmwDemo_output_message_header_t)
% Magic Word (8 bytes) - represented as a decimal uint8 array in MATLAB
MAGIC_WORD = uint8([2 1 4 3 6 5 8 7]);
HEADER_SIZE = 40; % Bytes

% --- DETECTED OBJECT STRUCTURES (from dpif_pointcloud.h) ---
% DPIF_PointCloudCartesian_t (x, y, z, velocity)
CARTESIAN_OBJ_SIZE = 16; % Bytes
% DPIF_PointCloudSideInfo_t (snr, noise) 
SIDE_INFO_OBJ_SIZE = 4; % Bytes

% --- TLV TYPE DEFINITIONS (from mmw_output.h) ---
MMWDEMO_OUTPUT_MSG_DETECTED_POINTS = 1; % TLV Type for Cartesian points
MMWDEMO_OUTPUT_MSG_DETECTED_POINTS_SIDE_INFO = 7; % TLV Type for SNR/Noise

% --- SERIAL PORT SETUP ---
global cfgPort;
global dataPort;
global allLoggedObjects; 

% Initialize serial port variables as empty. This ensures they always exist
% and prevents 'Reference to cleared variable' errors if port opening fails.
cfgPort = []; 
dataPort = [];

try
    % Configure and open the serial ports
    cfgPort = serialport(cfgPortName, cfgBaudRate);
    dataPort = serialport(dataPortName, dataBaudRate);
    % Set a shorter timeout for the data port for better responsiveness
    dataPort.Timeout = 0.1;
    disp(['Opened serial ports: ', cfgPortName, ' (CFG) and ', dataPortName, ' (DATA)']);
    fprintf('INFO: Data port baud rate is set to: %d\n', dataBaudRate);
catch ME
    disp(['Error opening serial port: ', ME.message]);
    disp('Please check if the ports are correct and not already in use.');
    % The 'cleanupObj' will handle any successfully opened ports upon script exit.
    return; % Exit script if ports cannot be opened
end

% --- SEND CONFIGURATION TO RADAR ---
disp(['Sending configuration from ', configFile, ' to radar...']);
try
    fileID = fopen(configFile, 'r');
    if fileID == -1
        error(['Could not open config file: ', configFile]);
    end
    % Read each line and send as a command
    tline = fgetl(fileID);
    while ischar(tline)
        if ~isempty(tline) && tline(1) ~= '%' % Ignore empty lines and comments
            writeline(cfgPort, tline);
            pause(0.01); % Small delay between commands
        end
        tline = fgetl(fileID);
    end
    fclose(fileID);
    disp('Configuration sent successfully.');
catch ME
    disp(['Error sending configuration: ', ME.message]);
    return;
end

% --- 3D PLOTTING SETUP ---
global hFigure; % Declare figure handle as global
ENABLE_REALTIME_PLOT = true; % Set to false to disable live plotting and save resources

if ENABLE_REALTIME_PLOT
    hFigure = figure('Name', 'mmWave Radar Real-time 3D Observer', 'NumberTitle', 'off', 'Position', [100 100 800 800]);
    ax = axes('Parent', hFigure, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
    view(ax, 3); % Set 3D view
    hold(ax, 'on');
    MAX_RANGE_DISPLAY = 5.0; % meters (Adjust this value as needed)
    % Initialize the 3D scatter plot for objects
    hScatter = plot3(ax, NaN, NaN, NaN, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    % Set plot limits
    ax.XLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
    ax.YLim = [0, MAX_RANGE_DISPLAY]; % Y-axis is forward from radar
    ax.ZLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
    xlabel(ax, 'X (meters)');
    ylabel(ax, 'Y (meters) - Forward');
    zlabel(ax, 'Z (meters) - Height');
    title(ax, ['mmWave Radar Detected Objects (3D View up to ', num2str(MAX_RANGE_DISPLAY), 'm)']);
    % Add a representation of the radar at the origin
    hRadar = plot3(ax, 0, 0, 0, 's', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 10);
    legend(ax, [hScatter, hRadar], {'Detected Objects', 'Radar'}, 'Location', 'northwest');
    % List to hold text annotations for objects (needed to clear them each frame)
    hTextLabels = gobjects(0); % Initialize as empty graphics object array
else
    hFigure = []; % No figure created
    hScatter = [];
    hTextLabels = [];
    disp('[INFO] Real-time plotting is DISABLED.');
end

% --- GLOBAL BUFFER ---
dataBuffer = uint8([]);

% --- Logging Data Setup ---
% Global variable to store all filtered objects for logging
allLoggedObjects = struct('x', {}, 'y', {}, 'z', {}, 'velocity', {}, 'snr_db', {}, 'noise_db', {}, 'range', {}, 'azimuth_deg', {});

% CRUCIAL GLOBAL DECLARATIONS FOR LOGGING:
global saveLogFile; 
global logDirectory; 
global finalLogFileName; 
global periodicLogFileName; 
global runScenario; 

runScenario = 'interference'; % <--- IMPORTANT: Change this to 'baseline' or 'interference' before running for data collection!

% Construct the full path to the 'radar_logs' directory.
logDirectory = fullfile(pwd, 'radar_logs'); 
% The final log file name will depend on the 'runScenario' setting.
finalLogFileName = fullfile(logDirectory, [runScenario, '_data.mat']); 
% Periodic save filename (for in-progress saves)
periodicLogFileName = fullfile(logDirectory, [runScenario, '_periodic_data.mat']);

% Ensure saveLogFile is explicitly logical
saveLogFile = logical(true);

% --- SETUP ONCLEANUP FOR GUARANTEED SERIAL PORT CLOSING ---
cleanupObj = onCleanup(@() cleanupFunction());

% --- PERIODIC SAVING CONFIGURATION ---
saveIntervalFrames = 20; % Save every 20 frames for robustness
frameCount = 0;           % Initialize frame counter

% --- NEW: Maximum Frames to Collect to Prevent Overload ---
MAX_FRAMES_TO_COLLECT = 500; 

% --- MAIN DATA ACQUISITION AND PLOTTING LOOP ---
disp('Starting data acquisition. Press Ctrl+C to stop.');
try
    if saveLogFile && ~exist(logDirectory, 'dir')
        mkdir(logDirectory);
        disp(['[INFO] Created logging directory: ', logDirectory]);
    end

    while frameCount < MAX_FRAMES_TO_COLLECT 
        % Read available bytes from serial port
        bytesAvailable = dataPort.NumBytesAvailable;
        if bytesAvailable > 0
            newData = read(dataPort, bytesAvailable, 'uint8');
            dataBuffer = [dataBuffer, newData];
        else
            % Add a tiny pause to yield CPU and allow MATLAB to process events (like Ctrl+C)
            pause(0.001);
            drawnow limitrate;
        end

        % Attempt to parse a full radar frame
        [parsedObjects, dataBuffer, frameParsed] = parseRadarFrame(dataBuffer, MAGIC_WORD, HEADER_SIZE, ...
                                                                 CARTESIAN_OBJ_SIZE, SIDE_INFO_OBJ_SIZE, ...
                                                                 MMWDEMO_OUTPUT_MSG_DETECTED_POINTS, ...
                                                                 MMWDEMO_OUTPUT_MSG_DETECTED_POINTS_SIDE_INFO);
        
        if frameParsed
            frameCount = frameCount + 1; % Increment frame counter for each parsed frame

            % Filter objects by display range
            filteredObjects = [];
            if ~isempty(parsedObjects)
                for i = 1:length(parsedObjects)
                    obj = parsedObjects(i);
                    % Basic filtering for non-zero/noise objects and within display range
                    if (abs(obj.x) > 1e-6 || abs(obj.y) > 1e-6 || abs(obj.z) > 1e-6 || ...
                        abs(obj.velocity) > 1e-6 || obj.snr_db > -50.0) && ...
                       (obj.range <= MAX_RANGE_DISPLAY)
                        filteredObjects = [filteredObjects; obj]; 
                        
                        % --- NEW: Print a statement with the velocity for each object
                        fprintf('INFO: Detected object with velocity: %.2f m/s\n', obj.velocity);
                        % --- END NEW ---
                    end
                end 
            end
            
            % Log filtered objects for this frame
            if ~isempty(filteredObjects)
                allLoggedObjects = [allLoggedObjects; filteredObjects]; 
                fprintf('DEBUG (Main Loop): Added %d objects to log. Total logged objects now: %d\n', length(filteredObjects), length(allLoggedObjects));
            end

            % --- PERIODIC SAVING ---
            if mod(frameCount, saveIntervalFrames) == 0 && saveLogFile && ~isempty(allLoggedObjects)
                % Create a temporary variable to save, so as not to affect allLoggedObjects during accumulation
                tempAllObjects = allLoggedObjects; 
                save(periodicLogFileName, 'tempAllObjects'); 
                disp(['Periodically saved data to ', periodicLogFileName, ' at frame ', num2str(frameCount)]);
            end

            % Update 3D plot ONLY if ENABLE_REALTIME_PLOT is true
            if ENABLE_REALTIME_PLOT && isvalid(hScatter) && isvalid(hFigure)
                if ~isempty(filteredObjects)
                    x = [filteredObjects.x];
                    y = [filteredObjects.y]; 
                    z = [filteredObjects.z];
                    set(hScatter, 'XData', x, 'YData', y, 'ZData', z);
    
                    % Clear previous text labels
                    if ~isempty(hTextLabels) && all(isvalid(hTextLabels))
                        delete(hTextLabels);
                    end
                    hTextLabels = gobjects(0); % Reset array
    
                    % Add new text labels for range and azimuth
                    for i = 1:length(filteredObjects)
                        obj = filteredObjects(i);
                        % Check if coordinates are finite before adding text to prevent errors
                        if all(isfinite([obj.x, obj.y, obj.z]))
                            labelText = sprintf('R:%.2fm\\nA:%.1f%c', obj.range, obj.azimuth_deg, char(176)); % char(176) is degree symbol, \n for new line in text
                            set(0, 'DefaultTextInterpreter', 'none');
                            hTextLabels(end+1) = text(ax, double(obj.x), double(obj.y), double(obj.z), labelText, 'Color', 'k', 'FontSize', 8, 'Interpreter', 'none'); %#ok<AGROW>
                        else
                            fprintf('WARNING: Skipping text label for object %d due to non-finite coordinates (x=%.2f, y=%.2f, z=%.2f).\n', i, obj.x, obj.y, obj.z);
                        end
                    end
                else
                    % Clear points if no objects detected
                    set(hScatter, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
                    if ~isempty(hTextLabels) && all(isvalid(hTextLabels))
                        delete(hTextLabels);
                    end
                    hTextLabels = gobjects(0);
                end
                drawnow limitrate; % Update plot efficiently
            end
        else
            % No valid frame parsed or no objects detected, clear plot to avoid old data lingering
            if ENABLE_REALTIME_PLOT && isvalid(hScatter) && isvalid(hFigure)
                 set(hScatter, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
            end
            if ENABLE_REALTIME_PLOT && ~isempty(hTextLabels) && all(isvalid(hTextLabels))
                delete(hTextLabels);
            end
            if ENABLE_REALTIME_PLOT
                hTextLabels = gobjects(0);
                drawnow limitrate;
            end
        end
    end
    disp(['[INFO] Maximum frame limit (', num2str(MAX_FRAMES_TO_COLLECT), ') reached. Stopping data acquisition.']);

    % --- NEW: FORCE FINAL SAVE IMMEDIATELY AFTER LOOP COMPLETION ---
    if saveLogFile && ~isempty(allLoggedObjects)
        disp('INFO (Main Script): Maximum frames reached. Performing final data save...');
        try
            % Ensure log directory exists before saving
            if ~exist(logDirectory, 'dir')
                mkdir(logDirectory);
                disp(['[INFO] Created logging directory: ', logDirectory, ' during final save.']);
            end
            % Saving the main data variable to the final file
            save(finalLogFileName, 'allLoggedObjects');
            
            disp(['INFO (Main Script): Final data saved to ', finalLogFileName]);
            
            % FIX: Explicitly clear the large data variable after saving to free up memory
            clear allLoggedObjects; 
            fprintf('INFO (Main Script): allLoggedObjects cleared from memory.\n');
        catch saveME % Catch any errors during the final save operation itself
            disp(['ERROR (Main Script): Failed to save final data: ', saveME.message]);
            disp(getReport(saveME)); % Display full error report for save error
        end
    else
        disp('INFO (Main Script): No final data to save as allLoggedObjects is empty or saveLogFile is false.');
    end
    % --- END NEW: FORCE FINAL SAVE ---

catch ME
    % This catch block handles errors during the main data acquisition loop.
    % The onCleanup function will still be called regardless of how the loop exits.
    if strcmp(ME.identifier, 'MATLAB:serialport:read:timeout')
        
        disp('Serial port read timeout.');
    else
        disp(['Error during data acquisition: ', ME.message]);
        disp(getReport(ME)); % Display full error report
    end
end

% --- RADAR FRAME PARSING FUNCTION ---
% This function parses raw byte data from the mmWave radar's data port
% to extract detected object information.
function [parsedObjects, remainingBuffer, frameParsed] = parseRadarFrame(inputBuffer, MAGIC_WORD, HEADER_SIZE, ...
                                                                      CARTESIAN_OBJ_SIZE, SIDE_INFO_OBJ_SIZE, ...
                                                                      MMWDEMO_OUTPUT_MSG_DETECTED_POINTS, ...
                                                                      MMWDEMO_OUTPUT_MSG_DETECTED_POINTS_SIDE_INFO)
    
    parsedObjects = [];          % Initialize empty output for detected objects
    remainingBuffer = inputBuffer; % Assume no frame parsed, return original buffer
    frameParsed = false;         % Flag indicating if a complete frame was successfully parsed
    
    % 1. Search for the Magic Word (frame synchronization)
    magicIdx = strfind(inputBuffer, MAGIC_WORD);
    if isempty(magicIdx)
        return; 
    end
    magicIdx = magicIdx(1); % Take the first occurrence
    
    % 2. Discard leading garbage bytes before the first magic word.
    if magicIdx > 1
        inputBuffer = inputBuffer(magicIdx:end);
    end
    
    % Re-check if enough data for the full header after alignment
    if length(inputBuffer) < HEADER_SIZE
        return;
    end

    % 3. Parse Main Header (40 bytes)
    headerBytes = inputBuffer(1:HEADER_SIZE);
    
    % Use typecast to convert byte arrays to numerical values based on the header structure (little-endian).
    totalPacketLen = typecast(headerBytes(13:16), 'uint32');
    frameNumber = typecast(headerBytes(21:24), 'uint32');
    numDetectedObj = typecast(headerBytes(29:32), 'uint32');
    numTLVs = typecast(headerBytes(33:36), 'uint32');
    
    fprintf('DEBUG (Header): Frame #%d, Total Packet Length: %d, Detected Objects: %d, TLVs: %d\n', frameNumber, totalPacketLen, numDetectedObj, numTLVs);

    % 4. Check for invalid or suspicious header values (robustness check)
    MAX_EXPECTED_PACKET_LEN = 20000; % Increased threshold to allow for more points
    MAX_EXPECTED_OBJECTS = 500; % A reasonable upper limit
    if totalPacketLen > MAX_EXPECTED_PACKET_LEN || totalPacketLen < HEADER_SIZE || numDetectedObj > MAX_EXPECTED_OBJECTS
        fprintf('WARNING (Parse): Corrupted header detected (PacketLen=%d, NumObj=%d). Discarding buffer.\n', totalPacketLen, numDetectedObj);
        remainingBuffer = uint8([]);
        return;
    end

    % 5. Check if enough data for the entire packet (header + all TLVs)
    if length(inputBuffer) < totalPacketLen
        return;
    end
    
    % Extract the complete radar frame packet.
    packet = inputBuffer(1:totalPacketLen);
    
    % Initialize containers for parsed TLV data.
    cartesianPoints = [];
    sideInfoPoints = [];
    
    % 6. Iterate and Parse TLVs (Type-Length-Value blocks)
    currentTlvOffset = HEADER_SIZE;
    for tlvIdx = 0:numTLVs-1
        if length(packet) < currentTlvOffset + 8
            fprintf('WARNING (Parse): Insufficient bytes for TLV %d header. Breaking.\n', tlvIdx);
            break; 
        end
        % Parse TLV header (Type and Length)
        tlvHeader = packet(currentTlvOffset + 1 : currentTlvOffset + 8);
        tlvType = typecast(tlvHeader(1:4), 'uint32');
        tlvLength = typecast(tlvHeader(5:8), 'uint32');
        
        fprintf('DEBUG (TLV): Found TLV %d. Type: %d, Length: %d\n', tlvIdx, tlvType, tlvLength);

        % 7. Validate TLV Length and ensure payload is within packet bounds
        MAX_EXPECTED_TLV_LENGTH = 10000;
        if tlvLength < 8 || tlvLength > MAX_EXPECTED_TLV_LENGTH || (currentTlvOffset + tlvLength > length(packet))
            fprintf('WARNING (Parse): TLV %d (Type: %d) reports invalid length %d or exceeds packet bounds. Skipping this TLV.\n', tlvIdx, tlvType, tlvLength);
            currentTlvOffset = currentTlvOffset + 8; % Move to the next TLV header
            continue;
        end
        
        if length(packet) < currentTlvOffset + tlvLength
             fprintf('WARNING: TLV %d (Type: %d) payload truncated. Skipping.\n', tlvIdx, tlvType);
             currentTlvOffset = currentTlvOffset + tlvLength;
             continue;
        end

        tlvPayloadData = packet(currentTlvOffset + 8 + 1 : currentTlvOffset + tlvLength);
        
        % 8. Parse TLV Payload based on Type
        if tlvType == MMWDEMO_OUTPUT_MSG_DETECTED_POINTS
            numObjectsInPayload = floor(length(tlvPayloadData) / CARTESIAN_OBJ_SIZE);
            if numObjectsInPayload == 0 && length(tlvPayloadData) > 0
                fprintf('WARNING: TLV Type 1 payload length (%d) too small for one object. Skipping.\n', length(tlvPayloadData));
                currentTlvOffset = currentTlvOffset + tlvLength;
                continue;
            end
            
            actualNumObjectsToParse = min(numDetectedObj, numObjectsInPayload);
            cartesianPoints = struct('x', cell(1, actualNumObjectsToParse), 'y', cell(1, actualNumObjectsToParse), 'z', cell(1, actualNumObjectsToParse), 'velocity', cell(1, actualNumObjectsToParse));
            
            for i = 0:actualNumObjectsToParse-1
                if length(tlvPayloadData) < (i+1)*CARTESIAN_OBJ_SIZE
                    fprintf('WARNING: TLV Type 1 payload truncated for object %d. Stopping parsing.\n', i+1);
                    break;
                end
                objData = typecast(tlvPayloadData(i*CARTESIAN_OBJ_SIZE + 1 : (i+1)*CARTESIAN_OBJ_SIZE), 'single');
                cartesianPoints(i+1).x = objData(1);
                cartesianPoints(i+1).y = objData(2);
                cartesianPoints(i+1).z = objData(3);
                cartesianPoints(i+1).velocity = objData(4);
            end
        elseif tlvType == MMWDEMO_OUTPUT_MSG_DETECTED_POINTS_SIDE_INFO
            numSideInfoInPayload = floor(length(tlvPayloadData) / SIDE_INFO_OBJ_SIZE);
            if numSideInfoInPayload == 0 && length(tlvPayloadData) > 0
                fprintf('WARNING: TLV Type 7 payload length (%d) too small for one object. Skipping.\n', length(tlvPayloadData));
                currentTlvOffset = currentTlvOffset + tlvLength;
                continue;
            end
            
            actualNumSideInfoToParse = min(numDetectedObj, numSideInfoInPayload);
            sideInfoPoints = struct('snr_db', cell(1, actualNumSideInfoToParse), 'noise_db', cell(1, actualNumSideInfoToParse));
            
            for i = 0:actualNumSideInfoToParse-1
                if length(tlvPayloadData) < (i+1)*SIDE_INFO_OBJ_SIZE
                    fprintf('WARNING: TLV Type 7 payload truncated for object %d. Stopping parsing.\n', i+1);
                    break;
                end
                infoData = typecast(tlvPayloadData(i*SIDE_INFO_OBJ_SIZE + 1 : (i+1)*SIDE_INFO_OBJ_SIZE), 'int16');
                sideInfoPoints(i+1).snr_db = double(infoData(1)) / 10.0;
                sideInfoPoints(i+1).noise_db = double(infoData(2)) / 10.0;
            end
        else
            fprintf('WARNING: Encountered unknown TLV type: %d. Skipping this TLV.\n', tlvType);
        end
        currentTlvOffset = currentTlvOffset + tlvLength;
    end
    
    % 9. Combine Cartesian points and Side Info into final parsedObjects
    numValidDetectionsForCombining = min(length(cartesianPoints), length(sideInfoPoints));
    parsedObjects = struct('x', {}, 'y', {}, 'z', {}, 'velocity', {}, 'snr_db', {}, 'noise_db', {}, 'range', {}, 'azimuth_deg', {});
    
    for i = 1:numValidDetectionsForCombining
        obj = struct();
        obj.x = cartesianPoints(i).x;
        obj.y = cartesianPoints(i).y;
        obj.z = cartesianPoints(i).z;
        obj.velocity = cartesianPoints(i).velocity;
        obj.snr_db = sideInfoPoints(i).snr_db;
        obj.noise_db = sideInfoPoints(i).noise_db;
        
        % Calculate derived values (range and azimuth angle)
        obj.range = sqrt(obj.x^2 + obj.y^2 + obj.z^2);
        obj.azimuth_deg = rad2deg(atan2(obj.x, obj.y));
        parsedObjects = [parsedObjects; obj]; 
    end
    
    % 10. Update remaining buffer and set frameParsed flag
    remainingBuffer = inputBuffer(totalPacketLen + 1 : end);
    frameParsed = true;
end

% --- CLEANUP FUNCTION (Called automatically by onCleanup) ---
function cleanupFunction()
    global cfgPort;
    global dataPort;
    global saveLogFile;
    global logDirectory;
    global finalLogFileName;
    global periodicLogFileName;
    global runScenario;
    global hFigure;
    
    disp('Stopping acquisition and closing serial ports (from cleanupFunction).');
    
    if ~isempty(cfgPort) && isvalid(cfgPort) && cfgPort.isopen
        try
            writeline(cfgPort, 'sensorStop');
            pause(0.1);
        catch ME_stop
            fprintf('WARNING: Could not send sensorStop command to cfgPort: %s\n', ME_stop.message);
        end
        fclose(cfgPort); 
        delete(cfgPort);
        disp('Closed CFG Port.');
    end
    if ~isempty(dataPort) && isvalid(dataPort) && dataPort.isopen
        fclose(dataPort); 
        delete(dataPort);
        disp('Closed DATA Port.');
    end
    
    if ~isempty(hFigure) && isvalid(hFigure)
        close(hFigure);
        disp('Closed plotting figure.');
    end
    
    if exist(periodicLogFileName, 'file')
        try
            delete(periodicLogFileName);
            disp(['INFO (Cleanup): Deleted temporary periodic log file: ', periodicLogFileName]);
        catch ME_del
            fprintf('WARNING (Cleanup): Could not delete periodic log file: %s\n', ME_del.message);
        end
    end
    
    disp('Cleanup complete.');
end
