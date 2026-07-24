%% 2D TMz Field Visualization
% This script visualizes Ez snapshots produced either by run_fdtd_mex_case.m
% or run_fdtd_file_case.m. It prefers workspace data when available and falls
% back to reading results/ez_field.csv from forward_solver.

%% User Options
fieldSource = "auto";      % "auto", "workspace", or "csv"
frameStride = 1;            % animate every Nth stored frame
frameDelay = 0.02;          % pause between frames during interactive playback
colorMode = "global";      % "global", "frame", or "fixed"
fixedColorLimit = 0.2;      % used only when colorMode == "fixed"
colormapName = turbo;
showSlice = true;           % overlay and plot one 1D slice through the 2D frame
sliceDirection = "y";       % "y" gives Ez(x, y0), "x" gives Ez(x0, y)
sliceIndex = [];            % [] selects the center index
usePhysicalAxes = true;     % use dx/dy when available

%% Video Export Options
exportVideo = false;        % set true to write an MP4 from the final section
videoFile = "ez_field_animation.mp4";
videoFrameRate = 24;
videoQuality = 95;

%% Load 2D Field Stack
[EzStack, meta] = loadEzVisualizationData(fieldSource);
[nx, ny, nt] = size(EzStack);

if isempty(sliceIndex)
    if sliceDirection == "y"
        sliceIndex = round(ny/2);
    else
        sliceIndex = round(nx/2);
    end
end

if usePhysicalAxes && isfield(meta, "dx") && isfield(meta, "dy") && ...
        isfinite(meta.dx) && isfinite(meta.dy)
    xAxis = (0:nx-1) * meta.dx;
    yAxis = (0:ny-1) * meta.dy;
    xLabelText = "x [m]";
    yLabelText = "y [m]";
else
    xAxis = 1:nx;
    yAxis = 1:ny;
    xLabelText = "x index";
    yLabelText = "y index";
end

frameList = 1:frameStride:nt;
if isempty(frameList)
    frameList = 1;
end

switch colorMode
    case "global"
        maxAbs = max(abs(EzStack), [], "all");
        if maxAbs == 0
            maxAbs = 1;
        end
        colorLimit = [-maxAbs maxAbs];
    case "fixed"
        colorLimit = [-fixedColorLimit fixedColorLimit];
    case "frame"
        colorLimit = [];
    otherwise
        error('Unknown colorMode: %s', colorMode);
end

%% Animate 2D Field
fig = figure("Name", "2D Ez field animation", "Color", "w");
if showSlice
    layout = tiledlayout(fig, 2, 1, "TileSpacing", "compact", "Padding", "compact");
    imageAxes = nexttile(layout, 1);
else
    imageAxes = axes(fig);
end

firstFrame = EzStack(:, :, frameList(1));
imageHandle = imagesc(imageAxes, xAxis, yAxis, firstFrame.');
axis(imageAxes, "image");
set(imageAxes, "YDir", "normal");
colormap(imageAxes, colormapName);
colorbar(imageAxes);
xlabel(imageAxes, xLabelText);
ylabel(imageAxes, yLabelText);
title(imageAxes, frameTitle(frameList(1), nt, meta));

if ~isempty(colorLimit)
    clim(imageAxes, colorLimit);
else
    setFrameColorLimit(imageAxes, firstFrame);
end

hold(imageAxes, "on");
if showSlice
    if sliceDirection == "y"
        sliceCoord = yAxis(sliceIndex);
        sliceOverlay = yline(imageAxes, sliceCoord, "w--", "LineWidth", 1.2);
        sliceX = xAxis;
        sliceY = firstFrame(:, sliceIndex);
        sliceXLabel = xLabelText;
        sliceTitle = sprintf("Ez(x, y = %g)", sliceCoord);
    else
        sliceCoord = xAxis(sliceIndex);
        sliceOverlay = xline(imageAxes, sliceCoord, "w--", "LineWidth", 1.2);
        sliceX = yAxis;
        sliceY = firstFrame(sliceIndex, :);
        sliceXLabel = yLabelText;
        sliceTitle = sprintf("Ez(x = %g, y)", sliceCoord);
    end

    sliceAxes = nexttile(layout, 2);
    sliceHandle = plot(sliceAxes, sliceX, sliceY, "LineWidth", 1.4);
    set(sliceAxes, "XGrid", "on", "YGrid", "on");
    xlabel(sliceAxes, sliceXLabel);
    ylabel(sliceAxes, "E_z");
    title(sliceAxes, sliceTitle);
else
    sliceAxes = matlab.graphics.axis.Axes.empty;
    sliceHandle = matlab.graphics.chart.primitive.Line.empty;
    sliceOverlay = matlab.graphics.GraphicsPlaceholder.empty;
end
hold(imageAxes, "off");

viewState = struct();
viewState.imageAxes = imageAxes;
viewState.imageHandle = imageHandle;
viewState.sliceAxes = sliceAxes;
viewState.sliceHandle = sliceHandle;
viewState.sliceOverlay = sliceOverlay;
viewState.xAxis = xAxis;
viewState.yAxis = yAxis;
viewState.showSlice = showSlice;
viewState.sliceDirection = sliceDirection;
viewState.sliceIndex = sliceIndex;
viewState.colorMode = colorMode;
viewState.colorLimit = colorLimit;
viewState.nt = nt;
viewState.meta = meta;

for frameNumber = frameList
    updateEzFrame(EzStack, frameNumber, viewState);
    drawnow;
    pause(frameDelay);
end

%% Export Animation To Video
% Set exportVideo = true in the options section, then run this section to
% write the same 2D animation to an MP4 file. The export reuses the current
% figure and redraws the selected frames at videoFrameRate.

if exportVideo
    videoPath = string(videoFile);
    if ~isfolder(fileparts(videoPath)) && strlength(fileparts(videoPath)) > 0
        mkdir(fileparts(videoPath));
    end

    writer = VideoWriter(videoPath, "MPEG-4");
    writer.FrameRate = videoFrameRate;
    writer.Quality = videoQuality;
    open(writer);

    for frameNumber = frameList
        updateEzFrame(EzStack, frameNumber, viewState);
        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    disp("Exported animation video: " + videoPath);
end

%% Local Functions
function [EzStack, meta] = loadEzVisualizationData(fieldSource)
    meta = struct();
    meta.source = "unknown";
    meta.dx = NaN;
    meta.dy = NaN;
    meta.dt = NaN;

    if evalin("base", "exist('dx','var')")
        meta.dx = evalin("base", "dx");
    end
    if evalin("base", "exist('dy','var')")
        meta.dy = evalin("base", "dy");
    end
    if evalin("base", "exist('dt','var')")
        meta.dt = evalin("base", "dt");
    end

    useWorkspace = fieldSource == "workspace" || fieldSource == "auto";
    if useWorkspace
        if evalin("base", "exist('Ez','var')")
            Ez = evalin("base", "Ez");
            EzStack = normalizeEzStack(Ez);
            meta.source = "workspace Ez";
            return;
        end

        if evalin("base", "exist('fdtdMexResult','var')")
            fdtdMexResult = evalin("base", "fdtdMexResult");
            if isfield(fdtdMexResult, "Ez") && ~isempty(fdtdMexResult.Ez)
                EzStack = normalizeEzStack(fdtdMexResult.Ez);
                meta.source = "workspace fdtdMexResult.Ez";
                if isfield(fdtdMexResult, "dx")
                    meta.dx = fdtdMexResult.dx;
                end
                if isfield(fdtdMexResult, "dy")
                    meta.dy = fdtdMexResult.dy;
                end
                if isfield(fdtdMexResult, "dt")
                    meta.dt = fdtdMexResult.dt;
                end
                return;
            end
        end

        if evalin("base", "exist('snapshots','var')")
            snapshots = evalin("base", "snapshots");
            [nx, ny] = workspaceGridSize(size(snapshots, 2));
            EzStack = snapshotsToStack(snapshots, nx, ny);
            meta.source = "workspace snapshots";
            return;
        end
    end

    if fieldSource == "workspace"
        error('No workspace Ez, fdtdMexResult.Ez, or snapshots variable was found.');
    end

    [csvPath, nx, ny] = findEzCsv();
    snapshots = readmatrix(csvPath, "NumHeaderLines", 1);
    EzStack = snapshotsToStack(snapshots, nx, ny);
    meta.source = "CSV " + string(csvPath);
end

function EzStack = normalizeEzStack(Ez)
    Ez = double(Ez);
    if ismatrix(Ez)
        EzStack = reshape(Ez, size(Ez, 1), size(Ez, 2), 1);
    elseif ndims(Ez) == 3
        EzStack = Ez;
    else
        error('Ez must be a 2D final field or a 3D frame stack.');
    end
end


function EzStack = snapshotsToStack(snapshots, nx, ny)
    snapshots = double(snapshots);
    nt = size(snapshots, 1);
    EzStack = zeros(nx, ny, nt);
    for frame = 1:nt
        EzStack(:, :, frame) = reshape(snapshots(frame, :), [ny, nx]).';
    end
end

function [nx, ny] = workspaceGridSize(numColumns)
    if evalin("base", "exist('gridSize','var')")
        gridSize = evalin("base", "gridSize");
        nx = gridSize(1);
        ny = gridSize(2);
    elseif evalin("base", "exist('Nx','var') && exist('Ny','var')")
        nx = evalin("base", "Nx");
        ny = evalin("base", "Ny");
    else
        n = sqrt(numColumns);
        assert(abs(n - round(n)) < eps(n), ...
            'Could not infer nx and ny from workspace. Define gridSize or Nx/Ny.');
        nx = round(n);
        ny = round(n);
    end
    assert(nx * ny == numColumns, 'Snapshot width does not match nx*ny.');
end

function [csvPath, nx, ny] = findEzCsv()
    candidates = strings(0);
    if evalin("base", "exist('outputPath','var')")
        candidates(end+1) = string(evalin("base", "outputPath"));
    end
    if evalin("base", "exist('repoRoot','var')")
        repoRoot = string(evalin("base", "repoRoot"));
        candidates(end+1) = fullfile(repoRoot, "results", "ez_field.csv");
    end
    candidates(end+1) = fullfile(pwd, "results", "ez_field.csv");
    candidates(end+1) = fullfile(pwd, "forward_solver", "results", "ez_field.csv");
    scriptDir = fileparts(mfilename('fullpath'));
    solverRoot = fileparts(scriptDir);
    candidates(end+1) = fullfile(solverRoot, "results", "ez_field.csv");

    csvPath = "";
    for candidate = candidates
        if isfile(candidate)
            csvPath = candidate;
            break;
        end
    end
    if csvPath == ""
        error('Could not find ez_field.csv. Run run_fdtd_file_case.m or set outputPath.');
    end

    fid = fopen(csvPath, "r");
    assert(fid > 0, 'Could not open %s.', csvPath);
    cleanup = onCleanup(@() fclose(fid));
    header = fgetl(fid);
    dims = str2double(split(string(header), ","));
    assert(numel(dims) >= 2 && all(isfinite(dims(1:2))), ...
        'The first line of %s must contain nx,ny.', csvPath);
    nx = dims(1);
    ny = dims(2);
end

function updateEzFrame(EzStack, frameNumber, viewState)
    frame = EzStack(:, :, frameNumber);
    viewState.imageHandle.CData = frame.';
    title(viewState.imageAxes, frameTitle(frameNumber, viewState.nt, viewState.meta));

    if viewState.colorMode == "frame"
        setFrameColorLimit(viewState.imageAxes, frame);
    elseif ~isempty(viewState.colorLimit)
        clim(viewState.imageAxes, viewState.colorLimit);
    end

    if viewState.showSlice
        if viewState.sliceDirection == "y"
            viewState.sliceHandle.YData = frame(:, viewState.sliceIndex);
        else
            viewState.sliceHandle.YData = frame(viewState.sliceIndex, :);
        end
        if viewState.colorMode == "frame"
            maxAbsSlice = max(abs(viewState.sliceHandle.YData), [], "all");
            if maxAbsSlice == 0
                maxAbsSlice = 1;
            end
            ylim(viewState.sliceAxes, [-maxAbsSlice maxAbsSlice]);
        end
    end
end

function setFrameColorLimit(ax, frame)
    maxAbsFrame = max(abs(frame), [], "all");
    if maxAbsFrame == 0
        maxAbsFrame = 1;
    end
    clim(ax, [-maxAbsFrame maxAbsFrame]);
end

function txt = frameTitle(frameNumber, nt, meta)
    if isfield(meta, "dt") && isfinite(meta.dt)
        txt = sprintf('2D E_z field, frame %d / %d, t = %.4g s (%s)', ...
            frameNumber, nt, (frameNumber - 1) * meta.dt, meta.source);
    else
        txt = sprintf('2D E_z field, frame %d / %d (%s)', ...
            frameNumber, nt, meta.source);
    end
end