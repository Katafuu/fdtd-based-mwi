function [scenes, index, timing] = buildFbtsStudyScenes(recordingSeconds, seed)
%buildFbtsStudyScenes Twelve fixed-size, single-target scenes in epsr=45.
if nargin < 1, recordingSeconds = 12e-9; end
if nargin < 2, seed = 42; end
validateattributes(recordingSeconds,{'numeric'},{'scalar','finite','positive'});
stream = RandStream('mt19937ar','Seed',seed);
started = tic; seconds = zeros(12,1);
scenes = cell(1,12); rows = cell(12,7);
shapes = {'circle','triangle','square','hexagon'};
[x,y] = ndgrid(1:400,1:400); centers = [x(:),y(:)];
relative = (centers-[200 200])*1e-3;
for s = 1:4
    shape = shapes{s};
    switch shape
        case 'circle', offsets = zeros(1,2);
        case 'triangle', angle = pi/2+(0:2)'*2*pi/3; offsets = .01*[cos(angle),sin(angle)];
        case 'square', offsets = .01*[-1 -1;1 -1;1 1;-1 1];
        case 'hexagon', angle = (0:5)'*pi/3; offsets = .01*[cos(angle),sin(angle)];
    end
    valid = true(size(x(:)));
    if strcmp(shape,'circle')
        valid = hypot(relative(:,1),relative(:,2))+.01 <= .075-1e-12;
    else
        for v = 1:size(offsets,1)
            valid = valid & hypot(relative(:,1)+offsets(v,1),relative(:,2)+offsets(v,2)) <= .075-1e-12;
        end
    end
    candidates = find(valid); before = stream.State;
    selected = candidates(randperm(stream,numel(candidates),3)); after = stream.State;
    for location = 1:3
        k = (s-1)*3+location; timer = tic;
        center = centers(selected(location),:);
        target = struct('name',shape,'properties',struct(), ...
            'material',struct('epsr',2,'cond_e',0,'cond_m',0,'murx',1,'mury',1));
        switch shape
            case 'circle', target.properties = struct('center',center,'radius',10);
            case 'square', target.properties = struct('bounds',[center(1)+[-10 10],center(2)+[-10 10]]);
            otherwise, target.properties = struct('vertices',center+offsets/1e-3);
        end
        options = struct('seed',seed,'numAntennas',8,'backgroundPermittivity',45, ...
            'doiDiameter',.15,'deltaF',1/recordingSeconds,'targetSpecs',target);
        cfg = buildFbtsConfig(options);
        cfg.study = struct('shape',shape,'location_number',location, ...
            'center_m',(center-1)*1e-3,'target_permittivity',2, ...
            'nominal_size_m',.02,'size_convention','square edge; otherwise circumdiameter', ...
            'recording_seconds_requested',recordingSeconds);
        cfg.generation.placement = struct('seed',seed,'rng_algorithm','mt19937ar', ...
            'before',before,'after',after,'candidate_count',numel(candidates), ...
            'selected_candidate_indices',selected,'sampling','uniform valid grid centers without replacement');
        scenes{k} = cfg; seconds(k) = toc(timer);
        rows(k,:) = {sprintf('scene_%06d',k),shape,location, ...
            cfg.study.center_m(1),cfg.study.center_m(2),.02, ...
            fullfile('scenes',sprintf('scene_%06d',k),'truth.mat')};
    end
end
index = cell2table(rows,'VariableNames',{'scene_id','shape','location_number', ...
    'center_x_m','center_y_m','nominal_size_m','truth_file'});
timing = struct('scene_seconds',seconds,'wall_seconds',toc(started));
end
