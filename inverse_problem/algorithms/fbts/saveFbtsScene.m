function references = saveFbtsScene(root, sceneId, cfg, options)
%saveFbtsScene Save shared system and exact truth as top-level numeric arrays.
[cfg,time,pulse,weight,source] = prepareFbtsSource(cfg);
template = fbtsPortableConfig(cfg);
remove = intersect({'grid','targets','generation'},fieldnames(template));
template = rmfield(template,remove);
template.antennas = rmfield(template.antennas,'doiMask');
template.source = rmfield(template.source,'samples');
gridMeta = cfg.grid;
dense = {'epsr','epsr_bg','cond_e','cond_e_bg','cond_m','cond_m_bg', ...
    'murx','murx_bg','mury','mury_bg','xIndex','yIndex','xPhysical','yPhysical'};
gridMeta = rmfield(gridMeta,intersect(dense,fieldnames(gridMeta)));
coarse = cfg; coarse.grid.epsr = cfg.grid.epsr_bg;
coarse = prepareCoarseSensitivityCfg(coarse,options.sensitivityDownsampleFactor);
coarseTime = (0:coarse.Nt-1)*coarse.dt;
model = struct('schema_version',1,'cfg_template',template,'grid_metadata',gridMeta, ...
    'array_order','x,y; receiver tensors tx,rx,time', ...
    'length_units','m','time_units','s','conductivity_units','S/m', ...
    'permittivity_units','relative, dimensionless','source',source, ...
    'estimated_properties',{{'epsr'}},'prescribed_properties',{{'cond_e','cond_m','murx','mury'}}, ...
    'effective_solver_constants',struct('eps0',8.854187817e-12,'mu0',4*pi*1e-7, ...
    'c0',1/sqrt(8.854187817e-12*4*pi*1e-7)));
model.antennas = template.antennas;
model.antennas.original_ids = 1:cfg.antennas.numAntennas;
model.antennas.positions_m = cfg.grid.originPhysical + (cfg.antennas.pos-1).*[cfg.dx cfg.dy];
model.antennas.nominal_angles_rad = (0:cfg.antennas.numAntennas-1)*2*pi/cfg.antennas.numAntennas;
model.sensitivity = struct('factor',options.sensitivityDownsampleFactor, ...
    'Nx',coarse.Nx,'Ny',coarse.Ny,'dx',coarse.dx,'dy',coarse.dy, ...
    'dt',coarse.dt,'Nt',coarse.Nt,'pml',coarse.pml, ...
    'antenna_positions',coarse.antennas.pos,'time_s',coarseTime, ...
    'source_pulse',cfg.source.func(coarseTime), ...
    'time_weight',interp1(time,weight,coarseTime,'linear',0), ...
    'resampling','fdtdmat.downsampleCfg; see source snapshot', ...
    'time_interpolation','linear, zero outside domain');
system = struct('model',model, ...
    'x_m',cfg.grid.originPhysical(1)+(0:cfg.Nx-1)*cfg.dx, ...
    'y_m',cfg.grid.originPhysical(2)+(0:cfg.Ny-1)*cfg.dy, ...
    'time_s',time,'sourcePulse',pulse,'timeWeight',weight, ...
    'epsr_background',cfg.grid.epsr_bg,'cond_e_background',cfg.grid.cond_e_bg, ...
    'cond_m_background',cfg.grid.cond_m_bg,'murx_background',cfg.grid.murx_bg, ...
    'mury_background',cfg.grid.mury_bg,'epsr_initial',cfg.grid.epsr_bg, ...
    'doi_mask',logical(cfg.antennas.doiMask),'doi_linear_indices',find(cfg.antennas.doiMask));
system.valid_domain_mask = true(cfg.Nx,cfg.Ny);
if isfield(cfg.pml,'condx') && isfield(cfg.pml,'condy')
    system.valid_domain_mask = cfg.pml.condx==0 & cfg.pml.condy==0;
end
systemId = ['system_' fbtsHash(system)];
systemPath = fullfile('systems',[systemId '.mat']);
commitImmutable(fullfile(root,systemPath),system);

targets = cfg.targets;
masks = false(cfg.Nx,cfg.Ny,numel(targets));
labels = zeros(cfg.Nx,cfg.Ny,'uint32');
assembled = cfg.grid.epsr_bg;
for k = 1:numel(targets)
    if isfield(targets(k),'mask') && ~isempty(targets(k).mask)
        mask = logical(targets(k).mask);
    else
        switch string(targets(k).name)
            case "circle"
                r = fdtdgeom.shape_circle(cfg.grid,targets(k).properties.center,targets(k).properties.radius,'index');
            case {"rectangle","square"}
                r = fdtdgeom.shape_rectangle(cfg.grid,targets(k).properties.bounds,'index');
            case {"triangle","hexagon","polygon"}
                r = fdtdgeom.shape_polygon(cfg.grid,targets(k).properties.vertices,'index');
            otherwise
                error('fbts:MissingTargetMask','Supply an actual simulation-grid mask for this shape.');
        end
        mask = r.mask;
    end
    assert(isequal(size(mask),[cfg.Nx cfg.Ny]),'fbts:InvalidTargetMask','Target mask size mismatch.');
    masks(:,:,k) = mask; labels(mask) = k;
    assembled(mask) = targets(k).material.epsr;
    targets(k).object_id = k;
    if ~isfield(targets(k),'physical') || isempty(targets(k).physical)
        targets(k).physical = physicalGeometry(targets(k),cfg);
    end
end
assert(isequal(assembled,cfg.grid.epsr),'fbts:TruthMaskMismatch', ...
    'Target masks/materials do not reproduce the actual FDTD permittivity map. Supply aligned rasterized masks.');
if isfield(targets,'mask'), targets = rmfield(targets,'mask'); end
generation = struct('seed',[],'provenance','imported scene; generation RNG unknown');
if isfield(cfg,'generation'), generation = cfg.generation; end
scene = struct('schema_version',1,'scene_id',sceneId,'system_id',systemId, ...
    'generation',generation,'overlap_policy','construction order; last object sets label/material', ...
    'config_hash',fbtsHash(fbtsPortableConfig(cfg)), ...
    'data_assumptions','synthetic, same forward model, lossless nonmagnetic');
truth = struct('scene',scene,'targets',targets,'epsr_true',cfg.grid.epsr, ...
    'cond_e_true',cfg.grid.cond_e,'cond_m_true',cfg.grid.cond_m, ...
    'murx_true',cfg.grid.murx,'mury_true',cfg.grid.mury, ...
    'target_masks',masks,'target_union_mask',any(masks,3),'label_image',labels);
truthPath = fullfile('scenes',sceneId,'truth.mat');
commitImmutable(fullfile(root,truthPath),truth);
references = struct('system_id',systemId,'system_file',systemPath,'truth_file',truthPath, ...
    'acquisition_file',fullfile('scenes',sceneId,'acquisition.mat'), ...
    'system_hash',fbtsHash(fullfile(root,systemPath),'file'), ...
    'truth_hash',fbtsHash(fullfile(root,truthPath),'file'));
end

function commitImmutable(path,payload)
hash = fbtsHash(payload);
if isfile(path)
    stored = load(path);
    assert(isfield(stored,'content_hash') && strcmp(hash,stored.content_hash) && ...
        strcmp(hash,fbtsHash(rmfield(stored,'content_hash'))), ...
        'fbts:DatasetMismatch','Existing shared dataset differs or is corrupt: %s.',path);
else
    payload.content_hash = hash;
    saveFbtsAtomic(path,payload);
end
end

function physical = physicalGeometry(target,cfg)
p = target.properties; origin = cfg.grid.originPhysical;
physical = struct('coordinate_units','m');
if isfield(p,'center'), physical.center_m = origin+(p.center-1).*[cfg.dx cfg.dy]; end
if isfield(p,'radius'), physical.radius_m = p.radius*cfg.dx; end
if isfield(p,'bounds')
    physical.bounds_m = (p.bounds-1).*[cfg.dx cfg.dx cfg.dy cfg.dy]+origin([1 1 2 2]);
end
if isfield(p,'vertices'), physical.vertices_m = origin+(p.vertices-1).*[cfg.dx cfg.dy]; end
if isfield(p,'orientation_deg'), physical.orientation_deg = p.orientation_deg; end
end
