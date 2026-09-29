function cfg = fbtsPortableConfig(cfg)
%fbtsPortableConfig Remove the reconstructed pulse handle, keep numeric inputs.
if isfield(cfg, 'source') && isfield(cfg.source, 'func')
    cfg.source = rmfield(cfg.source, 'func');
end
end
