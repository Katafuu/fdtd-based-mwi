function save_results(resultFile, cfg, results)
%save_results Save one FBTS run's configuration and results.
save(char(resultFile), 'cfg', 'results', '-v7.3');
end
