function [summary,root,manifest] = runFbtsStudy(directory,includeCuda)
%runFbtsStudy Calibrate, validate, benchmark locally, then run/resume 348 cases.
% No study is launched when a prerequisite fails. Use a new directory after
% any source/build change. GPU runs require the actual NVIDIA computer.
if nargin < 2, includeCuda=false; end
hardware=fbtsHardwareInfo();
if includeCuda
    assert(hardware.gpu_available && hardware.cuda_mex_available,'fbts:CudaUnavailable', ...
        'Run on the NVIDIA PC after building fdtd_cuda. No CPU fallback is performed.');
end
prepared=prepareFbtsStudy(directory);
benchmarkDirectory=fullfile(directory,'benchmarks');
benchmarkFile=fullfile(benchmarkDirectory,'benchmark_report.mat');
if isfile(benchmarkFile)
    saved=load(benchmarkFile,'report'); report=saved.report;
    assert(strcmp(report.code_hash,prepared.code_hash) && ...
        strcmp(report.hardware.machine_id,hardware.machine_id), ...
        'fbts:BenchmarkMismatch','Benchmark belongs to different source/build or computer.');
    if ~report.passed || (includeCuda && isempty(report.cuda))
        report=benchmarkFbtsStudy(prepared,benchmarkDirectory,includeCuda);
    end
else
    report=benchmarkFbtsStudy(prepared,benchmarkDirectory,includeCuda);
end
assert(report.passed && report.cpu_validation.passed,'fbts:StudyNotValidated','CPU validation is incomplete.');
chosen=report.recommended_options;
if ~includeCuda
    chosen=struct('solverBackend','cpu','solverThreads',report.cpu.threads,'maxWorkers',report.cpu.workers);
end
if strcmp(chosen.solverBackend,'cuda')
    assert(includeCuda && report.cuda_validation.passed,'fbts:StudyNotValidated','CUDA validation is incomplete.');
end
% The existing suite covers recovery, immutable scans and process-worker parity.
fbtsFolder=fileparts(mfilename('fullpath'));
repository=fileparts(fileparts(fileparts(fbtsFolder)));
tests=[runtests(fullfile(fbtsFolder,'tests')),runtests(fullfile(repository,'forward_solver','tests'))];
saveFbtsAtomic(fullfile(benchmarkDirectory,'regression_tests.mat'),struct('tests',tests));
assert(~any([tests.Failed]),'fbts:StudyNotValidated','Regression tests failed.');
if chosen.maxWorkers>1
    workers=tests(contains({tests.Name},'processWorkers'));
    assert(~isempty(workers) && all([workers.Passed]),'fbts:StudyNotValidated','Process-worker validation is incomplete.');
end
options=struct('numIterations',15,'sensitivityDownsampleFactor',prepared.sensitivity.factor, ...
    'solverBackend',chosen.solverBackend,'solverThreads',chosen.solverThreads, ...
    'outputDirectory',fullfile(directory,'dataset'));
batch=struct('symmetry','dihedral','seed',123,'maxWorkers',chosen.maxWorkers, ...
    'resume',isfile(fullfile(options.outputDirectory,'manifest.mat')));
if ~batch.resume
    initial=batch; initial.maxCases=0;
    runFbtsBatch(prepared.scenes,options,initial);
    batch.resume=true;
end
destination=fullfile(options.outputDirectory,'benchmarks');
if ~isfolder(destination), mkdir(destination); end
files=dir(benchmarkDirectory);
for k=1:numel(files)
    if ismember(files(k).name,{'.','..'}), continue; end
    copyfile(fullfile(benchmarkDirectory,files(k).name),fullfile(destination,files(k).name));
end
copyfile(fullfile(directory,'calibration'),fullfile(destination,'calibration'));
copyfile(fullfile(directory,'prepared.mat'),fullfile(destination,'prepared.mat'));
[summary,root,manifest]=runFbtsBatch(prepared.scenes,options,batch);
assert(height(summary)==348,'fbts:StudyCatalogMismatch','Expected 348 cases.');
fprintf('Study: %d successful, %d failed, %d pending.\n',sum(summary.status=="success"), ...
    sum(summary.status=="failed"),sum(~ismember(summary.status,["success","failed"])));
end
