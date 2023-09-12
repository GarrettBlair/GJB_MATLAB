% analysis to do:
% correlation - xcorr - grainger causality - canonical corr analysis
    % between ACC and HPC ipos
% behav var to look into - total path length, 1st/2nd half in water sess,
% inter-entrance times
clear
analysis_version = 'v1.43';
diarydir  = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\diaries\';
diaryname = sprintf('%stemp_reprocessing_%s', diarydir, analysis_version);
diary(diaryname)
params = [];
% APA_rat_imaging_params_current;
APA_PKCZmouse_imaging_params_current;

fprintf('rand shuff pcells = %d\n',params.num_random_shuffle_pcell)
fprintf('rand shuff decode = %d\n',params.num_random_shuffle_decode)

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
fprintf('%s\n', datetime())
%
animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
fprintf('\t\tanalysis ver = %s\n', analysis_version)
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
experiment_folder = [];
animals = {'mHPC23454', 'mHPC23459'};% animals = {'Acc20832', 'Acc19947'};
experiment_folder = {'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\PKCZ_imaging\'};
% experiment_folder{1} = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA/';
% experiment_folder{1} = 'D:/GarrettBlair/APA/';
% animals = {'Acc20832', 'Acc19947', 'Hipp18240'};
numAnimals = length(animals);
% experiment_folder{2} = 'D:/APA recordings/';

% DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

rerun_behav         = false;
rerun_place         = false;
resave_proccessed   = false;
resave_contours     = false;
rerun_processed     = true;
fit_contours_fullFOV= false;

%%
PRINT_ONLY = false;
this_ver = str2double(analysis_version(2:end));
skipfiles = []; updated_files = [];
timer_all = tic;
if rerun_processed == true
    for animalLoop = 1:numAnimals
        for expLoop = 1:length(experiment_folder)
            animal_name = animals{animalLoop};
            processed_dir = sprintf('%s%s/processed_files/', experiment_folder{expLoop}, animal_name);
            fprintf('\n\nREPROCESS files for %s in folder  \n%s...\n', animal_name, processed_dir)
%             temp = dir([processed_dir '*miniscope*']);
            temp = dir([processed_dir '*@placecells*']);
            nfiles = length(temp);
            for sessionLoop = 1:nfiles
                fname = temp(sessionLoop).name;
                processedFile   = [temp(sessionLoop).folder '/' fname];
                prev_version = load(processedFile, 'analysis_version');
                %
                prev_ver = str2double(prev_version.analysis_version(2:end));
                if prev_ver < this_ver
                    timer_session = tic;
                    if ~PRINT_ONLY
                        fprintf('~~~PROCESSED file rerunning for: %s...\n', fname)
                        fprintf('\t\tanalysis_version : %s\n', analysis_version)
                        fprintf('\t\tstart time: %s  \n', datetime())
                        clearvars tempms ms
                        tempms = load(processedFile, 'ms');
                        
%                         ms.spks = normalize_rows(ms.neuron.S_mat);
%                         ms.spks(isnan(ms.spks)) = 0;
                        [ms] = APA_generate_placemaps(tempms.ms, params);
                        
                        if ~isfield(tempms.ms.room, 'ensemble_pos_info')
                            fprintf('\t\trunning ipos...\n')
                            [ms.room.momentary_pos_info,  rtemp, ms.room.ensemble_pos_info]  = Fenton_ipos(ms, params.ipos_int_time, 'room', params);
                            [ms.arena.momentary_pos_info, atemp, ms.arena.ensemble_pos_info] = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
                        else
                            fprintf('\t\ttransferring ipos...\n')
                            ms.room.momentary_pos_info  = tempms.ms.room.momentary_pos_info;
                            ms.room.ensemble_pos_info   = tempms.ms.room.ensemble_pos_info;
                            ms.arena.momentary_pos_info = tempms.ms.arena.momentary_pos_info;
                            ms.arena.ensemble_pos_info  = tempms.ms.arena.ensemble_pos_info;
                        end
                        if ~isfield(tempms.ms.room, 'svm_decoding')
                            fprintf('\t\trunning svm decode...\n')
                            [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                        else
                            fprintf('\t\ttransferring svm decode...\n')
                            ms.room.svm_decoding        = tempms.ms.room.svm_decoding;
                            ms.arena.svm_decoding       = tempms.ms.arena.svm_decoding;
                        end
                        
                        save(processedFile, 'ms', 'params', 'analysis_version', '-v7.3');
                        updated_files{length(updated_files) + 1, 1} = processedFile;
                        %
                        fprintf(' Done! -> %s\n', processedFile)
                    else
                        fprintf('# debugging\n')
                        fprintf('\t\t%s...\n', fname)
                        updated_files{length(updated_files) + 1, 1} = processedFile;
                    end
                    fprintf('\t\ttime: %s  \n\t\tsession process time -- %.0f seconds\n\n', datetime(), toc(timer_session))
                else
%                     fprintf('~~~PROCESSED file RERUN skipped: %s\n', fname)
%                     fprintf('\tsame analysis_version : %s\n', prev_version.analysis_version)
                    skipfiles{length(skipfiles) + 1,1} = processedFile;
                end
%                 if length(updated_files) > 15
%                     delete(gcp)
%                     clear
%                     APA_analysis_batch_reprocess_files
%                 end
%                 length(updated_files)
            end % % RERUN PROCESSED FILES
            fprintf(' \t\t\t%s Done! \n\t\t#######################\n\t\t#######################\n\t\t#######################\n', animal_name)
            fprintf('total process time -- %.1f minutes\n\n\n', toc(timer_all)/60)
        end
    end
end
% delete(gcp)
diary off