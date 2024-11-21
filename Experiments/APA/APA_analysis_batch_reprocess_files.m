% analysis to do:
% correlation - xcorr - grainger causality - canonical corr analysis
    % between ACC and HPC ipos
% behav var to look into - total path length, 1st/2nd half in water sess,
% inter-entrance times
clear
analysis_version = 'v2.2';
diarydir  = 'C:\Users\gjb326\Documents\MATLAB\GJB_MATLAB\diaries\';
diaryname = sprintf('%stemp_reprocessing_%s', diarydir, analysis_version);
diary(diaryname)
params = [];
APA_rat_imaging_params_current;
% APA_PKCZmouse_imaging_params_current;

fprintf('rand shuff pcells = %d\n',params.num_random_shuffle_pcell)
fprintf('rand shuff decode = %d\n',params.num_random_shuffle_decode)

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
fprintf('%s\n', datetime())
%
% animals = {'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'Acc20832', 'Acc19947', 'Hipp18240'};% animals = {'Acc20832', 'Acc19947'};
animals = {'HPCACC24504', 'HPCACC24505', 'HPCACC24514'};% animals = {'Acc20832', 'Acc19947'};
fprintf('\t\tanalysis ver = %s\n', analysis_version)
% experiment_folder = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA_aquisition/';
experiment_folder = [];
% animals = {'mHPC23454', 'mHPC23459'};% animals = {'Acc20832', 'Acc19947'};
% experiment_folder = {'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\PKCZ_imaging\'};
% experiment_folder{1} = 'C:/Users/gjb326/Desktop/RecordingData/GarrettBlair/APA/';

animals = {'Hipp18240', 'Acc20832', 'Acc19947', 'HPCACC24500', 'HPCACC24502', 'HPCACC24503', 'HPCACC24504', 'HPCACC24505', 'HPCACC24514'};% animals = {'Acc20832', 'Acc19947'};

experiment_folder{1} = 'D:/GarrettBlair/APA/';
% animals = {'Acc20832', 'Acc19947', 'Hipp18240'};
numAnimals = length(animals);
% experiment_folder{2} = 'D:/APA recordings/';

% DAT_Dir             = sprintf('%sDAT_files/', experiment_folder);

rerun_behav         = false;
rerun_place         = false;
resave_proccessed   = false;
resave_contours     = false;
resave_simple       = true;
rerun_processed     = false;
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
            simple_dir = sprintf('%s%s/simple_files/', experiment_folder{expLoop}, animal_name);
            fprintf('\n\nREPROCESS files for %s in folder  \n%s...\n', animal_name, processed_dir)
%             temp = dir([processed_dir '*miniscope*']);
            temp = dir([processed_dir '*@placecells*']);
%             temp2 = dir([processed_dir '*@placecells*']);
            nfiles = length(temp);
            for sessionLoop = 1:nfiles
                fname = temp(sessionLoop).name;
                processedFile   = [temp(sessionLoop).folder '/' fname];
                simpleFile      = [simple_dir '/' fname];
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
                                                
                        %                         [ms] = APA_generate_placemaps(tempms.ms, params);
                        ms = tempms.ms;
                        raw_traces = ms.neuron.C + ms.neuron.YrA;
                        if sessionLoop==1
                            plotthisshit = true;
                        else
                            plotthisshit = false;
                        end
                        [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, plotthisshit);
                        
                        ms.spks = normalize_rows(dff_filt);
                        ms.spks(isnan(ms.spks)) = 0;
                        clearvars dff_filt raw_filt raw_traces
                        
                        fprintf('\t\t   [placefields]: \t\t%s  \n', datetime())
%                         [ms] = APA_generate_polar_placemaps(ms, params);
                        [ms] = APA_generate_placemaps(ms, params);
                        
                        fprintf('\t\t   [ipos]: \t\t%s  \n', datetime())
                        ms.room.momentary_pos_info = []; ms.room.conjoint_ipos_min = []; ms.room.conjoint_ipos_av = [];
                        ms.arena.momentary_pos_info = []; ms.arena.conjoint_ipos_min = []; ms.arena.conjoint_ipos_av = [];
%                         [ms.room.momentary_pos_info, ~, ms.room.conjoint_ipos_min, ms.room.conjoint_ipos_av] = ...
%                             Fenton_ipos_testing(ms, params.ipos_int_time, 'room_polar',  params);
                        [ms.room.momentary_pos_info, ~, ms.room.conjoint_ipos_min, ms.room.conjoint_ipos_av] = ...
                            Fenton_ipos_testing(ms, params.ipos_int_time, 'room',  params);
%                         [ms.arena.momentary_pos_info, ~, ms.arena.conjoint_ipos_min, ms.arena.conjoint_ipos_av] = ...
%                             Fenton_ipos_testing(ms, params.ipos_int_time, 'arena_polar', params);
                        [ms.arena.momentary_pos_info, ~, ms.arena.conjoint_ipos_min, ms.arena.conjoint_ipos_av] = ...
                            Fenton_ipos_testing(ms, params.ipos_int_time, 'arena', params);
                        ir = abs(ms.room.momentary_pos_info); 
                        ir(isnan(ir)) = 0;
                        ia = abs(ms.arena.momentary_pos_info); 
                        ia(isnan(ia)) = 0;
%                         ipos = abs(ms.room.momentary_pos_info) - abs(ms.arena.momentary_pos_info);
                        ipos = ir - ia; 
                        ipos(isnan(ms.room.momentary_pos_info) & isnan(ms.arena.momentary_pos_info)) = NaN;
                        ipos(isnan(ipos))=0;
                        iposm = nanmean(ipos,1);
                        if plotthisshit==true
                        figure; imagesc(ipos, [-.1 .1]); colormap redblue;
                        end
                        fprintf('\t\t   [svm decode]: \t\t%s  \n', datetime())
                        [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                        
                        
                        fprintf('\t\t   [seg corr]: \t\t%s  \n', datetime())
                        spks_1secbin = bin_spks_time(ms.spks>0, 1, ms.timestamps./1000, false);
                        spks_1secbin(isnan(spks_1secbin)) = 0;
                        corrtype = 'Kendall';
                        [c,p] = corr(spks_1secbin','Type', corrtype);
                        segs_corr = [];
                        segs_corr.c = c;
                        segs_corr.p = p;
                        segs_corr.type = corrtype;
                        segs_corr.tau = '1sec';
                        segs_corr.binspks = spks_1secbin;
                        
                        spks1 = spks_1secbin;
                        spks2 = spks1(:, floor(size(spks1,2)/2)+1:end);
                        spks1 = spks1(:, 1:floor(size(spks1,2)/2));
                        
                        [c1,p1] = corr(spks1', 'Type', corrtype);
                        [c2,p2] = corr(spks2', 'Type', corrtype);
                        segs_corr.split1_corr = c1;
                        segs_corr.split1_p = p1;
                        segs_corr.split2_corr = c2;
                        segs_corr.split2_p = p2;
                        
                        
%                         fprintf('\t\trunning ipos...\n')
%                         [ms.room.momentary_pos_info,  rtemp, ms.room.conjoint_ipos_min, ms.room.conjoint_ipos_av]  =  Fenton_ipos(ms, params.ipos_int_time, 'room', params);
%                         [ms.arena.momentary_pos_info, atemp, ms.arena.conjoint_ipos_min, ms.arena.conjoint_ipos_av] = Fenton_ipos(ms, params.ipos_int_time, 'arena', params);
                        
%                         fprintf('\t\trunning svm decode...\n')
%                         [ms.room.svm_decoding, ms.arena.svm_decoding] = APA_within_sess_decoding(ms, params);
                        
                        save(processedFile, 'ms', 'segs_corr', 'params', 'analysis_version', '-v7.3');
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
            end % % RERUN PROCESSED FILES
            fprintf(' \t\t\t%s Done! \n\t\t#######################\n\t\t#######################\n\t\t#######################\n', animal_name)
            fprintf('total process time -- %.1f minutes\n\n\n', toc(timer_all)/60)
        end
    delete(gcp)
    end
end
% delete(gcp)



%%
PRINT_ONLY = false;
this_ver = str2double(analysis_version(2:end));
skipfiles = []; updated_files = [];
timer_all = tic;
if resave_simple == true
    for animalLoop = 1:numAnimals
        for expLoop = 1:length(experiment_folder)
            animal_name = animals{animalLoop};
%             processed_dir = sprintf('%s%s/processed_files/', experiment_folder{expLoop}, animal_name);
%             simple_dir = sprintf('%s%s/simpleFiles/', experiment_folder{expLoop}, animal_name);
            simple_dir = sprintf('%s%s/simple_files/', experiment_folder{expLoop}, animal_name);
            fprintf('\n\nREPROCESS simple files for %s in folder  \n%s...\n', animal_name, simple_dir)
            temp = dir([simple_dir '*@placecells*']);
            nfiles = length(temp);
            for sessionLoop = 1:nfiles
                fname = temp(sessionLoop).name;
                simpleFile   = [temp(sessionLoop).folder '/' fname];
                timer_session = tic;
                if ~PRINT_ONLY
                    fprintf('~~~SIMPLE file rerunning for: %s...\n', fname)
                    fprintf('\t\tanalysis_version : %s\n', analysis_version)
                    fprintf('\t\tstart time: %s  \n', datetime())

                    old = load(simpleFile, 'craw');
                    raw_traces = old.craw;
                    if mod(sessionLoop, 5) == 1
                        plotthisshit = true;
                    else
                        plotthisshit = false;
                    end
                    [dff_filt, ~] = ca_trace_dff_filter(raw_traces, 30, .1, 1.5, plotthisshit);

                    spks = normalize_rows(dff_filt);
                    spks(isnan(spks)) = 0;
                    save(simpleFile, 'spks', '-append');
                    updated_files{length(updated_files) + 1, 1} = simpleFile;

                    clearvars spks dff_filt raw_traces old

                else
                    fprintf('# debugging\n')
                    fprintf('\t\t%s...\n', fname)
                    updated_files{length(updated_files) + 1, 1} = simpleFile;
                end
                fprintf('\t\ttime: %s  \n\t\tsession process time -- %.0f seconds\n\n', datetime(), toc(timer_session))
            end % % RERUN PROCESSED FILES
            fprintf(' \t\t\t%s Done! \n\t\t#######################\n\t\t#######################\n\t\t#######################\n', animal_name)
            fprintf('total process time -- %.1f minutes\n\n\n', toc(timer_all)/60)
        end
    end
end

%%




diary off