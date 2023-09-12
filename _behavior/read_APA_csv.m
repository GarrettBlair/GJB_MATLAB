function [behav, behav_params] = read_APA_csv(fname_position, fname_params)
python_posfile = any(contains(fname_position, 'behav_position_data.csv')) == 1;
apaDAT = false(1,4);
apaDAT(1) = any(contains(fname_position, '_arena.dat')) == 1;
apaDAT(2) = any(contains(fname_position, '_Arena.dat')) == 1;
apaDAT(3) = any(contains(fname_position, '_room.dat')) == 1;
apaDAT(4) = any(contains(fname_position, '_Room.dat')) == 1;
BSG_posfile_DAT = any(apaDAT==1);
biosig_time_fmt = 'dd.MM.yyyy hh:mm a';

apaMAT = false(1,4);
apaMAT(1) = any(contains(fname_position, '_arena.mat')) == 1;
apaMAT(2) = any(contains(fname_position, '_Arena.mat')) == 1;
apaMAT(3) = any(contains(fname_position, '_room.mat')) == 1;
apaMAT(4) = any(contains(fname_position, '_Room.mat')) == 1;
BSG_posfile_MAT = any(apaMAT==1);
if python_posfile
    %% Read in the extracted position from APA_extract_pos_batch.ipynb
    % fname_position = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_20\17_24_27\experiment\behav_position_data.csv';
    behav = readtable(fname_position);
    %% Decode the .json file for parameters
    % fname_params = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_20\17_24_27\experiment\behav_ext_params.json';
    behav_params = readJSON(fname_params);
elseif BSG_posfile_MAT
    %%
    % preprocessed files from APA_DAT_preprocessing.py script
    temp = load(fname_position);
    behav               = [];
    behav.fileName      = fname_position;
    behav.frameNum      = temp.data(:,1);  % FrameCount
    behav.timestamps    = double(temp.data(:,2));  % 1msTimeStamp
    if behav.timestamps(1) ~= 0
        error('FrameNum not starting at 0; check header init')
    end
    behav.rawx          = double(temp.data(:,3));  % RoomX
    behav.rawy          = double(temp.data(:,4));  % RoomY
    behav.Sector        = temp.data(:,5);  % Sectors
    behav.State         = temp.data(:,6);  % State
    behav.CurrentLevel  = temp.data(:,7);  % CurrentLevel
    behav.MotorState    = temp.data(:,8);  % MotorState
    if size(temp.data,2) == 9
        behav.frameInfo     = temp.data(:,9); % FrameInfo
    elseif size(temp.data,2) == 10
        behav.frameInfo     = temp.data(:,10); % FrameInfo
    end
    behav_params            = [];
    behav_params.headerSize = NaN;
elseif BSG_posfile_DAT
    %% Read in the extracted position from APA_extract_pos_batch.ipynb
    %     fname_position = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\datfiles\TM2_pretrain_arena.dat';
    %     fname_position = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\datfiles\V4Min\GRAU_V4Miniscope_TestMouse1_CaImaging_Training1_Arena.dat';
    %     temp = readtable(fname_position);
    %     try
    opts = detectImportOptions(fname_position);
    temp = readtable(fname_position, 'HeaderLines', 0, 'ReadVariableNames', false);
    heads = temp.Var1(1:opts.DataLines(1)*5);
    headerSize = opts.DataLines(1)-1;
    temp = readtable(fname_position, 'HeaderLines', headerSize, 'ReadVariableNames', false);
    counts=0;
    DAT_fileinfo = [];
    for ii = 1:length(heads)
        counts;
        s = heads{ii};
        d_start = strfind(s, '(')+2;
        d_stop = strfind(s, ')')-2;
        d = s(d_start:d_stop);
        if contains(s, '//')
            % ignore comment
        elseif contains(s, '%%END_HEADER')
            if counts<8
                warning('Not all expected fields were found (only found n=%d', counts)
            end
            break
        elseif contains(s, '%%BEGIN_HEADER')
            
            badfile=false;
        elseif contains(s, '%Date.0') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.Date = d;
        elseif contains(s, '%Time.0') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.Time = d;
            datevec = [DAT_fileinfo.Date ' ' DAT_fileinfo.Time];
            DAT_fileinfo.datetime = datetime(datevec, 'InputFormat', biosig_time_fmt);
        elseif contains(s, '%ElapsedTime_ms.0') %%%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.ElapsedTime_ms = str2double(d);
        elseif contains(s, '%ShockParameters.0') %%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.ShockParameters = spaced_string2double_vec(d);
        elseif contains(s, '%ElapsedPath_m.0') %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.ElapsedPath_m = str2double(d);
        elseif contains(s, '%ArenaDiameter_m.0') %%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.ArenaDiameter_m = str2double(d);
        elseif contains(s, '%TrackerResolution_PixPerCM.0') %%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.TrackerResolution_PixPerCM = str2double(d);
        elseif contains(s, '%ArenaCenterXY.0') %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.ArenaCenterXY = spaced_string2double_vec(d);
        elseif contains(s, '%Frame.0') %%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            DAT_fileinfo.Frame = d;
        elseif contains(s, '%ReinforcedSector.0') %%%%%%%%%%%%%%%%%%%%%%%%%
            counts=counts+1;
            d = d(1:end-1);
            DAT_fileinfo.ReinforcedSector = spaced_string2double_vec(d);
        elseif contains(s, '%RoomTrackReinforcedSector.0') %%%%%%%%%%%%%%%%
            counts=counts+1;
            d = d(1:end-1);
            DAT_fileinfo.RoomTrackReinforcedSector = spaced_string2double_vec(d);
        elseif contains(s, '%ArenaTrackReinforcedSector.0') %%%%%%%%%%%%%%%%
            counts=counts+1;
            d = d(1:end-1);
            DAT_fileinfo.RoomTrackReinforcedSector = spaced_string2double_vec(d);
        elseif contains(s, '%Sample.0 (') %%%%%%%%%%%%%%%%
            d = s(strfind(s, '(')+2:strfind(s, ')')-2);
            isspc = strfind(d, ' ');
            varnames = {d(1:isspc(1)-1)};
            for i2 = 1:length(isspc)-1
                ss = d(isspc(i2)+1:isspc(i2+1)-1);
                if strcmp(ss, 'RoomX') || strcmp(ss, 'ArenaX')
                    ss = 'rawx';
                elseif strcmp(ss, 'RoomY') || strcmp(ss, 'ArenaY')
                    ss = 'rawy';
                elseif strcmp(ss, '1msTimeStamp')
                    ss = 'timestamps';
                elseif strcmp(ss, 'FrameCount')
                    ss = 'frameNum';
                end
                varnames{i2+1,1} = ss;
                
            end
            nvar = length(varnames);
        end
    end
    %     catch
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     headerSize = 28; % 28
    %     maxheader = headerSize+60;
    % %     fname_position = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\datfiles\V4Min\GRAU_V4Miniscope_TestMouse1_CaImaging_Training1_Arena.dat';
    %     temp = readtable(fname_position, 'HeaderLines', headerSize, 'ReadVariableNames', false);
    %     sizechk = size(temp,2)==10 | size(temp,2)==9;
    %     okgo = sizechk && temp.Var1(1) == 1;
    %     while ~okgo && headerSize<=maxheader
    %         headerSize = headerSize+1;
    %         try
    %             temp = readtable(fname_position, 'HeaderLines', headerSize, 'ReadVariableNames', false);
    %         catch
    % %             warning('Couldn''t read, trying again with %d header lines', headerSize+1)
    %         end
    %         sizechk = size(temp,2)==10 | size(temp,2)==9;
    %         okgo = sizechk && temp.Var1(1) == 1;
    %     end
    %     badfile = ~okgo && headerSize>maxheader;
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if badfile == true
        error('Couldn''t read in file \n %s \n using header line up to %d', fname_position, headerSize)
    else
        behav               = [];
        behav.fileName      = fname_position;
        behav.DAT_fileinfo  = DAT_fileinfo;
        
        if exist('nvar', 'var')
            for i2 = 1:nvar
                eval(sprintf('behav.%s = temp.Var%d;', varnames{i2}, i2))
            end
        else
            warning('using old method without header Samle.0 names')
            behav.frameNum      = temp.Var1;  % FrameCount
            behav.timestamps    = temp.Var2;  % 1msTimeStamp
            if behav.timestamps(1) ~= 0
                error('FrameNum not starting at 0; check header init')
            end
            behav.rawx          = temp.Var3;  % RoomX
            behav.rawy          = temp.Var4;  % RoomY
            behav.Sector        = temp.Var5;  % Sectors
            behav.State         = temp.Var6;  % State
            behav.CurrentLevel  = temp.Var7;  % CurrentLevel
            behav.MotorState    = temp.Var8;  % MotorState
            if size(temp,2) == 9
                behav.frameInfo     = temp.Var9; % FrameInfo
            elseif size(temp,2) == 10
                behav.frameInfo     = temp.Var10; % FrameInfo
            end
            
        end
        behav_params            = [];
        behav_params.headerSize = headerSize;
    end
end

    function [double_vec] = spaced_string2double_vec(string_vec)
        spc = strfind(string_vec, ' ');
        b = [spc-1 length(string_vec)];
        a = [1 spc+1];
        double_vec=[];
        for iii=1:length(spc)+1
            double_vec = [double_vec, str2double(string_vec(a(iii):b(iii)))];
        end
        
    end
end