function [behav, behav_params] = read_APA_csv(fname_position, fname_params)
python_posfile = any(contains(fname_position, 'behav_position_data.csv')) == 1;
apaDAT = false(1,4); 
apaDAT(1) = any(contains(fname_position, '_arena.dat')) == 1;
apaDAT(2) = any(contains(fname_position, '_Arena.dat')) == 1;
apaDAT(3) = any(contains(fname_position, '_room.dat')) == 1;
apaDAT(4) = any(contains(fname_position, '_Room.dat')) == 1;
BSG_posfile_DAT = any(apaDAT==1);

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
    headerSize = 28; % 28
    maxheader = headerSize+60;
%     fname_position = 'C:\Users\gjb326\Desktop\RecordingData\AlejandroGrau\datfiles\V4Min\GRAU_V4Miniscope_TestMouse1_CaImaging_Training1_Arena.dat';
    temp = readtable(fname_position, 'HeaderLines', headerSize);
    sizechk = size(temp,2)==10 | size(temp,2)==9;
    okgo = sizechk && temp.Var1(1) == 1;
    while ~okgo && headerSize<=maxheader
        headerSize = headerSize+1;
        try
            temp = readtable(fname_position, 'HeaderLines', headerSize);
        catch
%             warning('Couldn''t read, trying again with %d header lines', headerSize+1) 
        end
        sizechk = size(temp,2)==10 | size(temp,2)==9;
        okgo = sizechk && temp.Var1(1) == 1;
    end
    if ~okgo && headerSize>maxheader
        error('Couldn''t read in file \n %s \n using header line up to %d', fname_position, headerSize)
    else
    behav               = [];
    behav.fileName      = fname_position;
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
    behav_params            = [];
    behav_params.headerSize = headerSize;
    end
end

end