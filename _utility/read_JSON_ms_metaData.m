function [data_out] = read_JSON_ms_metaData(json_filename)
%% Read in .json file into a structure
% json_filename = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_20\17_24_27\experiment\behav_ext_params.json';
fid = fopen(json_filename); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data_out = jsondecode(str);
