function [behav, behav_params] = read_APA_csv(fname_position, fname_params)
    %% Read in the extracted position from APA_extract_pos_batch.ipynb
    % fname_position = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_20\17_24_27\experiment\behav_position_data.csv';
    behav = readtable(fname_position);
    %% Decode the .json file for parameters
    % fname_params = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_20\17_24_27\experiment\behav_ext_params.json';
    behav_params = readJSON(fname_params);
%     fid = fopen(fname_params); % Opening the file
%     raw = fread(fid,inf); % Reading the contents
%     str = char(raw'); % Transformation
%     fclose(fid); % Closing the file
%     P2 = struct2table(jsondecode(str));
% 
%     params = jsondecode(str); % Using the jsondecode function to parse JSON from string
%     names = P2.Properties.VariableNames;
%     for i = 1:length(names)
%         eval(sprintf('temp = str2double(P2.%s);',  names{i}))
%         if isnan(temp) % its prob a tuple
%             eval(sprintf('temp = P2.%s;',  names{i}))
%             a = strfind(temp, '(');
%             b = strfind(temp, ',');
%             c = strfind(temp, ')');
%             if ~isempty(a) % it's a tuple
% %                 warning('\ntuple found  -  %s = %s', names{i}, temp)
%                 inds = [a b c]; 
%                 temp2 = [str2double(temp(inds(1)+1:inds(2)-1))];
%                 for j = 1:length(b)
%                     temp2 = [temp2, str2double(temp(inds(j+1)+1:inds(j+2)-1))];
%                 end
%                 temp = temp2;
%             else % not a tuple
%                 warning('\nBad value!  -  %s', names{i})
%                 temp = NaN;
%             end
%         end
%         eval(sprintf('params.%s = temp;', names{i}))
%     end
end