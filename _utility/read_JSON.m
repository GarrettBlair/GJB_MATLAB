function [json_data] = read_JSON(json_filename)
%% Read in .json file into a structure
% json_filename = 'C:\Users\gjb326\Desktop\RecordingData\GarrettBlair\APA_aquisition\Hipp16942\2022_06_20\17_24_27\experiment\behav_ext_params.json';
fid = fopen(json_filename); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
J_struct = struct2table(jsondecode(str));

json_data = jsondecode(str); % Using the jsondecode function to parse JSON from string
names = J_struct.Properties.VariableNames;
for i = 1:length(names)
    strName = sprintf('J_struct.%s', names{i});
    eval(sprintf('temp = str2double(%s);',  strName))
    if isnan(temp) % its prob a tuple or structure
        eval(sprintf('temp = J_struct.%s;',  names{i}))
        if ischar(temp) % because I saved tuples as strings...
        a = strfind(temp, '(');
        b = strfind(temp, ',');
        c = strfind(temp, ')');
        if ~isempty(a) % it's a tuple
            %                 warning('\ntuple found  -  %s = %s', names{i}, temp)
            inds = [a b c];
            temp2 = [str2double(temp(inds(1)+1:inds(2)-1))];
            for j = 1:length(b)
                temp2 = [temp2, str2double(temp(inds(j+1)+1:inds(j+2)-1))];
            end
            temp = temp2;
%         else % not a tuple, prob just a string
%             warning('\nBad value!  -  %s', names{i})
%             temp = NaN;
        end
        elseif isstruct(temp)
            J_struct2 = struct2table(temp);
            names2 = J_struct2.Properties.VariableNames;
            for j = 1:length(names2)
                eval(sprintf('json_data.%s.%s = J_struct2.%s;', names{i}, names2{j}, names2{j}))
            end
        end
    end
    eval(sprintf('json_data.%s = temp;', names{i}))
end


