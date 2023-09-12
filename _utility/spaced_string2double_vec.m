function [double_vec] = spaced_string2double_vec(string_vec)
spc = strfind(string_vec, ' ');
b = [spc-1 length(string_vec)];
a = [1 spc+1];
double_vec=[];
for iii=1:length(spc)+1
    double_vec = [double_vec, str2double(string_vec(a(iii):b(iii)))];
end
