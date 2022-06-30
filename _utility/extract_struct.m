function [varargout] = extract_struct(struct_in)
    fields = fieldnames(struct_in);
    varargout = cell(size(fields,1), 1);
    for i = 1:size(fields,1)
        eval(sprintf('varargout{i} = struct_in.%s', fields{i}));
    end
end