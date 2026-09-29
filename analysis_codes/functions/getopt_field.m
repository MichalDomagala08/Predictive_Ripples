function val = getopt_field(s, field, default)
% GETOPT_FIELD  Return s.(field) if it exists and is non-empty, else default.
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end
