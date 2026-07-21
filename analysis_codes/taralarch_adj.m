function M3x4 = read_talxfm(xfmFile)
% READ_TALXFM Read a talairach .xfm file and return the 3x4 numeric matrix.
%   M3x4 = read_talxfm(xfmFile)
% Looks specifically for "Linear_Transform = <numbers>;" and parses those numbers.
% Returns 3x4 matrix in the same row-wise order as in the file.

txt = fileread(xfmFile);
txt = strrep(txt, ',', '.'); % safety for comma decimals

% try to capture the block after "Linear_Transform =" up to the next semicolon
tok = regexp(txt, 'Linear_Transform\s*=\s*([^;]+)', 'tokens', 'once', 'ignorecase');

if ~isempty(tok)
    num_block = tok{1};
else
    % fallback: use entire file if the pattern not found
    num_block = txt;
end

% extract numeric tokens (supports decimals, signs, scientific notation)
numTokens = regexp(num_block, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
nums = str2double(numTokens);

if numel(nums) < 12
    error('read_talxfm: need at least 12 numeric entries; found %d', numel(nums));
end

% Take the first 12 numbers in file order and arrange row-wise into 3x4
v = nums(1:12);
M3x4 = [ v(1:4); v(5:8); v(9:12) ];

end


found = 'C:\Users\barak\Documents\Predictive_Ripples_2025\data\elec_coordinates\NS128_02_talairach.xfm';
M3x4 = read_talxfm(found);
disp(M3x4);