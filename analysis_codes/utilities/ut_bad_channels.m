
%%% Quick script that has BAD Channels based on Ripple char: 


%%% Bad Channels

artfi_subj_chan = cell(1,19);
artfi_subj_chan{1}  = [];   %NS127_02 
artfi_subj_chan{2}  = [];   %NS128_02 
artfi_subj_chan{3}  = {'RDa2-RDa3','RDa3-RDa4','RDa4-RDa5'};  
artfi_subj_chan{4}  = {'RDa2-RDa3',};   %NS135
artfi_subj_chan{5}  = {'RDh2-RDh3'};  %NS136
artfi_subj_chan{6}  = {'LDh3-LDh4','LDh4-LDh5'}; %NS137
artfi_subj_chan{7}  = []; %NS138
artfi_subj_chan{8}  = {'LDa2-LDa3','LDa5-LDa6'}; %NS140_2
artfi_subj_chan{9}  = {'LDa2-LDa3'}; %NS140
artfi_subj_chan{10} = {'RDa3-RDa4','RDa5-RDa6','RDh3-RDh4','RDh4-RDh5','RDh5-RDh6','LDa2-LDa3','LDa3-LDa4','LDa4-LDa5'};  %NS142
artfi_subj_chan{11} = {'LDa1-LDa2'};  %NS144_02
artfi_subj_chan{12} = {'RDa3-RDa4','RDa4-RDa5'};  %NS149
artfi_subj_chan{13} = [];  %NS153
artfi_subj_chan{14} = [];  %NS155
artfi_subj_chan{15} = [];  %NS166
artfi_subj_chan{16} = [];  %NS167
artfi_subj_chan{17} = {'LDh2-LDh3'};  %NS173
artfi_subj_chan{18} = [];  %NS174_02
artfi_subj_chan{19} = [];  %NS174_03

% artfi_subj_chan{4} = {'RDa2-RDa3'};
% artfi_subj_chan{5} = {'LDa1-LDa2','LDa3-LDa4','LDh2-LDh3','LT1-LTi2','RDh2-RDh3'}; 
% artfi_subj_chan{7} = {'LDa1-LDa2'};
% artfi_subj_chan{8} = {'LDa2-LDa3'};
% artfi_subj_chan{10} = {'LDa2-LDa3','RDa3-RDa4','RDa5-RDa6','RDh3-RDh4','RDh4-RDh5','RDh5-RDh6','LDa4-LDa5','LDa3-LDa4'}; % ALL! WTF!!!
% artfi_subj_chan{11} = {'LDa2-LDa3''LDa3-LDa4'};
% artfi_subj_chan{12} = {'RDa3-RDa4','RDa4-RDa5','RDh3-RDh4'}; %Qute ok! Dont need too much
% artfi_subj_chan{13} = {'RDa2-RDa3'};
% artfi_subj_chan{14} = {'RDa1-RDa2''RDa5-RDa6'};
% artfi_subj_chan{17} = {'LDh1-LDh2','LDh2-LDh3'}; % Important!
% artfi_subj_chan{18} = {'LDh1-LDh2','LDh2-LDh3'};
% artfi_subj_chan{19} = {'LDa3-LDa4'};




%%% Bad Trials
artfi_subj_chan_tr = cell(1,19);
% NS127_02
artfi_subj_chan_tr{1}{2} =  [35,50,54,55,56,82,83,91,124,125];
artfi_subj_chan_tr{1}{3} =  [5,19,42,79,114,124,152,];
artfi_subj_chan_tr{1}{4} =  [5,19,42,79,114,124,152,];
% NS135
artfi_subj_chan_tr{4}{1} =  [2,4,16,17,26,30,37,44,84,...
                             86,95,134,135,151,188];
artfi_subj_chan_tr{4}{2} =  [34,45,50,93,107,108,117,118,...
                             124,136,137138,150];
artfi_subj_chan_tr{4}{3} =  [28,30,101,103,118,130,133];
% NS136
artfi_subj_chan_tr{5}{1} =  [35,36,44,47,49,60,61,69,76];
artfi_subj_chan_tr{5}{2} =  [47,76,77,79,123,131,132,167];
artfi_subj_chan_tr{5}{5} =  [2,18,17,27,37,41,43,95,96,103,...
                            104,107,117,129,123,125,141];
artfi_subj_chan_tr{5}{7} = [41,42,43,56,55,57,58,63,66,69,79,...
                            84,86,89,90,93,94,125,127,130,156,157];
artfi_subj_chan_tr{5}{8} = [9,79,80,86,106,108,133];
artfi_subj_chan_tr{5}{9} = [9,20,21,79,80,118];
% NS137
artfi_subj_chan_tr{6}{1} = [32,34,35,62,69,70,71,72,73,74,75,77,78,...
                            81,83,84,86,95,99,101,105,110,111,113,...
                            114,37,139,160,162,165,171,174,198,205,...
                            210,212,213,214,216,218,223,224,233,234,...
                            235,236,237,243,245,246]; 
artfi_subj_chan_tr{6}{2} = [19,20,21,23,26,39,55,56,57,58,59,60,63,...
                            68,74,77,78,85,87,88,89,91,92,93,95,96,...
                            111,112,114,115,118,119,126,127,170,174,...
                            175,1776,177,184,188,193,204,206,210,213,...
                            214,218,219,223,224,226,233,234,236,237]; 
artfi_subj_chan_tr{6}{3} = [2,19,20,21,25,26,29,33,34,35,36,37,38,39,...
                            51,52,53,59,63,72,77,78,80,94,96,103,106,...
                            108,133,132,134,135,137,138,139,144,146,...
                            153,155,157,163,175,200,203,205,212,215,...
                            223,228,230,231,232,233,235,236,237,251,...
                            282,288,302,303,304,310,311,312]; 
artfi_subj_chan_tr{6}{4} = [17,19,21,22,23,24,30,31,34,46,54,68,76,...
                            90,94,111,114,117,118,122,128,130,131,...
                            133,135,136,138,141,142,143,146,147,173,...
                            175,178,188,191,195,204,214,220,223,229,...
                            232,237,238,240,243,244,255,269,274,275,...
                            279,281,284,287,290,297,300,302,303,305,...
                            307,308,326,330,]; 
artfi_subj_chan_tr{7}{1} = [59,61,90,93,102,103,109,131,134]; 
artfi_subj_chan_tr{8}{2} = [40,65]; 
artfi_subj_chan_tr{10}{1} = [9,10,14,15,17,18,40,41];
artfi_subj_chan_tr{10}{2} = [12,13,20,22];
artfi_subj_chan_tr{10}{3} = [2,3,4,16,18,19,24];
artfi_subj_chan_tr{10}{4} = [5,6,7,16,17,24,26,29,31,34,...
                            37,38,39,40,41,42,45];
artfi_subj_chan_tr{10}{5} = [11,13,27,29,33,36,40,42,56,57,59,64,66];
artfi_subj_chan_tr{10}{6} = [1,2,3,4,6,7,8,9,10,11,13,....
                            15,22,23,24,26,27,29,31,35,...
                            36,38,39,40,41,42,43]; % Remove
artfi_subj_chan_tr{10}{7} = [3,4,8,9,14,25,27,28,31,32,35];
artfi_subj_chan_tr{10}{8} = [8,12,13,16,20,21,23,25,26,28,29,31,...
                            37,38,39,41,42,43,44,45,46,47,48,49];
artfi_subj_chan_tr{11}{2} = [11,17,27,28,42,48,49,51,52,80,...
                             92,93,94,95,96,97,98,137,141,...
                             142,153,154,155,156,157,158,186,187,188];
artfi_subj_chan_tr{11}{1} = [15,10,16,67,70,96,112,114,116,118,134,...
                            135,142,1143,147,148,164,168,173174,192,193];
artfi_subj_chan_tr{11}{3} = [27,33,54,62,82,93,106,107,108,110,...
                            111,142,143,144,153,163,164,165,166,192];
artfi_subj_chan_tr{12}{1} = [58,123,129,141,143,150];
artfi_subj_chan_tr{12}{2} = [6,13,69,97,113,];
artfi_subj_chan_tr{12}{3} = [25,26,37];
artfi_subj_chan_tr{12}{4} = [28,29,64,73];
artfi_subj_chan_tr{12}{5} = [15,16,19,100,];
artfi_subj_chan_tr{13}{3} = [16,20,23,24,90,91,151];
artfi_subj_chan_tr{14}{2} = [11,14,39,90];
artfi_subj_chan_tr{14}{4} = [48,131];
artfi_subj_chan_tr{14}{6} = [64];
artfi_subj_chan_tr{17}{1} = [4,11,12,13,14,16,17,19,20,24,27,33,...
                            39,41,50,60,63,64,67,69,71];
artfi_subj_chan_tr{17}{2} = [3,4,5,12,14,22,24,25,27,29,31,33];
artfi_subj_chan_tr{18}{2} = [74,76,104,109];
artfi_subj_chan_tr{18}{3} = [40,42,81,103,105];
artfi_subj_chan_tr{19}{2} = [57,58,59,61,6263,64,67];