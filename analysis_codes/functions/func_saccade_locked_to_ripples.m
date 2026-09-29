function ripplewise_table_s = func_saccade_locked_to_ripples(saccades_table,ripplewise_table,anmode)
    % This scriopts supplement ripple info with saccade or fixation
    % informaton such as region, its AWS/luminance, durations and so on! 

    %cols = [39,1,4,5,6,7,12:15,18:23,27:38,48,49];

    cols  = {'subject_new','name','trial_number','indx_sac','latency','duration','Xpx_original','Ypx_original','blinkSac','blinkTrial','firstInTrial','lastInTrial',...
             'amplitude','peakVelocity','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff','lum_r50px','aws_r20px','lum_r20px'...
             'aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs','fixation_beg','fixation_end'};
    % entropy_table = unique(aa(:, entropy_vars));
    temp_table = saccades_table(1:0,cols);
    temp_table.r_num = zeros(height(temp_table),1);
    temp_table.is_saccades = zeros(height(temp_table),1);
    temp_table.is_fixation = zeros(height(temp_table),1);
    
    for i = 1:height(ripplewise_table)
    
         curr_trial = saccades_table(saccades_table.subject_new == ripplewise_table(i,:).subject &  saccades_table.trial_number == ripplewise_table(i,:).trial ,cols);
         
         switch anmode
             case 0 % In exact saccade/ time
                 is_in_saccade = find(curr_trial.saccade_onset_time*500+2000<ripplewise_table(i,:).median & ripplewise_table(i,:).median < curr_trial.saccade_offset_time*500+2000);
                 is_in_fix = find(curr_trial.fixation_beg*500+2000<ripplewise_table(i,:).median & ripplewise_table(i,:).median < curr_trial.fixation_end*500+2000);
             case 1 % In mean saccade time (13s) 
                is_in_saccade = find(curr_trial.saccade_onset_time*500+2000 <ripplewise_table(i,:).median & ripplewise_table(i,:).median < 13+curr_trial.saccade_onset_time*500+2000);
                is_in_fix = find(curr_trial.fixation_beg*500+2000+13<ripplewise_table(i,:).median & ripplewise_table(i,:).median < curr_trial.fixation_end*500+2000);
             case 2 % In perisaccadic vs no Peri-saccadic 
                is_in_saccade = find(curr_trial.saccade_onset_time*500+2000 -25 <ripplewise_table(i,:).median & ripplewise_table(i,:).median < curr_trial.saccade_offset_time*500+2000+25);
                is_in_fix = find(curr_trial.fixation_beg*500+2000+25<ripplewise_table(i,:).median & ripplewise_table(i,:).median < curr_trial.fixation_end*500+2000-25);
         end
         curr_inf = curr_trial([is_in_fix,is_in_saccade],:);
    
         if ~isempty(curr_inf)
             curr_inf.r_num = ripplewise_table(i,:).rip_num;
              if is_in_fix; curr_inf.is_fixation =1;  curr_inf.is_saccades =0; end;      if is_in_saccade; curr_inf.is_saccades =1;  curr_inf.is_fixation =0; end
    
             temp_table = [temp_table; curr_inf];
         else
             newRow = repmat({NaN}, 1, width(temp_table));  % wypełnij NaNami
            newRow{1} = ripplewise_table.subject(i);        % subject_new
            newRow{2} = curr_trial.name(1);                  % Image
            newRow{3} = ripplewise_table.trial(i);          % trial_number
            newRow{end-2} = ripplewise_table.rip_num(i);      % r_num
            newRow{end-1} = 0;    newRow{end} = 0;    
    
             temp_table = [temp_table; table(newRow{:}, 'VariableNames', temp_table.Properties.VariableNames)];
         end
    end
    
    
    ripplewise_table_sac = outerjoin(ripplewise_table,temp_table,'LeftKeys',{'subject','trial','rip_num'},'RightKeys',{'subject_new','trial_number','r_num'});
    ripplewise_table_sac.is_before_trial = ripplewise_table_sac.median < 2000;
    ripplewise_table_sac.is_after_trial = ripplewise_table_sac.median > 5000;
    ripplewise_table_sac.is_during_trial = ~ripplewise_table_sac.is_after_trial & ~ripplewise_table_sac.is_before_trial;
    
    ripplewise_table_sac.is_fixation = logical(ripplewise_table_sac.is_fixation); ripplewise_table_sac.is_saccades = logical(ripplewise_table_sac.is_saccades); 
    
    ripplewise_table_sac.occular_beg = zeros(height(temp_table),1); ripplewise_table_sac.occular_end = zeros(height(temp_table),1);
    
    ripplewise_table_sac(ripplewise_table_sac.is_fixation,:).occular_beg = ripplewise_table_sac(ripplewise_table_sac.is_fixation,:).fixation_beg;
    ripplewise_table_sac(ripplewise_table_sac.is_fixation,:).occular_end = ripplewise_table_sac(ripplewise_table_sac.is_fixation,:).fixation_end;
    ripplewise_table_sac(ripplewise_table_sac.is_saccades,:).occular_beg = ripplewise_table_sac(ripplewise_table_sac.is_saccades,:).saccade_onset_time;
    ripplewise_table_sac(ripplewise_table_sac.is_saccades,:).occular_end = ripplewise_table_sac(ripplewise_table_sac.is_saccades,:).saccade_offset_time;
    ripplewise_table_sac.occular_duration = ripplewise_table_sac.occular_end - ripplewise_table_sac.occular_beg;
    
    ripplewise_table_sac(ripplewise_table_sac.is_fixation,[16,19,23,24]) = num2cell(nan(sum(ripplewise_table_sac.is_fixation),4));
    
    ripplewise_table_s = ripplewise_table_sac(:,[1:7,9,8,10,44:46,42,43,19:22,14,47:49,17,18,15,16,23,24,27:38]);

end