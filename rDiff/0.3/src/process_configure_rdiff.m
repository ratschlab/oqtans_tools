function CFG = process_configure_rdiff(CFG)
  %  process_configure_rDiff(CFG)

for i=1:length(CFG.BAM_FILES)
    if isempty(CFG.NAMES{i})
        CFG.NAMES{i}=CFG.BAM_FILES{i};
    end
end

%Get the paths of the filenames
for i=find(CFG.SAMPLES==1)
  CFG.BAM_FILES{i}=fullfile(CFG.data_dir,CFG.BAM_FILES{i}); 
end

for i=find(CFG.SAMPLES==2)
  CFG.BAM_FILES{i}=fullfile(CFG.data_dir,CFG.BAM_FILES{i}); 
end


%Check that the variance functions are available if necessary
if or(CFG.perform_nonparametric,CFG.perform_parametric) &&  length(CFG.predefined_variance_function1)<3 && isempty(CFG.variance_function_1)
    CFG.compute_variance_function_1=1;
end
    
if or(CFG.perform_nonparametric,CFG.perform_parametric) &&  length(CFG.predefined_variance_function2)<3 && isempty(CFG.variance_function_2)
    CFG.compute_variance_function_2=1;
end

%Check that there are sufficient samples to estimate the variance funtions
if CFG.compute_variance_function_1==1
  if and(length(find(CFG.SAMPLES==1))<2,not(CFG.merge_sample1))
        error('Not sufficient samples to estimate variance function for sample 1');
    end
end

if CFG.compute_variance_function_2==1
      if and(length(find(CFG.SAMPLES==2))<2,not(CFG.merge_sample2))
        error('Not sufficient samples to estimate variance function for sample 2');
    end
end



%If a variance function can be loaded use this function and do not
%compute a new one
if or(CFG.perform_nonparametric,CFG.perform_parametric) &&  length(CFG.predefined_variance_function1)<3 && not(isempty(CFG.variance_function_1))
    CFG.compute_variance_function_1=0;
end
    
if or(CFG.perform_nonparametric,CFG.perform_parametric) &&  length(CFG.predefined_variance_function2)<3 && not(isempty(CFG.variance_function_2))
    CFG.compute_variance_function_2=0;
end

   
%If a variance function can be loaded use this function and do not
%compute a new one
if or(CFG.perform_nonparametric,CFG.perform_parametric) &&  length(CFG.predefined_variance_function1)<3 && not(isempty(CFG.variance_function_1))
    CFG.compute_variance_function_1=0;
end
    
if or(CFG.perform_nonparametric,CFG.perform_parametric) &&  length(CFG.predefined_variance_function2)<3 && not(isempty(CFG.variance_function_2))
    CFG.compute_variance_function_2=0;
end
    


