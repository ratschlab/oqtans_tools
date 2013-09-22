function []=rdiff(ARGS)
% rdiff()
%

if isempty(ARGS)
    usage() ;
    exit(-1) ;
end
%ARGS=[ARGS ':.:.'] ;

%%% Add paths  %%%
fprintf('Set the paths\n')
CFG.paths = set_rdiff_paths();


%%% Read configuration file %%%
fprintf('Load configuration\n')
CFG = configure_rdiff(CFG);
CFG = process_command_line_args(CFG,ARGS);
CFG = process_configure_rdiff(CFG);

%%% Get read counts %%%

%load the gene structure
load(CFG.genes_path , 'genes');

% mask the regions which overlap with other genes
fprintf('Compute regions common to multiple genes\n')
[genes]=detect_overlapping_regions(genes);

%Precompute testing regions
fprintf('Compute alternative regions\n')
[genes]=compute_testing_region(CFG,genes);

%Get the gene expression
if CFG.estimate_gene_expression
    fprintf('Measure gene expression\n')
    get_read_counts(CFG,genes);
end

%%% Estimate variance function %%%
if CFG.perform_nonparametric
    variance_function_nonparametric_1=[];
    variance_function_nonparametric_2=[];
    [variance_function_nonparametric_1, variance_function_nonparametric_2]=estimate_variance_nonparametric(CFG,genes);
end

if CFG.perform_poisson
    variance_function_parametric_1=[];
    variance_function_parametric_2=[];
end
if CFG.perform_parametric
    variance_function_parametric_1=[];
    variance_function_parametric_2=[];
    [variance_function_parametric_1, variance_function_parametric_2]=estimate_variance_parametric(CFG,genes); 
end

%If only gene expression is needed, stop here
if CFG.only_gene_expression
    return
end


%%% Perform tests &  Write output %%%

%Run the prametric tests
if or(CFG.perform_parametric,CFG.perform_poisson)
    perform_parametric_tests(CFG,genes,variance_function_parametric_1, variance_function_parametric_2)
end
%Run the nonparametric tests
if or(CFG.perform_nonparametric,CFG.perform_mmd)
    perform_nonparametric_tests(CFG,genes,variance_function_nonparametric_1, variance_function_nonparametric_2)
end

return



