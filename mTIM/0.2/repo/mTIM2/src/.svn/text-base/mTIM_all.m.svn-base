function mTIM_train(exp_name)

warning('This tool is currently under development (do not use it yet)');

CFG = general_settings(exp_name,[],1,1);

fh = fopen('out.txt','w');
fprintf(fh,'exp_name is: %s\n',exp_name);


fprintf(fh,'Generate data.\n');
main_generate_data(CFG);

fprintf(fh,'Run training.\n');
main_run_training(CFG);

fprintf(fh,'Run prediction.\n');
main_predict(CFG);


fprintf(fh,'Finished.\n');





