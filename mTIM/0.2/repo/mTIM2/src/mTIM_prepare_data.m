function mTIM_prepare_data(exp_name)

warning('This tool is currently under development (do not use it yet)');

CFG = general_settings(exp_name,[],1,0);
main_generate_data(CFG);

fprintf('Data Generated.\n');

