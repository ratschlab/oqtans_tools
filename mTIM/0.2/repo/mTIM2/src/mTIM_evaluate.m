function mTIM_predict(exp_name, time_stamp)
% time_stamp : corresponds to a previous trained CFG.out_train_dir

warning('This tool is currently under development (do not use it yet)');

CFG = general_settings(exp_name, time_stamp, 0,0);
main_run_evaluation(CFG);

