function mTIM_predict(exp_name, time_stamp)
% time_stamp : corresponds to a previous trained CFG.out_train_dir

CFG = general_settings(exp_name, time_stamp, 0,0);
main_run_evaluation(CFG);

