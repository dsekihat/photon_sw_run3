#!/bin/bash


###for LHC23zzf
#python run_nm.py -i AnalysisResults_HL_165787.root -c configs/config_PbPb_5.36TeV_LHC23zzf_pi0.yml -t data -s _V0withPCB_20240212;

###for LHC23zs
#python run_nm.py -i AnalysisResults_HL_166514.root -t data -c configs/config_pp_13.6TeV_LHC23zs_pi0.yml -s _V0withPCB_20240213;
#python run_nm.py -i AnalysisResults_HL_157194.root -t data -c configs/config_pp_13.6TeV_LHC23zs_pi0.yml -s _V0withPCB_20240124;
#python run_hbt.py -i AnalysisResults_HL_155758.root -t data -c configs/config_pp_13.6TeV_LHC23zs_photon_hbt.yml -s _V0withPCB_20240121;
#python run_nm.py -i AnalysisResults_HL_155758.root -t data -c configs/config_pp_13.6TeV_LHC23zs_pi0.yml -s _V0withPCB_20240121;
#python run_nm.py -i AnalysisResults_HL_153470.root -t data -c configs/config_pp_13.6TeV_LHC23zs_pi0.yml -s _V0withPCB_20240116;
#python run_hbt.py -i AnalysisResults_HL_153470.root -t data -c configs/config_pp_13.6TeV_LHC23zs_photon_hbt.yml -s _V0withPCB_20240116;
#python run_tagging_pi0.py -i AnalysisResults_HL_153470.root -t data -c configs/config_pp_13.6TeV_LHC23zs_tagging_pi0.yml -s _V0withPCB;

###for LHC23zzk_apass1_relval
#python run_nm.py -i AnalysisResults_alien20_new.root -c configs/config_PbPb_5.36TeV_LHC23zzk_pi0.yml -t data -s _V0withPCB_20231204;
#python run_nm.py -i AnalysisResults_HL_148948.root -c configs/config_PbPb_5.36TeV_LHC23zzk_pi0.yml -t data -s _V0withPCB_20240108;

###for LHC23zzh
#python run_nm.py -i AnalysisResults_HL_134949.root -c configs/config_PbPb_5.36TeV_LHC23zzh_pi0.yml -t data -s _V0withPCB_20231114;
#python run_nm.py -i AnalysisResults_HL_133926.root -c configs/config_PbPb_5.36TeV_LHC23zzh_pi0.yml -t data -s _V0withPCB_20231110;
#python run_nm.py -i AnalysisResults_HL_132818.root -c configs/config_PbPb_5.36TeV_LHC23zzh_pi0.yml -t data -s _V0withPCB_20231108;


###for LHC23zo, zp low B
#python run_nm.py -i AnalysisResults_HL_133022.root -c configs/config_pp_13.6TeV_LHC23zo_eta.yml -t data -s _V0withPCB_20231108_1bigbin;
#python run_nm.py -i AnalysisResults_HL_133022.root -c configs/config_pp_13.6TeV_LHC23zo_pi0.yml -t data -s _V0withPCB_20231108_1bigbin;
#python run_nm.py -i AnalysisResults_HL_133022.root -c configs/config_pp_13.6TeV_LHC23zo_pi0.yml -t data -s _V0withPCB_20231108;
#python run_nm.py -i AnalysisResults_HL_141865.root -c configs/config_pp_13.6TeV_LHC23zo_pi0.yml -t data -s _V0withPCB_20231205;
#python run_tagging_pi0.py -i AnalysisResults_HL_133022.root -c configs/config_pp_13.6TeV_LHC23zo_tagging_pi0.yml -t data -s _V0withPCB_20231108;
#python run_tagging_pi0.py -i AnalysisResults_HL_141865.root -c configs/config_pp_13.6TeV_LHC23zo_tagging_pi0.yml -t data -s _V0withPCB_20231205;

###for LHC23zz
#python run_nm.py -i AnalysisResults_HL_125314.root -c configs/config_PbPb_5.36TeV_LHC23zz_pi0.yml -t data -s _V0withAnyTrack_125314;

###for LHC22s
#python run_nm.py -i AnalysisResults_HL_94251.root -c configs/config_PbPb_5.36TeV_LHC22s_pi0.yml -t data -s _V0withAnyTrack_KF_PV_M0;
#python run_nm.py -i AnalysisResults_HL_92433.root -c configs/config_PbPb_5.36TeV_LHC22s_pi0.yml -t data -s _V0withAnyTrack_KF;
#python run_nm.py -i AnalysisResults_HL_89016.root -c configs/config_PbPb_5.36TeV_LHC22s_pi0.yml -t data -s _V0withAnyTrack;

####for LHC22o

python run_nm.py -i AnalysisResults_HL_197579.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _20240413;
python run_nm.py -i AnalysisResults_HL_197579.root -c configs/config_pp_13.6TeV_LHC22o_eta.yml -t data -s _20240413;
python run_nm.py -i AnalysisResults_HL_197593.root -c configs/config_pp_13.6TeV_LHC24b1_pi0.yml -t mc -s _20240413;
python run_nm.py -i AnalysisResults_HL_197593.root -c configs/config_pp_13.6TeV_LHC24b1_eta.yml -t mc -s _20240413;
python run_tagging_pi0.py -i AnalysisResults_HL_197579.root -t data -c configs/config_pp_13.6TeV_LHC22o_tagging_pi0.yml -s _20240413;
python run_tagging_pi0.py -i AnalysisResults_HL_197593.root -t mc -c configs/config_pp_13.6TeV_LHC24b1_tagging_pi0.yml -s _20240413;

#python run_nm.py -i AnalysisResults_HL_195791.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _20240410;
#python run_nm.py -i AnalysisResults_HL_195791.root -c configs/config_pp_13.6TeV_LHC22o_eta.yml -t data -s _20240410;
#python run_nm.py -i AnalysisResults_HL_196011.root -c configs/config_pp_13.6TeV_LHC24b1_eta.yml -t mc -s _20240410_TPCPID;
#python run_nm.py -i AnalysisResults_HL_196011.root -c configs/config_pp_13.6TeV_LHC24b1_pi0.yml -t mc -s _20240410_TPCPID;
#python run_nm.py -i AnalysisResults_HL_195792.root -c configs/config_pp_13.6TeV_LHC24b1_eta.yml -t mc -s _20240410_TPCPIDNN;
#python run_nm.py -i AnalysisResults_HL_195792.root -c configs/config_pp_13.6TeV_LHC24b1_pi0.yml -t mc -s _20240410_TPCPIDNN;
#python run_tagging_pi0.py -i AnalysisResults_HL_195791.root -t data -c configs/config_pp_13.6TeV_LHC22o_tagging_pi0.yml -s _20240410;
#python run_tagging_pi0.py -i AnalysisResults_HL_195792.root -t mc -c configs/config_pp_13.6TeV_LHC24b1_tagging_pi0.yml -s _20240410_TPCPIDNN;
#python run_tagging_pi0.py -i AnalysisResults_HL_196011.root -t mc -c configs/config_pp_13.6TeV_LHC24b1_tagging_pi0.yml -s _20240410_TPCPID;

#python run_tagging_pi0.py -i AnalysisResults_HL_191046.root -t mc -c configs/config_pp_13.6TeV_LHC24b1_tagging_pi0.yml -s _20240401;
#python run_tagging_pi0.py -i AnalysisResults_HL_190549.root -t data -c configs/config_pp_13.6TeV_LHC22o_tagging_pi0.yml -s _20240401;
#python run_nm.py -i AnalysisResults_HL_191046.root -c configs/config_pp_13.6TeV_LHC24b1_eta.yml -t mc -s _20240401;
#python run_nm.py -i AnalysisResults_HL_190549.root -c configs/config_pp_13.6TeV_LHC22o_eta.yml -t data -s _20240401;
#python run_nm.py -i AnalysisResults_HL_191046.root -c configs/config_pp_13.6TeV_LHC24b1_pi0.yml -t mc -s _20240401;
#python run_nm.py -i AnalysisResults_HL_190549.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _20240401;
#python run_hbt.py -i AnalysisResults_HL_153950.root -t data -c configs/config_pp_13.6TeV_LHC22o_photon_hbt.yml -s _V0withPCB_20240117;
#python run_tagging_pi0.py -i AnalysisResults_HL_153950.root -t data -c configs/config_pp_13.6TeV_LHC22o_tagging_pi0.yml -s _V0withPCB_20240117;
#python run_nm.py -i AnalysisResults_HL_141697.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withPCB_20231205;
#python run_nm.py -i AnalysisResults_HL_133366.root -c configs/config_pp_13.6TeV_LHC22o_eta.yml -t data -s _V0withPCB_20231108_1bigbin;
#python run_nm.py -i AnalysisResults_HL_133366.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withPCB_20231108_1bigbin;
#python run_nm.py -i AnalysisResults_HL_133366.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withPCB_20231108;
#python run_nm.py -i AnalysisResults_HL_132808.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withPCB_20231108;
#python run_nm.py -i AnalysisResults_HL_117656.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withAnyTrack_20230828_nsw1;
#python run_nm.py -i AnalysisResults_HL_105667.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withAnyTrack_nsw5_20230721;
#python run_nm.py -i AnalysisResults_HL_105511.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withAnyTrack_nsw1_20230721;
#python run_nm.py -i AnalysisResults_HL_103430.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withAnyTrack_nsw5_20230718;
#python run_nm.py -i AnalysisResults_HL_88823.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withLKB;
#python run_nm.py -i AnalysisResults_HL_88830.root -c configs/config_pp_13.6TeV_LHC22o_pi0.yml -t data -s _V0withAnyTrack;


##for LHC22m
#python run_nm.py -i AnalysisResults_HL_117657.root -c configs/config_pp_13.6TeV_LHC22m_pi0.yml -t data -s _V0withAnyTrack_20230828_nsw1;
#python run_nm.py -i AnalysisResults_HL_105666.root -c configs/config_pp_13.6TeV_LHC22m_pi0.yml -t data -s _V0withAnyTrack_nsw5_20230721;
#python run_nm.py -i AnalysisResults_HL_105513.root -c configs/config_pp_13.6TeV_LHC22m_pi0.yml -t data -s _V0withAnyTrack_nsw1_20230721;
#python run_nm.py -i AnalysisResults_HL_103429.root -c configs/config_pp_13.6TeV_LHC22m_pi0.yml -t data -s _V0withAnyTrack_nsw5_20230718;
#python run_nm.py -i AnalysisResults_HL_96456.root -c configs/config_pp_13.6TeV_LHC22m_pi0.yml -t data -s _V0withAnyTrack_nsw1_20230718;
#python run_nm.py -i AnalysisResults_nsw1_nopv_ap_rcut.root -c configs/config_pp_13.6TeV_test_nsw.yml -t data -s _V0withAnyTrack_test_nsw1;
#python run_nm.py -i AnalysisResults_nsw5_nopv_ap_rcut.root -c configs/config_pp_13.6TeV_test_nsw.yml -t data -s _V0withAnyTrack_test_nsw5;

#python run.py -i AnalysisResults_HL_78448.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTracksBar;
#python run.py -i AnalysisResults_HL_78449.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTracks;
#python run.py -i AnalysisResults_HL_78450.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;
#python run.py -i AnalysisResults_HL_78451.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withAnyTracks;
#python run.py -i AnalysisResults_HL_78314.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB_msas;

####for LHC22q
#python run_nm.py -i AnalysisResults_HL_117651.root -c configs/config_pp_13.6TeV_LHC22q_eta.yml -t data -s _V0withAnyTrack_20230828_nsw1;
#python run_nm.py -i AnalysisResults_HL_117651.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack_20230828_nsw1;
#python run_tagging_pi0.py -i AnalysisResults_HL_111542.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack_20230806_nsw1;
#python run_nm.py -i AnalysisResults_HL_111542.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack_20230806_nsw1;
#python run_nm.py -i AnalysisResults_HL_111542.root -c configs/config_pp_13.6TeV_LHC22q_eta.yml -t data -s _V0withAnyTrack_20230806_nsw1;
#python run_single_photon.py -i AnalysisResults_HL_106685.root -c configs/config_pp_13.6TeV_LHC22q_photon.yml -t data -s _V0withAnyTrack_20230725_nsw1;
#python run_tagging_pi0.py -i AnalysisResults_HL_106685.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack_20230725_nsw1;
#python run_tagging_pi0.py -i AnalysisResults_HL_103791.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack_20230721_nsw1;
#python run_tagging_pi0.py -i AnalysisResults_HL_103290.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack_20230718_nsw1;
#python run_nm.py -i AnalysisResults_HL_103290.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack_20230718_nsw1;
#python run_tagging_pi0.py -i AnalysisResults_HL_10313.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack_20230711;
#python run_nm.py -i AnalysisResults_HL_10313.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack_20230711;
#python run_nm.py -i AnalysisResults_HL_98868.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack_LHC22q_pass4_lowIR;
#python run_nm.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack;
#python run_nm.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_LHC22q_eta.yml -t data -s _V0withAnyTrack;
#python run_tagging_pi0.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack;
#python run_single_photon.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_LHC22q_photon.yml -t data -s _V0withAnyTrack;
#python run_hbt.py -i AnalysisResults_hbt_test_all.root -t data -c configs/config_pp_13.6TeV_LHC22q_photon_hbt.yml -s _V0withAnyTrack;
#python run_tagging_pi0.py -i AnalysisResults_HL_87321.root -c configs/config_pp_13.6TeV_LHC22q_tagging_pi0.yml -t data -s _V0withAnyTrack;
#python run_nm.py -i AnalysisResults_HL_87321.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTrack_onlyPHOS;
#python run_nm.py -i AnalysisResults_HL_87321.root -c configs/config_pp_13.6TeV_LHC22q_eta.yml -t data -s _V0withAnyTrack_onlyPHOS;

##python run.py -i AnalysisResults_HL_75289.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withAnyTrack;
##python run.py -i AnalysisResults_HL_75291.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTrack;
#python run.py -i AnalysisResults_HL_75290.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;
#python run.py -i AnalysisResults_HL_79577.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;

#python run_nm.py -i AnalysisResults_HL_82428.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withAnyTracks_LHC22q;
#python run_nm.py -i AnalysisResults_HL_82428.root -c configs/config_pp_13.6TeV_LHC22q_eta.yml -t data -s _V0withAnyTracks_LHC22q;
#python run_nm.py -i AnalysisResults_HL_82427.root -c configs/config_pp_13.6TeV_LHC22q_pi0.yml -t data -s _V0withLKB;
#python run_nm.py -i AnalysisResults_HL_82427.root -c configs/config_pp_13.6TeV_LHC22q_eta.yml -t data -s _V0withLKB;

##
####for LHC22f
#python run_nm.py -i AnalysisResults_HL_124838.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_20231006;
#python run_nm.py -i AnalysisResults_HL_103289.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_20230718_nsw1;
#python scan_cb_param.py -i AnalysisResults_HL_94252.root -c configs/config_pp_13.6TeV_LHC22f_pi0_cb.yml -t data -s _V0withAnyTrack_KF_PV_M0;
#python run_nm.py -i AnalysisResults_HL_97028.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_KF_PV_M0_PID3_ITSonly;
#python run_nm.py -i AnalysisResults_HL_96457.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_KF_PV_M0_PID3;
#python run_nm.py -i AnalysisResults_HL_94252.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_KF_PV_M0;
#python run_nm.py -i AnalysisResults_HL_94252.root -c configs/config_pp_13.6TeV_LHC22f_eta.yml -t data -s _V0withAnyTrack_KF_PV_M0;
#python run_nm.py -i AnalysisResults_HL_92593.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_KF_X83;
#python run_nm.py -i AnalysisResults_HL_92432.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_KF;
#python run_nm.py -i AnalysisResults_HL_88822.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withLKB;
#python run_nm.py -i AnalysisResults_HL_88933.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack;
#python run_nm.py -i AnalysisResults_HL_88933.root -c configs/config_pp_13.6TeV_LHC22f_eta.yml -t data -s _V0withAnyTrack;

##python run.py -i AnalysisResults_HL_74251.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withAnyTrack;
##python run.py -i AnalysisResults_HL_74252.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withTPConlyTrack;
###python run.py -i AnalysisResults_HL_74253.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withLKB;
#python run_nm.py -i AnalysisResults_HL_85104.root -c configs/config_pp_13.6TeV_LHC22f_pi0.yml -t data -s _V0withAnyTrack_LHC22f;


##### for MC
#python run_nm.py -i AnalysisResults_HL_169607.root -c configs/config_PbPb_5.36TeV_LHC23k6_pi0.yml -t mc -s _V0withPCM_20240221;
#python run_nm.py -i AnalysisResults_HL_117572.root -c configs/config_pp_13.6TeV_pi0_mc.yml -t mc -s _V0withAnyTrack_20230828_nsw1;
#python run_single_photon.py -i AnalysisResults_LHC23d10_10kHz.root -c configs/config_pp_13.6TeV_photon_mc.yml -t mc -s _V0withAnyTracks_10kHz;
#python run_nm.py -i AnalysisResults_LHC23d10_10kHz_20230606.root -c configs/config_pp_13.6TeV_eta_mc.yml -t mc -s _V0withAnyTrack_eta_xsection;
#python run_nm.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_eta_mc.yml -t data -s _V0withAnyTrack_eta_xsection;
#python run_nm.py -i AnalysisResults_mc_pi0_20230606.root -c configs/config_pp_13.6TeV_pi0_mc.yml -t mc -s _V0withAnyTrack_xsection;
#python run_nm.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_pi0_mc.yml -t data -s _V0withAnyTrack_xsection;
#python run_single_photon.py -i AnalysisResults_HL_85687.root -c configs/config_pp_13.6TeV_photon_mc.yml -t data -s _V0withAnyTrack_xsection;
#python run_single_photon.py -i AnalysisResults_LHC23d10_10kHz_20230606.root -c configs/config_pp_13.6TeV_photon_mc.yml -t mc -s _V0withAnyTrack_xsection;
#python run_single_photon.py -i AnalysisResults_HL_132941.root -c configs/config_pp_13.6TeV_LHC23d1k_photon.yml -t mc -s _V0withPCB_20231108;
