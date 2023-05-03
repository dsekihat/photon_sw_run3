#!/bin/bash


##for LHC22m
python run.py -i AnalysisResults_HL_78448.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTracksBar;
python run.py -i AnalysisResults_HL_78449.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTracks;
#python run.py -i AnalysisResults_HL_78450.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;
python run.py -i AnalysisResults_HL_78451.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withAnyTracks;


####for LHC22q
##python run.py -i AnalysisResults_HL_75289.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withAnyTrack;
##python run.py -i AnalysisResults_HL_75291.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTrack;
###python run.py -i AnalysisResults_HL_75290.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;
##
##
####for LHC22f
##python run.py -i AnalysisResults_HL_74251.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withAnyTrack;
##python run.py -i AnalysisResults_HL_74252.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withTPConlyTrack;
###python run.py -i AnalysisResults_HL_74253.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withLKB;
