#!/bin/bash


##for LHC22q
python run.py -i AnalysisResults_HL_75289.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withAnyTrack;
python run.py -i AnalysisResults_HL_75291.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTrack;
#python run.py -i AnalysisResults_HL_75290.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;


##for LHC22f
python run.py -i AnalysisResults_HL_74251.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withAnyTrack;
python run.py -i AnalysisResults_HL_74252.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withTPConlyTrack;
#python run.py -i AnalysisResults_HL_74253.root -c configs/config_pp_13.6TeV_LHC22f.yml -t data -s _V0withLKB;
