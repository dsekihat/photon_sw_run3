#!/bin/bash

python run.py -i AnalysisResults_HL_75289.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLooseTrack;
python run.py -i AnalysisResults_HL_75290.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withLKB;
python run.py -i AnalysisResults_HL_75291.root -c configs/config_pp_13.6TeV.yml -t data -s _V0withTPConlyTrack;
