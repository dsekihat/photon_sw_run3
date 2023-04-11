import os
import sys
import argparse
import json
import subprocess
import logging

logger = logging.getLogger('LoggingTest')
logger.setLevel(10)
sh = logging.StreamHandler()
logger.addHandler(sh)
fh = logging.FileHandler('test.log')
logger.addHandler(fh)
formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
fh.setFormatter(formatter)
sh.setFormatter(formatter)

parser = argparse.ArgumentParser(description="argparse for skimming O2 data table")
parser.add_argument("--json", default="config_skimming.json", type=str, help="json file for configuration", required=True)
parser.add_argument('--pcm', action='store_true')
parser.add_argument('--phos', action='store_true')
parser.add_argument('--emc', action='store_true')
args = parser.parse_args()

def run(args):
    print(sys._getframe().f_code.co_name);
    print("skimming PCM :", args.pcm);
    print("skimming PHOS :", args.phos);
    print("skimming EMC :", args.emc);

    config = {};
    with open(args.json) as configFile:
        config = json.load(configFile)
    #print(config);

    #dependencies
    common_dep = [
#        "o2-analysis-collision-converter",
        "o2-analysis-timestamp",
        "o2-analysis-event-selection",
        "o2-analysis-multiplicity-table",
        "o2-analysis-trackselection",
        "o2-analysis-track-propagation",
    ]
    pcm_dep = [
        "o2-analysis-pid-tpc-base",
        "o2-analysis-pid-tpc-full",
        "o2-analysis-em-create-pcm",
        "o2-analysis-em-skimmergammaconversions"
    ]
    phos_dep = [
        "o2-analysis-calo-clusters",
        "o2-analysis-em-skimmer-phos"
    ]

    emc_dep = [
        "o2-analysis-je-emcal-correction-task",
        "o2-analysis-em-skimmergammacalo"
    ]

    emreducedevent_dep = [
        "o2-analysis-em-create-emreduced-event"
    ]

    tables = {
        "EMReducedEvents": {"table": "AOD/EMREDUCEDEVENT/0", "treename": "EMEvents"},
        "V0Photons": {"table": "AOD/V0PHOTON/0", "treename": "V0Photons"},
        "V0Legs": {"table": "AOD/V0LEG/0", "treename": "V0Legs"},
        "V0KFs": {"table": "AOD/V0RECALCANDKF/0", "treename": "V0KFs"},
        "PHOSClusters": {"table": "AOD/PHOSCLUSTERS/0", "treename": "PHOSClusters"},
        "EMCClusters": {"table": "AOD/SKIMEMCCLUSTERS/0", "treename": "EMCClusters"},
    }

    # Tables to be written
    common_tables = ["EMReducedEvents"]
    pcm_tables = ["V0Photons","V0Legs", "V0KFs"]
    phos_tables = ["PHOSClusters"]
    emc_tables = ["EMCClusters"]

    tables_needed = {};
    for it in common_tables:
        tables_needed[it]=1;
    if args.pcm:
        for it in pcm_tables:
            tables_needed[it]=1;
    if args.phos:
        for it in phos_tables:
            tables_needed[it]=1;
    if args.emc:
        for it in emc_tables:
            tables_needed[it]=1;

    writer_config_name = "writer_photon_table.json";
    generate_aod_writer(tables_needed, tables, writer_config_name);

    config["create-emreduced-event"]["processDummy"]="false";
    config["create-emreduced-event"]["process_PCM"]="false";
    config["create-emreduced-event"]["process_PHOS"]="false";
    config["create-emreduced-event"]["process_EMC"]="false";
    config["create-emreduced-event"]["process_PCM_PHOS"]="false";
    config["create-emreduced-event"]["process_PCM_EMC"]="false";
    config["create-emreduced-event"]["process_PHOS_EMC"]="false";
    config["create-emreduced-event"]["process_PCM_PHOS_EMC"]="false";
    if args.pcm and not args.phos and not args.emc:
        config["create-emreduced-event"]["process_PCM"]="true";
    if args.phos and not args.pcm and not args.emc:
        config["create-emreduced-event"]["process_PHOS"]="true";
    if args.emc and not args.pcm and not args.phos:
        config["create-emreduced-event"]["process_PHOS"]="true";
    if args.pcm and args.phos and not args.emc:
        config["create-emreduced-event"]["process_PCM_PHOS"]="true";
    if args.pcm and args.emc and not args.phos:
        config["create-emreduced-event"]["process_PCM_EMC"]="true";
    if args.phos and args.emc and not args.pcm:
        config["create-emreduced-event"]["process_PHOS_EMC"]="true";
    if args.pcm and args.phos and args.emc:
        config["create-emreduced-event"]["process_PCM_PHOS_EMC"]="true";

    # Write the updated configuration file into a temporary file
    updatedConfigFileName = "tmp_config_skimming_photon.json"
    with open(updatedConfigFileName, "w") as outputFile:
        json.dump(config, outputFile, indent = 4)

    command = "";
    for idp in common_dep:
        command += idp + " --configuration json://" + updatedConfigFileName + " -b | ";
    if args.pcm:
        for idp in pcm_dep:
            command += idp + " --configuration json://" + updatedConfigFileName + " -b | ";
    if args.phos:
        for idp in phos_dep:
            command += idp + " --configuration json://" + updatedConfigFileName + " -b | ";
    if args.emc:
        for idp in emc_dep:
            command += idp + " --configuration json://" + updatedConfigFileName + " -b | ";

    if args.pcm or args.phos or args.emc:
        for idp in emreducedevent_dep:
            #print(cdp);
            command += idp + " --configuration json://" + updatedConfigFileName + " -b ";

    command += " --aod-writer-json " + writer_config_name + " --severity error;";
    logger.info(command);
    logger.info(tables_needed.keys());
    subprocess.run(command,shell=True);

#________________________________________________
def generate_aod_writer(tables_needed:dict, tables:dict, writer_config_name):
    writerConfig = {}
    writerConfig["OutputDirector"] = {
        "debugmode": True,
        "resfile": "emAO2D",
        "resfilemode": "RECREATE",
        "ntfmerge": 1,
        "OutputDescriptors": [],
        }

    iTable = 0
    for table in tables_needed.keys():
        writerConfig["OutputDirector"]["OutputDescriptors"].insert(iTable, tables[table])
        iTable += 1

    writerConfigFileName = writer_config_name;
    with open(writerConfigFileName, "w") as writerConfigFile:
        json.dump(writerConfig, writerConfigFile, indent = 4)
    #print(writerConfig)
#________________________________________________
if __name__ == "__main__":
    sys.exit(run(args));
#________________________________________________
#________________________________________________
#________________________________________________
