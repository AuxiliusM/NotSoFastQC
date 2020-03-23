from statistics import start_analysis, BLOCK_DIST
from cut import cut_file
import json
import os
import argparse

#shebang for unix sys




ILUMINA_1_5 = "KLMNOPQRSTUVWXYZ[\]^_`abcdefghi"
ILUMINA_1_8 = "!\"#$%&'()*+,-./0123456789:;<=>?@"

CONFIG_PATH = "./fastq_config.json"

def main():

    fastq_config = load_json(CONFIG_PATH)

    fastq_config_exists = len(fastq_config) != 0

    if not fastq_config_exists:
        create_config()
    
    

    fastq_config = load_json(CONFIG_PATH)


    consoleinput = argparse_fastq()
    consoleinput = vars(consoleinput)
    
    file_path = consoleinput["file"]
    should_slice = consoleinput["slice"]
 

    #edit_config_consol(consolinput,fastq_config) | see argparse_fastq why | use doc opt maybe
    

    offset = get_version_offset(file_path) 
    
    print("For more settings and customization please use the generated .json file\n")
    print("Info: Cut Mode 1: Discard whole read when avg phred score under threshold \nCut Mode 2: Sliding Window \nCut Mode 3: Discard whole read when under length requirement\n")
    if int(should_slice) == 1:
        cut_file(f"./{file_path}", offset, fastq_config)
    else:
        start_analysis(f"./{file_path}", fastq_config, offset,  fastq_config["show_gui"])
    
    

def get_version_offset(file_path):
    with open(file_path, "r") as f: 
        for i, line in enumerate(f):
            if (i + 1) % BLOCK_DIST == 0:
                for c in line:
                    if c in ILUMINA_1_5:
                        return 66
                    elif c in ILUMINA_1_8:
                        return 33
                if i > 1000:
                    print("Couldnt determine version!")
                    exit()
    
def create_config(out_path="cut", excel_path="statistics", cut_mode=1, threshold=30, window_size=12, window_startcut=100, desired_len=300, show_gui=True):
    config = {
        
        "out_path": out_path,
        "excel_out": excel_path,
        "cut_mode": cut_mode,
        "threshold": threshold,
        "window_size": window_size,
        "window_startcut": window_startcut,
        "desired_len": desired_len,
        "show_gui": show_gui
    } 

    
    save_json(CONFIG_PATH, config)


def edit_config_consol(consolinput,fastq_config):
    print(consolinput)
    print("test")
   
    for key,value in consolinput.items():
        if key != "file" and value != None:
            fastq_config[key] = value

    save_json(CONFIG_PATH, fastq_config)


def argparse_fastq():

    # maybewriteahelphere
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-f', dest="file", required=True, help="path for file goes here")

    # throws errors that i can't figure out i suspect cross interactions | use gui settings or json for change
    
    # parser.add_argument('-g', dest="gui", required=True)
    # parser.add_argument('-o', dest="out_path")
    # parser.add_argument('-e', dest="exce_lOut")
    # parser.add_argument('-c', dest="cut_mode")
    # parser.add_argument('-t', dest="threshold")
    # parser.add_argument('-w', dest="window_size")
    # parser.add_argument('-s', dest="windowstartcut")

    parser.add_argument('-s', dest="slice", default = 0 , help="should the program cut 1 or analyise 0")
    


    args=parser.parse_args()
        

   
    return args
     

def load_json(path):
    if os.path.isfile(path): 
        f = open(path, "r")
        json_data = json.load(f)
        f.close()
    else:
        json_data = []

    return json_data


def save_json(path, json_data):
    f = open(path, "w")
    json.dump(json_data, f)
    f.close()

class Data:
    
    def __init__(self):
        # not sure yet what to do with classes
        pass


if __name__ == "__main__":
    main()