import datetime
import time
import re
import os


BLOCK_DIST = 4

def cut_file(file_path, offset, fastq_config): 

    print(file_path)
    print(fastq_config)

    date_time = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")

    #cut_newnewtestfile 2020_02_20_035843.fastq

    file_name = os.path.basename(file_path)
    file_name = file_name.split(".")[0]

    res = re.findall(fastq_config["out_path"] + "_.+\s[0-9]{4}_[0-9]{2}_[0-9]{2}_[0-9]{6}", file_name)
    is_match = len(res) > 0 and res[0] == file_name

    if is_match:
        file_name = re.sub(fastq_config["out_path"] + "_", "", file_name)
        file_name = re.sub(" [0-9]{4}_[0-9]{2}_[0-9]{2}_[0-9]{6}", "", file_name)


    # if ":" in file_path:
    #     oldname = file_path.split("/")
    #     removed_some = oldname[-1].split("/")
    #     print(removed_some)
    #     removed_all_packed = removed_some[0].split(".")
    #     removed_all = removed_all_packed[0]
    #     print(removed_all_packed)
    #     print(removed_all)
        
        

    # else:
    #     oldname = file_path.split(".")
    #     removed_sl = oldname[1].split("/")
    #     removed_all = removed_sl[1]

    print("Cut started")
    file_line_count = file_len(file_path)
    

    with open(f"{fastq_config['out_path']}_{file_name} {date_time}.fastq", "w") as out:
        with open(file_path, "r") as f:
            block = [""] * 4
            identifier = ""
            written_count = 0


            for i, line in enumerate(f):
                index = i % 4

                if i == 0:
                    identifier = line.split(".")[0]
                    identifier = identifier.split("@")[1]
                    

                if index == 0 and i > 0:
                    
                    
                    if fastq_config["cut_mode"] == 1:
                        
                        
                        discarded = quality_discard_block(block, fastq_config["threshold"], offset)
                        if not discarded:

                            written_count += 1
                            s = f"{identifier}.{written_count} {written_count} lenght={len(block[1])-1}\n"
                            block[0] = "@" + s
                            block[2] = "+" + s

                            out.writelines(block)
                    elif fastq_config["cut_mode"] == 2:

                        written_count += 1

                        sliding_window(block, fastq_config["window_size"], fastq_config["threshold"], offset, fastq_config["window_startcut"], identifier)

                        s = f"{identifier}.{written_count} {written_count} lenght={len(block[1])-1}\n"
                        block[0] = "@" + s
                        block[2] = "+" + s

                        out.writelines(block)

                    elif fastq_config["cut_mode"] == 3:
                        
                        if len(block[3]) >= fastq_config["desired_len"]:

                            written_count += 1
                            s = f"{identifier}.{written_count} {written_count} lenght={len(block[1])-1}\n"
                            block[0] = "@" + s
                            block[2] = "+" + s

                            out.writelines(block)

                    else:
                        pass

                if i/BLOCK_DIST % 10000 == 0:
                    percent_done = ((i/BLOCK_DIST*4) / file_line_count)*100
                    print(f"I am working on line {i//BLOCK_DIST*4} ({round(percent_done)}%) atm")

                block[index] = line

    print("Cut done (100%)")

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1  

def sliding_window(block, size, threshold, offset, startcut, identifier):
    window_sum = 0
    

    for i, char in enumerate(block[3]):

        if i < startcut:
            continue

        elif i - startcut < size:
            window_sum += ord(char) - offset
            
            continue
        
        elif window_sum / size < threshold:


            block[1] = block[1][:i - 1] + "\n"
            block[3] = block[3][:i - 1] + "\n"
            break
        
        window_sum += ord(char) - ord(block[3][i - size - 1])


def quality_discard_block(block, threshold, offset):
    amount_scores = (len(block[3])) - 1
    total_line_score = 0

    for c in block[3]:
        total_line_score += ord(c) - offset

    total_line_score = total_line_score / amount_scores

    return total_line_score <= threshold





    

