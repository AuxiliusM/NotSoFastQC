import numpy as np
import matplotlib.pyplot as plt
import time
import math
from appJar import gui
import json
import pandas as pd
import datetime
from cut import cut_file


GC_ACCURACY = 100
BLOCK_DIST = 4
CONFIG_PATH = "./fastq_config.json"


def start_analysis(file_path, fastq_config, offset, GUI=True): 
    if GUI:
        print(file_path)
        loading_meter_primegui(file_path, fastq_config, offset)
       
        #loading_meter()

    else:
        print(file_path)
        analyze_file(file_path, fastq_config, offset, GUI)
        pass
                
def analyze_file(file_path, fastq_config, offset, GUI, app=None):
    total_scores = []
    total_count = []

    length_count = [] 

    max_read_length = 0
    min_read_length = math.inf

    distribution_GC = [0] * GC_ACCURACY

    poor_quality = 0

    version_name = get_version_name(offset)

    print("Analysis started")
    print(version_name)
    print("")
    
    file_line_count = file_len(file_path)
    

    with open(file_path, "r") as f: 

        for i, line in enumerate(f):
            if (i + 1) % BLOCK_DIST == 0:
                scores = get_scores(line, offset)


                score_average = np.average(scores)
                poor_quality += 1 if score_average <= fastq_config["threshold"] else 0
                
                l = len(scores)
                
                if l > max_read_length:
                    max_read_length = l

                if l < min_read_length:
                    min_read_length = l

                expand_list(total_scores, l)
                expand_list(total_count, l)
                expand_list(length_count, l + 1)

                count_scores(total_scores, total_count, scores)
                evaluate_length(scores, length_count)


            elif (i + 3) % BLOCK_DIST == 0:
                count_GC(line, distribution_GC)

            if i/BLOCK_DIST % 10000 == 0:
                percent_done = ((i/BLOCK_DIST*4) / file_line_count)*100
                print(f"I am working on line {i//BLOCK_DIST*4} ({round(percent_done)}%) atm")
                
                if GUI:
                    #update direct from other thread uff / change at some point!
                    app.setMeter("progress", percent_done, text=None)
                    
    avg = np.array(avg_scores(total_scores, total_count))

    if GUI:
        app.setMeter("progress", 100, text=None)
    
    print("Analysis done (100%)")
    
    if not GUI:    
        matplotlib_plots(distribution_GC,avg,min_read_length,max_read_length,length_count,version_name,poor_quality) 
        plt.show()

    else:
        create_main_gui(file_path,fastq_config,offset,distribution_GC,avg,min_read_length,max_read_length,length_count,version_name,poor_quality,app)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1   


def expand_list(li, min_length):
    if len(li) < min_length:
        diff = min_length - len(li)
        li.extend([0] * diff)


def get_scores(chars, offset):
    scores = []
    for c in chars:
        if c != "\n":
            scores.append(ord(c) - offset)

    return scores


def count_scores(total_scores, total_count, new_scores):
    for i in range(0, len(new_scores)):
        total_scores[i] += new_scores[i]
        total_count[i] += 1


def avg_scores(total_scores, total_count):
    avg = [0] * len(total_scores)
    
    for i in range(0, len(total_scores)):
        avg[i] = total_scores[i] / total_count[i] if total_count[i] > 0 else 0

        
    return avg



def count_GC(chars, distribution_GC):
    chars = chars.upper()
    amount_base = (len(chars)) - 1
    
    ammount_gc = chars.count('C') + chars.count('G')
    percent_gc = ammount_gc / amount_base
    
    index = round(percent_gc * GC_ACCURACY)
    distribution_GC[index] += 1


def evaluate_length(scores, length_count):
    length_count[len(scores)] += 1


def std_avg(ev_array):

    total_elements = np.sum(ev_array)

 

    ev_array_mult = ev_array.copy()

    for i, num in enumerate(ev_array): 
        ev_array_mult[i] =  num * i

    avg = np.sum(ev_array_mult) / total_elements

    numerator = 0
    for i, num in enumerate(ev_array):
        numerator += ((i - avg) ** 2) * num

    variance = numerator / total_elements    
    
    return {"avg": avg, "var": variance, "std": variance**(1/2)} 


def matplotlib_plots(distribution_GC,avg,min_read_length,max_read_length,length_count,version_name,poor_quality):
    
    gcplot = plt.figure()
    phredplot = plt.figure()
    lengthplot = plt.figure()


    ax1 = gcplot.add_subplot()
    ax1.plot(distribution_GC)
    ax1.set_title("GC distribution over all sequences")
    ax1.set_xlabel(f"Mean GC content in {100/GC_ACCURACY} - % steps")

    returndic = std_avg(distribution_GC)
    print(f"Statistics: Mean GC Content{returndic}")


    ax2 = phredplot.add_subplot()
    ax2.plot(avg)
    ax2.set_title(f"Quality scores across all bases ({version_name} encoding)")
    ax2.set_xlabel("Position in read (bp)")

    returndic["avg"] = np.average(avg) 
    returndic["std"]= np.std(avg) 
    
    print(f"Statistics: Quality Scores{returndic}")


    
    ax3 = lengthplot.add_subplot()
    ax3.plot(range(min_read_length, max_read_length), length_count[min_read_length:max_read_length])
    ax3.set_title(f"Distribution of sequence lengths over all sequences")
    ax3.set_xlabel("Sequences Length (bp)")

    

    returndic = std_avg(length_count)
    print(f"Statistics: Sequence Lengths{returndic}")


    print(f"'Poor quality' sequences: {poor_quality}")

def multithreading():

    # Comming in ver 2.0 / is it even viable for this task? | (hard)
    # idea: multipe files / same content 
    # for first try: 2 cores, one starting at start of file, other at end both going to middel
    # then bring results together
    # text could be split into available_cores parts | bandwidth limit and threading problems ahead

    pass

def ncbi_blast():

    # If i have time for it with biopython
    # rip because of argparse
    # Bio.Blast.NCBIWWW module , API Stuff here | return organism
    pass
    
def loading_meter_primegui(file_path, fastq_config, offset):
    GUI = True
    app = gui("Results","500x250")

    app.addLabel("plswait", "Analysing please wait")
    app.getLabelWidget("plswait").config(font="Arial 30 bold underline")
    
    app.addMeter("progress")
    app.setMeterFill("progress", "green")
    app.thread(analyze_file, file_path, fastq_config, offset, GUI, app) 

    app.go()

    
def create_main_gui(file_path,fastq_config,offset,distribution_GC,avg,min_read_length,max_read_length,length_count,version_name,poor_quality,app):

    

    def press_decider(button):
        if button == "Plots & Ø,σ,M":
            
            show_stats_window()

            matplotlib_plots(distribution_GC,avg,min_read_length,max_read_length,length_count,version_name,poor_quality)
            plt.show()

        elif button == "Settings":
            config_click()
        
        elif button == "Excel Export":

            date_time = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")

            oldname = file_path.split(".")
            removed_sl = oldname[1].split("/")

            excel_out_string = f"{fastq_config['excel_out']}_{removed_sl[1]} {date_time} {version_name}.xlsx"

            excel_export(excel_out_string,distribution_GC,avg,min_read_length,max_read_length,length_count)

        elif button == "Cut Fastq":
            cut_click()

    def show_stats_window():

        #should make function out of this
        returndic = {}
        returndic["avg"] = np.average(avg) 
        returndic["std"]= np.std(avg) 
        stats1 = f"Quality Scores: {returndic}"

        returndic = std_avg(length_count)
        stats2 = f"Sequence Lengths: {returndic}"
            

        returndic = std_avg(distribution_GC)
        stats3 = f"Mean GC Content: {returndic}"
        #----------------------------------

        app.setLabel("stat1", stats1)
        app.setLabel("stat2", stats2)
        app.setLabel("stat3", stats3)

        

        app.showSubWindow("stat_window")

    def config_click():
            app.setEntry("Out Prefix", fastq_config["out_path"])
            app.setEntry("Excel Prefix", fastq_config["excel_out"])
            app.setEntry("Cut Mode", fastq_config["cut_mode"])
            app.setEntry("Quality Threshold", fastq_config["threshold"])
            app.setEntry("Sliding Window Size", fastq_config["window_size"])
            app.setEntry("Sliding Window Start", fastq_config["window_startcut"])
            app.setEntry("Desired Read Length", fastq_config["desired_len"])

            app.showSubWindow("config_window")
            app.config_window_shown = True



    def config_save_click():

        
        fastq_config["out_path"] = app.getEntry("Out Prefix")
        fastq_config["excel_out"] = app.getEntry("Excel Prefix")
        fastq_config["cut_mode"] = int(app.getEntry("Cut Mode"))
        fastq_config["threshold"] = int(app.getEntry("Quality Threshold"))
        fastq_config["window_size"] = int(app.getEntry("Sliding Window Size"))
        fastq_config["window_startcut"] = int(app.getEntry("Sliding Window Start"))
        fastq_config["desired_len"] = int(app.getEntry("Desired Read Length"))
        
        f = open(CONFIG_PATH, "w")
        json.dump(fastq_config, f)
        f.close()

        app.hideSubWindow("config_window", useStopFunction=False)

    def cut_click():
        
        app.setEntry("File Path", file_path)
        

        app.showSubWindow("cut_window")
        app.config_window_shown = True

    def explorer_win():
        box_entry = app.openBox(title=None, dirName=None, fileTypes=None, asFile=False, parent="cut_window", multiple=False, mode='r')
        app.setEntry("File Path", str(box_entry))

    def cut_file_gui():
        to_cut = app.getEntry("File Path")
        
        
        try:

            # no loading bar for the cut_function, this is already convoluted enough and i wont implement more multithreading for such a small thing

            print("This part is not threaded / GUI will not respond until cut is finished")
            cut_file(to_cut, offset, fastq_config)

        except Exception as e:
            print("File does not follow standard fastq format")
            print(e)
            app.bell()

        

    
    
    
    app.setFont(16)
    app.removeAllWidgets(True)



    app.addLabel("title", "Results:")
    
    
    app.getLabelWidget("title").config(font="Arial 30 bold")

   

    app.addLabel("f1", f"Filename: {file_path}")
    app.addLabel("f2", f"Encoding: {version_name}")
    app.addLabel("f3", f"Total sequences {np.sum(distribution_GC)}")
    app.addLabel("f4", f"'Poor quality' sequences: {poor_quality}")
    app.addLabel("f5", f"Sequence lenght: {min_read_length}-{max_read_length} (bp)")

    returndic = std_avg(distribution_GC)
    
    app.addLabel("f6", f"%GC: {round(returndic['avg'],2)}")
    
    
    app.addButtons(["Plots & Ø,σ,M", "Excel Export", "Cut Fastq", "Settings"], press_decider)

    #sub window cut
    app.startSubWindow("cut_window", title="Cut Window")
    app.addLabelEntry("File Path")
    app.addButton("Select file", explorer_win)
    app.addButton("Cut", cut_file_gui)
    
    
    app.config_window_shown = False
    app.setOnTop(stay=True)
    app.stopSubWindow()

    #sub window settings
    app.startSubWindow("config_window", title="Config Window")
    app.addLabelEntry("Out Prefix")
    app.addLabelEntry("Excel Prefix")
    app.addLabelEntry("Cut Mode")
    app.addLabelEntry("Quality Threshold")
    app.addLabelEntry("Sliding Window Size")
    app.addLabelEntry("Sliding Window Start")
    app.addLabelEntry("Desired Read Length")
    app.addButton("Save", config_save_click)
    app.config_window_shown = False
    app.setOnTop(stay=True)

    #sub window stats
    app.startSubWindow("stat_window", title="Stat Window")

    app.addLabel("detailed", "Detailed Statistics:")
    app.getLabelWidget("detailed").config(font="Arial 24 bold")
    
    app.addLabel("stat1", "")
    app.addLabel("stat2", "")
    app.addLabel("stat3", "")
    
    
    app.config_window_shown = False
    app.setOnTop(stay=True)
    app.stopSubWindow()
    

    app.stopSubWindow()

    
    




def get_version_name(offset):
    
    if offset >= 64:
        version_name = "Illumina 1.5"
    else:
        version_name = "Illumina 1.8"

    return version_name


def excel_export(excel_out_string,distribution_GC,avg,min_read_length,max_read_length,length_count):

    returndic = {}



    df1 = pd.DataFrame(avg)
    df2 = pd.DataFrame(length_count)
    df3 = pd.DataFrame(distribution_GC)


    #should make function out of this
    returndic["avg"] = np.average(avg) 
    returndic["std"]= np.std(avg) 
    stats1 = f"Statistics: Quality Scores{returndic}"

    returndic = std_avg(length_count)
    stats2 = f"Statistics: Sequence Lengths{returndic}"
        

    returndic = std_avg(distribution_GC)
    stats3 = f"Statistics: Mean GC Content{returndic}"
    #----------------------------------

    
    writer = pd.ExcelWriter(excel_out_string, engine='xlsxwriter')



    
    df1.to_excel(writer, sheet_name='Phred q')
    df2.to_excel(writer, sheet_name='Read l')
    df3.to_excel(writer, sheet_name='GC d')
    

    workbook = writer.book

    bold = workbook.add_format({'bold': True})


    worksheet = writer.sheets['Phred q']

    worksheet.write(38, 3, stats1, bold)

    chart = workbook.add_chart({'type': 'line'})

    chart.add_series({
        'categories': ['Phred q', 1, 0, max_read_length, 0],
        'values':     ['Phred q', 1, 1, max_read_length, 1],
    })

    chart.set_x_axis({'name': 'Position in read (bp)', 'position_axis': 'on_tick'})
    chart.set_y_axis({'name': 'Phred quality score', 'major_gridlines': {'visible': False}})
    chart.set_size({'width': 1280, 'height': 720})

    chart.set_legend({'position': 'none'})

    worksheet.insert_chart('D2', chart)




    worksheet = writer.sheets['Read l']

    worksheet.write(38, 3, stats2, bold)

    chart = workbook.add_chart({'type': 'line'})

    chart.add_series({
        'categories': ['Read l', 1, 0, max_read_length, 0],
        'values':     ['Read l', 1, 1, max_read_length, 1],
        'line':   {'color': 'orange'}
    })

    chart.set_x_axis({'name': 'Sequences Length (bp)', 'position_axis': 'on_tick'})
    chart.set_y_axis({'name': 'Ammount of reads', 'major_gridlines': {'visible': False}})
    chart.set_size({'width': 1280, 'height': 720})

    chart.set_legend({'position': 'none'})

    worksheet.insert_chart('D2', chart)




    worksheet = writer.sheets['GC d']

    worksheet.write(38, 3, stats3, bold)

    chart = workbook.add_chart({'type': 'line'})

    chart.add_series({
        'categories': ['GC d', 1, 0, 100, 0],
        'values':     ['GC d', 1, 1, 100, 1],
        'line':   {'color': 'green'}
    })

    chart.set_x_axis({'name': 'GC Content (%)', 'position_axis': 'on_tick'})
    chart.set_y_axis({'name': 'Ammount of reads', 'major_gridlines': {'visible': False}})
    chart.set_size({'width': 1280, 'height': 720})

    chart.set_legend({'position': 'none'})

    worksheet.insert_chart('D2', chart)


    
    writer.save()
 

        

        



        
        