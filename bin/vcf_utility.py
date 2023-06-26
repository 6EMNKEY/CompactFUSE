import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored

def is_float(string):
    if string.replace(".", "").isnumeric():
        return True
    else:
        return False

def read_data(data):
    x = sub.run("grep '##' "+data+" | wc -l", shell=True, capture_output=True)
    df = pd.read_table(data,index_col=None ,skiprows=int(x.stdout.decode())+1, header=None)
    lines = sub.run("grep '##' "+data, shell=True, capture_output=True)
    return df, lines.stdout.decode().split("\n")
def get_metadata_on_vcf(lines):
    #!<hdsu> : Description
    alter = {}
    #!<id> : (type, description, number)
    info = {}
    #!<id> : Description
    formatd = {}
    for line in lines:
        if "##ALT" in line:
            firstspl= line.strip(">").split("<")[1].split(",")
            item= []
            for i in firstspl:
                i = i.split("=")
                item.append(i[0])
                item.append(i[1])
            alter["<"+item[1]+">"] = item[3].strip('"')

        elif "##INFO" in line:
            firstspl= line.strip(">").split("<")[1].split(",")
            item= []
            if len(firstspl) > 4:
                tmp = ",".join(firstspl[3:])
                firstspl = firstspl[0:4]
                firstspl[3] = tmp
            
            for i in firstspl:
                i = i.split("=")
                item.append(i[0])
                item.append(i[1])

            info[item[1]]= [item[7], item[5], item[3]]
        elif "##FORMAT" in line:
            firstspl= line.strip(">").split("<")[1].split(",")
            if len(firstspl) > 4:
                tmp = ",".join(firstspl[3:])
                firstspl = firstspl[0:4]
                firstspl[3] = tmp
            item= []
            for i in firstspl:
                i = i.split("=")
                item.append(i[0])
                item.append(i[1])
            formatd[item[1]]= [item[7]]     
    return alter,info,formatd
    
def write_info_part(info, information):
    #!THIS SHOULD BE THE FINAL FORMAT:
    #! CHROM INIPOS DESCRIPTION-FULL
    #! fast_search fast_search Already_formated("\n")
    # print(information)
    # print(info)
    numerals = {}
    text = "INFORMATION: \n"
    flags = "FLAGS: "
    for item in info.split(";"):
        #print(item)
        if "=" in item:
            items = item.split("=")
            if items[0] in information.keys():  
                text += information[items[0]][0] +": " +items[1] + "\n"
            else:
                text += items[0] +": " +items[1] + "\n"
        else:
            flags+=item +"; "
    text+= flags +"\n\n"

    return text


def write_format_part(formata, formatc):
    text="FORMAT: \n\n"
    flags = []
    for i in list(formata)[0].split(":"):
        flags.append([i])
    for test in list(formata[1:]):
        try:
            if str(test) != "nan":
                spli = test.split(":")
                for val in range(len(spli)):
                    flags[val].append(spli[val])
            elif "" == test:
                pass
            else:
                flags[0].append("nan")
        except:
            print("TEST",str(test))
            print("FORMATA",list(formata))
            print("FLAGS",flags)
            time.sleep(50)

    for i in flags:
        text += formatc[i[0]][0]+": "
        for j in i[1:]:
            text+=j+"; "
        text +="\n"

    return text

def write_basics_part(basics, alters):
    text = "GENERAL INFO:\n\n"
    if "chr" in basics[0]:
        Chrom = basics[0]
    else:
        Chrom = "chr"+basics[0]
    Pos = basics[1]
    alter = basics[4]
    text+="Chromosome: "+str(Chrom)+"\n"
    text+="Position: "+ str(Pos)+ "\n"
    text+="ID:" + str(basics[2])+ "\n"
    text+="Base/bases of Reference:" + str(basics[3])+ "\n"
    text+="Alteration: "
    if basics[4] in alters.keys():
        text += alters[basics[4]]+ "\n"
    else:
        text += basics[4]+ "\n"
    text+="Quality: "+ str(basics[5])+ "\n"
    text+="Filter:" + basics[6]+ "\n"

    return [Chrom,Pos, alter], text

def redefine_information(df,information):
    for row in df.iterrows():
        info = row[1][7]
        for item in info.split(";"):
            if "=" in item:
                items = item.split("=")
                try:
                    if information[items[0]]:
                        pass
                except:
                    if items[0] not in information.keys():
                        if is_float(items[1]):
                            information[items[0]] = [items[0]+" (undefined in vcf header): ", "Float", 1]
                        else:
                            information[items[0]] = [items[0]+" (undefined in vcf header): ", "String","A"]
    return information


def get_numeric_columns(df, information):
    numerals = {}
    #!if number = A --> means that instead of number = AGATATCGGAGACGCAGCATGCAT
    for key in information.keys():
        numerals[information[key][0]] = []

    for row in df.iterrows():
        info = row[1][7]
        to_find = list(numerals.keys())
        for item in info.split(";"):
            if "=" in item:
                items = item.split("=")
            if information[items[0]][0] in to_find:
                numerals[information[items[0]][0]].append(items[1])
                to_find.remove(information[items[0]][0])
        for missing in to_find:
            numerals[missing].append("Na")
    return numerals


def main():
    files = sub.run("ls | grep .vcf ", capture_output=True, shell=True)
    for file in files.stdout.decode().split("\n"):
        if file =="":
            break
        print("FILE:", file)
        df , lines = read_data(file)
        alteration,information,formatc = get_metadata_on_vcf(lines)
        c = 0
        information = redefine_information(df, information)
        final_dict= get_numeric_columns(df, information)
        # for i in final_dict.keys():
        #     print(len(final_dict[i]))
        text = []
        chroms = []
        positions = []
        alters = []
        c = 1
        for row in df.iterrows():
            # print("INFO", row[1][7])
            # print("FORMAT",row[1][8:] )
            # print("BASICS",row[1][0:7])
            search,basics_for_line = write_basics_part(row[1][0:7],alteration)
            info_for_line = write_info_part(row[1][7],information)
            format_for_line = write_format_part(row[1][8:],formatc)
            text.append(basics_for_line + info_for_line + format_for_line)
            chroms.append(search[0])
            positions.append(search[1])
            alters.append(search[2])
            print ("formated "+ str(c)+" lines out of "+ str(df.shape[0]), end="\r")
            c+=1
        print ("formated "+ str(c-1)+" lines out of "+ str(df.shape[0]))
        final_dict["Chrom"] = chroms
        final_dict["Base"] = positions
        final_dict["Alter"] = alters
        final_dict["Information"] = text
        # print("Chrom", len(chroms))
        # print("Base", len(positions))
        # print("Alter", len(alters))
        # print("Info", len(text))
        pd.DataFrame(final_dict).to_csv('FORMATED_'+file.strip(".vcf")+'.csv', index=False)
        #print(numerals)
        # print("ALTERATION:\n")
        # print(alteration)
        # print("INFORMATION:\n")
        # print(information)
        # print("FORMAT:\n")
        # print(formatc)







if __name__ == "__main__":
    main()