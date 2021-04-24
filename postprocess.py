#! /usr/bin/python

import os

cases = (["Prism", 
         "Unit Cube", 
         "Triangulated Prism", 
         "Triangulated Hexahedron", 
         "Symmetric Prism", 
         "Symmetric Hexahedron", 
         "Stellated Dodecahedron", 
         "Stellated Icosahedron"])

tools = (["irl","r3d","voftools"])

def readFile(filename):
    read_file = open(filename,'r')
    lines = []
    for line in read_file:
        lines.append(line)
    read_file.close()
    return lines

def processFileOutput(file_dict):
    # Scaling to convert from ncases and seconds to milliseconds per 10^6 cases.
    print("Number of cases: "+ file_dict["irl"][0].split(" ")[0].strip())
    scaling = 1.0
    max_planes = int(float(file_dict["irl"][0].split(" ")[1].strip()))
    for tool_name in file_dict:
        if(int(float(file_dict[tool_name][0].split(" ")[1].strip())) != max_planes):
            print("Different max planes length provided in the result files.\n")
            exit()    

    initialization_comp = {}
    intersection_comp = {}
    volume_comp = {}
    total_comp = {}
    for tool_name in file_dict:
        initialization_comp[tool_name] = {}
        intersection_comp[tool_name] = {}
        volume_comp[tool_name] = {}
        total_comp[tool_name] = {}        
        current_line = 2
        for case in cases:
            initialization_comp[tool_name][case] = {}
            intersection_comp[tool_name][case] = {}
            volume_comp[tool_name][case] = {}
            total_comp[tool_name][case] = {}            
            for line in file_dict[tool_name][current_line:current_line+max_planes]:
                results = line.strip().split(" ")
                plane = int(float(results[0].strip()))
                initialization_comp[tool_name][case][plane] = float(results[1].strip())*scaling
                intersection_comp[tool_name][case][plane] = float(results[2].strip())*scaling
                volume_comp[tool_name][case][plane] = float(results[3].strip())*scaling
                total_comp[tool_name][case][plane] = float(results[4].strip())*scaling       
            current_line += max_planes+1
    return  max_planes, initialization_comp, intersection_comp, volume_comp, total_comp


def getTimingHeader():
    csize = 16
    return "NPlane".center(11) + \
        "IRL".center(csize) + "IRL/IRL".center(csize) +\
        "R3D".center(csize) + "R3D/IRL".center(csize) + \
        "VOFTOOLS".center(csize) + "VOFTOOLS/IRL".center(csize)

def getTimes(timing_dict, case, max_planes):
    s = ""
    for p in range(1, max_planes+1):
        s+= str(p).ljust(9)
        for tool in tools:
            s += ("{0:15.4E}".format(timing_dict[tool][case][p])+" ")
            s += ("{0:15.4E}".format(timing_dict[tool][case][p]/timing_dict["irl"][case][p])+" ")            
        s += "\n"
    return s

def prepareOutput(max_planes, initialization_comp, intersection_comp, volume_comp, total_comp):
    output = {}
    for case in cases:
        output[case] = {}
        output[case]["header"] = "\n\n*********************************************************************\n"
        output[case]["header"] += "Time comparison for Random Planes Intersecting a "+case+"\n"
        output[case]["header"] += "********************************************************************\n"        
        output[case]["initialization"] = "Initialization Time Comparison\n"
        output[case]["initialization"] += (getTimingHeader()+"\n")
        output[case]["initialization"] += (getTimes(initialization_comp, case, max_planes))
        
        output[case]["intersection"] = "Intersection Time Comparison\n"
        output[case]["intersection"] += (getTimingHeader()+"\n")
        output[case]["intersection"] += (getTimes(intersection_comp, case, max_planes))

        output[case]["volume"] = "Volume Time Comparison\n"
        output[case]["volume"] += (getTimingHeader()+"\n")
        output[case]["volume"] += (getTimes(volume_comp, case, max_planes))

        output[case]["total"] = "Total Time Comparison\n"
        output[case]["total"] += (getTimingHeader()+"\n")
        output[case]["total"] += (getTimes(total_comp, case, max_planes))         
    return output

def printToScreen(case_output):
    "\n**Times are given in milliseconds per 10^6 random cases.**\n"
    for case in cases:
        print(case_output[case]["header"])
        print(case_output[case]["initialization"])
        print(case_output[case]["intersection"])
        print(case_output[case]["volume"])
        print(case_output[case]["total"])    
        
def printToResultsFolder(case_output):
    cat = ["initialization","intersection","volume","total"]
    for case in cases:
        for stage in cat:
            filename = ("./results/"+case+"_"+stage+".txt").replace(" ","_")
            f = open(filename,'w')
            f.write(case_output[case][stage])
            f.close()

def fT(a_number):
    return "{0:8.3F}".format(a_number)

def fTR(a_number):
    return "{0:10.2F}".format(a_number)
            
def printLatexTables(max_planes, initialization_comp, intersection_comp, volume_comp, total_comp):
    tables = ["total","total_rel","initialization", "intersection","volume","multiple_plane"]
    output = {}
    output["initialization"] = {}
    output["intersection"] = {}
    output["volume"] = {}
    output["total"] = {}
    output["total_rel"] = {}    
    for case in cases:

        output["total"][case] = (case + " & " +
                                 fT(total_comp["irl"][case][1]) + " & " +
                                 fT(total_comp["r3d"][case][1]) + " & " +
                                 fT(total_comp["voftools"][case][1]) +
                                 "\\\\ \hline")                               
        
        output["total_rel"][case] = (case + " & " +
                                     fTR(total_comp["irl"][case][1]/total_comp["irl"][case][1]) + " & " +
                                     fTR(total_comp["r3d"][case][1]/total_comp["irl"][case][1]) + " & " +
                                     fTR(total_comp["voftools"][case][1]/total_comp["irl"][case][1]) +
                                     "\\\\ \hline")        
        
        output["initialization"][case] = (case + " & " +
                                          fT(initialization_comp["irl"][case][1]) + " & " +
                                          fT(initialization_comp["r3d"][case][1]) + " & " +
                                          fT(initialization_comp["voftools"][case][1]) +
                                          "\\\\ \hline")                                        
        
        output["intersection"][case] = (case + " & " +
                                        fT(intersection_comp["irl"][case][1]) + " & " +
                                        fT(intersection_comp["r3d"][case][1]) + " & " +
                                        fT(intersection_comp["voftools"][case][1]) +
                                        "\\\\ \hline")                                        
        
        output["volume"][case] = (case + " & " +
                                  fT(volume_comp["irl"][case][1]) + " & " +
                                  fT(volume_comp["r3d"][case][1]) + " & " +
                                  fT(volume_comp["voftools"][case][1]) +
                                  "\\\\ \hline")                                        
    
    nplane_cases = ["Unit Cube","Triangulated Hexahedron", "Symmetric Hexahedron"]
    output["multiple_plane"] = {}
    for case in nplane_cases:
        output["multiple_plane"][case] = ("\multirow{5}{*}{"+case+"} & 1 " +
                                          " & " + fT(total_comp["irl"][case][1]) +
                                          " & " + fT(total_comp["r3d"][case][1]) +
                                          " &  " +fT(total_comp["voftools"][case][1]) +
                                          "\\\\ \n" )
        for p in range (2,max_planes+1):
            output["multiple_plane"][case] += (" & " + str(p) + " & " +
                                               fT(total_comp["irl"][case][p])+ " & " +
                                               fT(total_comp["r3d"][case][p]) + " & " +
                                               fT(total_comp["voftools"][case][p])) + "\\\\ \n"            
        output["multiple_plane"][case] += "\hline\n"

    
    
    f = open("./results/latex_output.txt",'w')
    for table in tables:
        f.write("\n\nCreating table for "+table+"\n")        
        if(table != "multiple_plane"):
            for case in cases:
                f.write(output[table][case]+"\n")
                
        else:
            for case in nplane_cases:
                f.write(output[table][case]+"\n")                
    f.close()
        

    
if __name__ == "__main__":
    file_dict = ({
        "irl": readFile("irl_timing.txt"),
        "r3d": readFile("r3d_timing.txt"),
        "voftools": readFile("voftools_timing.txt"),        
        })
    max_planes, initialization_comp, intersection_comp, volume_comp, total_comp  = processFileOutput(file_dict)

    case_output = prepareOutput(max_planes, initialization_comp, intersection_comp, volume_comp, total_comp)

    printToScreen(case_output)

    if not os.path.exists("./results"):
        os.mkdir("./results")    

    printToResultsFolder(case_output)

    printLatexTables(max_planes, initialization_comp, intersection_comp, volume_comp, total_comp)    
       
        


