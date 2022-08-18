import numpy as np
import argparse
import os

def openListFile(fileName, sex):
    dictFiles = {}

    listFile = open(fileName)

    for fileLine in listFile:
        split = fileLine.strip().split()
        if not sex:
            if len(split) == 3:
                dictFiles[split[0]] = {}
                dictFiles[split[0]]["hybrid"] = split[1]
                dictFiles[split[0]]["pgen"] = split[2]
            else:
                print(f"Ignoring line {fileLine} because it has only ({len(split)} columns)")
        else:
            if len(split) == 4:
                dictFiles[split[0]] = {}
                dictFiles[split[0]]["hybrid"] = split[1]
                dictFiles[split[0]]["pgen"] = split[2]
                dictFiles[split[0]]["sex"] = split[3]
            else:
                print(f"Ignoring line {fileLine} because it has only ({len(split)} columns)")

    return dictFiles

def execute(command):
    print(command)
    os.system(command)

def calculateFrq(fileDict, plink2, folder):
    freqDict = {}
    for pop in fileDict:
        pgenFile = fileDict[pop]["pgen"]

        command = f"{plink2} --pfile {pgenFile} --freq --out {folder}/{pop}_freq"
        execute(command)

        file = open(f"{folder}/{pop}_freq.afreq")
        for line in file:
            split = line.strip().split()
            if split[1] not in freqDict:
                freqDict[split[1]] = {}
            if pop not in freqDict[split[1]]:
                freqDict[split[1]][pop] = {}
                freqDict[split[1]][pop]["ALT"] = split[3]
                freqDict[split[1]][pop]["ALT_FREQ"] = split[4]
    return freqDict

def getStatus(fileDict):
    statusDict = {}
    for pop in fileDict:
        pgenFile = fileDict[pop]["pgen"]
        pvarFile = open(f"{pgenFile}.pvar")

        header = True
        for line in pvarFile:
            if header:
                if "#CHROM" in line:
                    header = False
            else:
                split = line.strip().split()
                if split[2] not in statusDict:
                    statusDict[split[2]] = {}
                if pop not in statusDict[split[2]]:
                    if "TYPED" in line:
                        statusDict[split[2]][pop] = 0
                    else:
                        statusDict[split[2]][pop] = 1
    return statusDict

def getPositionCol(headerLine):
    dictFields = {}
    split = headerLine.strip().split()

    interest = ["#CHROM", "POS", "ID", "REF", "ALT", "A1", "OR", "LOG(OR)_SE", "P", "OBS_CT", "U95", "L95"]

    for i in range(0, len(split)):
        if split[i] in interest:
            print(f"Col {split[i]} -> index {i}")
            dictFields[split[i]] = i

    return dictFields

def prepareInputGWAMAOR(fileDict, freqDict, statusDict, folder, name, sex):
    fileToGWAMA = open(f"{folder}/input{name}.in", 'w')

    for pop in fileDict:
        fileOut = open(f"{folder}/input{pop}.in", 'w')
        if sex:
            fileToGWAMA.write(f"{folder}/input{pop}.in\t{fileDict[pop]['sex']}\n")
        else:
            fileToGWAMA.write(f"{folder}/input{pop}.in\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tBETA\tSE\n")
        fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tOR\tOR_95L\tOR_95U\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tOR\tOR_95L\tOR_95U\n")

        hybrid = fileDict[pop]["hybrid"]

        print(f"Open {hybrid} file")

        fileHybrid = open(hybrid)
        header = True

        for line in fileHybrid:
            if header:
                dictColHeader = getPositionCol(line)
                header = False
            else:
                split = line.strip().split()
                if split[dictColHeader["OR"]] != "NA":
                    ID = split[dictColHeader["ID"]]
                    CHR = split[dictColHeader["#CHROM"]]
                    POS = split[dictColHeader["POS"]]
                    IMPUTED = statusDict[ID][pop]
                    N = split[dictColHeader["OBS_CT"]]
                    EA = split[dictColHeader["A1"]]
                    if EA == split[dictColHeader["ALT"]]:
                        NEA = split[dictColHeader["REF"]]
                    else:
                        NEA = split[dictColHeader["ALT"]]

                    if EA == freqDict[ID][pop]["ALT"]:
                        EAF = freqDict[ID][pop]["ALT_FREQ"]
                    else:
                        EAF = 1-float(freqDict[ID][pop]["ALT_FREQ"])

                    OR = split[dictColHeader["OR"]]
                    OR_95L = split[dictColHeader["L95"]]
                    OR_95U = split[dictColHeader["U95"]]

                    fileOut.write(f"{ID}\t{CHR}\t{POS}\t{IMPUTED}\t{N}\t{EA}\t{NEA}\t{EAF}\t{OR}\t{OR_95L}\t{OR_95U}\n")
        fileOut.close()
    return f"{folder}/input{name}.in"


def prepareInputGWAMABeta(fileDict, freqDict, statusDict, folder, name, sex):
    fileToGWAMA = open(f"{folder}/input{name}.in", 'w')

    for pop in fileDict:
        fileOut = open(f"{folder}/input{pop}.in", 'w')
        if sex:
            fileToGWAMA.write(f"{folder}/input{pop}.in\t{fileDict[pop]['sex']}\n")
        else:
            fileToGWAMA.write(f"{folder}/input{pop}.in\n")
        fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tBETA\tSE\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tOR\tOR_95L\tOR_95U\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tOR\tOR_95L\tOR_95U\n")

        hybrid = fileDict[pop]["hybrid"]

        print(f"Open {hybrid} file")

        fileHybrid = open(hybrid)
        header = True

        for line in fileHybrid:
            if header:
                dictColHeader = getPositionCol(line)
                header = False
            else:
                split = line.strip().split()
                if split[dictColHeader["OR"]] != "NA":
                    ID = split[dictColHeader["ID"]]
                    CHR = split[dictColHeader["#CHROM"]]
                    POS = split[dictColHeader["POS"]]
                    IMPUTED = statusDict[ID][pop]
                    N = split[dictColHeader["OBS_CT"]]
                    EA = split[dictColHeader["A1"]]
                    if EA == split[dictColHeader["ALT"]]:
                        NEA = split[dictColHeader["REF"]]
                    else:
                        NEA = split[dictColHeader["ALT"]]

                    if EA == freqDict[ID][pop]["ALT"]:
                        EAF = freqDict[ID][pop]["ALT_FREQ"]
                    else:
                        EAF = 1-float(freqDict[ID][pop]["ALT_FREQ"])

                    BETA = np.log(float(split[dictColHeader["OR"]]))
                    SE = split[dictColHeader["LOG(OR)_SE"]]
                    fileOut.write(f"{ID}\t{CHR}\t{POS}\t{IMPUTED}\t{N}\t{EA}\t{NEA}\t{EAF}\t{BETA}\t{SE}\n")
        fileOut.close()
    return f"{folder}/input{name}.in"

def runGWAMABeta(fileGWAMA, gwama, folder, name, randomEffectCorrection, genomicControl, sex):
    line = f"{gwama} -i {fileGWAMA} -o {folder}/{name} -qt"
    if randomEffectCorrection:
        line = f"{line} -r"
    if genomicControl:
        line = f"{line} -gc"
    if sex:
        line = f"{line} --sex"

    execute(line)

def runGWAMAOR(fileGWAMA, gwama, folder, name, randomEffectCorrection, genomicControl, sex):
    line = f"{gwama} -i {fileGWAMA} -o {folder}/{name}"
    if randomEffectCorrection:
        line = f"{line} -r"
    if genomicControl:
        line = f"{line} -gc"
    if sex:
        line = f"{line} --sex"


    execute(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GWAMA automatic')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-l', '--list', help='List of hybrid files. Format: <NAME> <HYBRID> <PGEN prefix file>',
                          required=True)
    required.add_argument('-f', '--folder', help='Output folder name', required=True)
    required.add_argument('-n', '--name', help='Name to output file', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-o', '--odds', help='Use OR instead BETA (Default: False)', default=False, action="store_true")
    optional.add_argument('-s', '--sex', help='Run gender-differentiated and gender-heterogeneity analysis (Default: False)', default=False,
                          action="store_true")
    optional.add_argument('-r', '--random', help='Random effect correction from GWAMA (Default: False)', default=False,
                          action="store_true")
    optional.add_argument('-g', '--genomic', help='Use genomic control for adjusting studies result files from GWAMA (Default: False)', default=False,
                          action="store_true")

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-G', '--gwama', help='Path to gwama', required=True)
    programs.add_argument('-P', '--plink2', help='Path to PLINK2', required=True)
    args = parser.parse_args()

    execute(f"mkdir {args.folder}")

    fileDict = openListFile(args.list, args.sex)
    freqDict = calculateFrq(fileDict, args.plink2, args.folder)
    statusDict = getStatus(fileDict)
    if args.odds:
        fileGWAMA = prepareInputGWAMAOR(fileDict, freqDict, statusDict, args.folder, args.name, args.sex)
        runGWAMAOR(fileGWAMA, args.gwama, args.folder, args.name, args.random, args.genomic, args.sex)
    else:
        fileGWAMA = prepareInputGWAMABeta(fileDict, freqDict, statusDict, args.folder, args.name, args.sex)
        runGWAMABeta(fileGWAMA, args.gwama, args.folder, args.name, args.random, args.genomic, args.sex)
