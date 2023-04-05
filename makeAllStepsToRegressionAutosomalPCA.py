# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:07:10 2022

@author: PEIXOTT
"""

import os
import argparse
import numpy as np
import gzip
import time

def addPCAOnCovarDict(PCAFile, covarDict, sex):
    fileMale = open(PCAFile)

    dataInfo = {}
    indList = []
    
    header = True
    for line in fileMale:
        if header:
            header = False
            headerLine = line.strip().split()
            
        else:
            split = line.strip().split()
            indList.append(split[0])
            for i in range(1,len(split)):
                if headerLine[i] not in covarDict[split[0]]:
                    covarDict[split[0]][headerLine[i]] = {}
                covarDict[split[0]][headerLine[i]][sex] = float(split[i])
                
                if headerLine[i] not in dataInfo:
                    dataInfo[headerLine[i]] = []
                dataInfo[headerLine[i]].append(float(split[i]))
            covarDict[split[0]]["Outlier"] = False
    
    return covarDict, dataInfo, indList
                
    

def createOutlierFileBySex(PCAMale, PCAFemale, folder, name, covarDict, N):
    covarDict, dataInfoMale, indListMale = addPCAOnCovarDict(PCAMale, covarDict, "Male")
    covarDict, dataInfoFemale, indListFemale = addPCAOnCovarDict(PCAFemale, covarDict, "Female")
    indListOutlierMale = calculateOutilier(dataInfoMale, indListMale, folder, name, covarDict, "Male", N)
    indListOutlierFemale = calculateOutilier(dataInfoFemale, indListFemale, folder, name, covarDict, "Female", N)

    return indListOutlierMale, indListOutlierFemale, covarDict

def createOutlierFile(PCA, folder, name, covarDict, N, sex):
    covarDict, dataInfo, indList = addPCAOnCovarDict(PCA, covarDict, sex)
    indListOutlier = calculateOutilier(dataInfo, indList, folder, name, covarDict, "Both", N)
    
    return indListOutlier, covarDict
    

def calculateOutilier(dataInfo, indList, folder, name, covarDict, sex, N):
    dataInfo["SD"] = {}
    for i in range(1,11):
        PC = f'PC{i}'
        dataInfo['SD'][PC] = np.std(dataInfo[PC])
        
    for i in range(1,10,2):
        lower = f'PC{i}'
        higher = f'PC{i+1}'
        fileOut = open(f"{folder}/{name}_{sex}_{lower}_{higher}_toPlot.txt", 'w')
        fileOut.write(f'IND\t{lower}\t{higher}\tOUTLIER\tSTATUS\tCOUNTRY\tSEX\n')
        
        for ind in indList:
            PCLower = covarDict[ind][lower][sex]
            PCHigher = covarDict[ind][higher][sex]
            
            sdLower = dataInfo['SD'][lower]
            sdHigher = dataInfo['SD'][higher]
            
            statusText = "Case"
            if covarDict[ind]["DISEASE"] == "1":
                statusText = "Control"        
                    
            sexText = "Male"
            if covarDict[ind]["SEX"] == "2":
                sexText = "Female"
            
            if PCLower < (-1*N*sdLower) or PCLower > (N*sdLower) or PCHigher < (-1*N*sdHigher) or PCHigher > (N*sdHigher):
                covarDict[ind]["Outlier"] = True
                fileOut.write(f'{ind}\t{PCLower}\t{PCHigher}\tTRUE\t{statusText}\t{covarDict[ind]["COUNTRY"]}\t{sexText}\n')
            else:
                fileOut.write(f'{ind}\t{PCLower}\t{PCHigher}\tFALSE\t{statusText}\t{covarDict[ind]["COUNTRY"]}\t{sexText}\n')
        fileOut.close()
        
    fileOut = open(f"{folder}/{name}_{sex}_outliersToRemove", 'w')
    for ind in indList:
        if covarDict[ind]["Outlier"]:
           fileOut.write(ind+"\n")
    fileOut.close()
    
    return f"{folder}/{name}_{sex}_outliersToRemove"


def execute(line, run=True):
    print(" ==================== ======================== ===================")
    print(line)
    if run:
        os.system(f"{line}")
    print("******************************************************************")
    #input()
    

def filterData(desiredCountry, folder, name, covarDict):    
            
    if desiredCountry[0].upper() == "ALL":
        print(f'We will use all countries')
    else:
        print(f"We will select: {' '.join(desiredCountry)}")
    
            
    fileMale = open(f'{folder}/toExtract{name}_Male', 'w')
    fileFemale = open(f'{folder}/toExtract{name}_Female', 'w')
    for ID in covarDict:
        keep = True
        
        if desiredCountry[0].upper() != "ALL":
            sameCountry = False
            for country in desiredCountry:
                if country.upper() == covarDict[ID]['COUNTRY'].upper():
                    sameCountry = True
            if not sameCountry:
                keep = False
                
        if keep:
            if covarDict[ID]['SEX'].upper() == "FEMALE" or covarDict[ID]['SEX'].upper() == "2":
                fileFemale.write(f'{ID}\n')
            else:
                if covarDict[ID]['SEX'].upper() == "MALE" or covarDict[ID]['SEX'].upper() == "1":
                    fileMale.write(f'{ID}\n')
                    
    fileFemale.close()
    fileMale.close()
    return f'{folder}/toExtract{name}_Male', f'{folder}/toExtract{name}_Female'
    

def readInformationAboutSamples(covarTable, countryFile):
    print('We are reading the covar table. We are asssuming that the ID is the first col')
    print('We are also assuming that there is the column SEX and the Phenotype column is named DISEASE')
    print('We are also checking if there is any covar field that is empty or NA')
    
    file = open(covarTable)
    
    covarDict = {}
    
    header = True
    for line in file:
        if header:
            header = False
            splitHeader = line.strip().split()
        else:
            split = line.strip().split()
            
            toInclude = True
            for i in range(len(split)):
                if split[i] == "NA" or split[i] == "" or split[i] == " " or split[i] == "nan":
                    toInclude = False
                    print(f'Removing the ind {split[0]} because there is missing data ({split[i]}) on the field {splitHeader[i]}')
            
            if toInclude:
                covarDict[split[0]] = {}
                for i in range(1,len(split)):
                    if splitHeader[i].upper() == "DISEASE":
                        if split[i] == "0" or split[i] == 0:
                            split[i] = "1"
                        else:
                            split[i] = "2"
                    covarDict[split[0]][splitHeader[i].upper()] = split[i]
                    
    
    print('We are reading the ID country file. We are asssuming that the ID is the first col and country is the second')
    
    file = open(countryFile)
    
    header = True
    for line in file:
        if header:
            header = False
        else:
            split = line.strip().split()
            
            if split[0] in covarDict:
                covarDict[split[0]]["COUNTRY"] = split[1]
    
    return covarDict

def getIndListFromVCF(vcfFile, fileOutName):
    fileOut = open(fileOutName, "w")

    file = gzip.open(vcfFile)
    for line in file:
        line = line.decode("utf-8")
        if "#CHROM" in line:
            split = line.strip().split()

            for i in range(9, len(split)):
                fileOut.write(split[i] + "\n")
            break
    fileOut.close()

    return fileOutName

def bcftoolsExtractAutosomal(fileNameX, fileNameAutosomal, fileToExtractMale, fileToExtractFemale, folder, name, run=True):
    bcftoolsIndex(fileNameX, run)
    bcftoolsIndex(fileNameAutosomal, run)
    execute(f'bcftools view -S {fileToExtractMale} -Oz -o {folder}/{name}_Male_X_step1.vcf.gz {fileNameX} --force-samples', run)
    execute(f'bcftools view -S {fileToExtractFemale} -Oz -o {folder}/{name}_Female_X_step1.vcf.gz {fileNameX} --force-samples', run)

    malesFromX = getIndListFromVCF(f"{folder}/{name}_Male_X_step1.vcf.gz", f"{folder}/extractMaleFromX.txt")
    femalesFromX = getIndListFromVCF(f"{folder}/{name}_Female_X_step1.vcf.gz", f"{folder}/extractFemaleFromX.txt")

    execute(f'bcftools view -S {malesFromX} -Oz -o {folder}/{name}_Male_A_step1.vcf.gz {fileNameAutosomal} --force-samples', run)
    execute(f'bcftools view -S {femalesFromX} -Oz -o {folder}/{name}_Female_A_step1.vcf.gz {fileNameAutosomal} --force-samples', run)


    
    return f'{folder}/{name}_Male_A_step1.vcf.gz', f'{folder}/{name}_Female_A_step1.vcf.gz'


def bcftoolsIndex(fileName, run):
    execute(f'bcftools index {fileName}', run)

def removeOutiler(listToRemove, vcfFile, name, folder, sex, run = True):
    bcftoolsIndex(vcfFile, run)
    execute(f'bcftools view -S ^{listToRemove} -Oz -o {folder}/{name}_{sex}_withoutOutlier.vcf.gz {vcfFile} --force-samples', run)
    
    return f'{folder}/{name}_{sex}_withoutOutlier.vcf.gz'
    

def removeOutliersMalesAndFemales(outlierFileMale, outlierFileFemale, fileNameMale, fileNameFemale, name, folder, run = True):
    maleWithoutOutlier = removeOutiler(outlierFileMale, fileNameMale, name, folder, "Male", run)
    femaleWithoutOutlier = removeOutiler(outlierFileFemale, fileNameFemale, name, folder, "Female", run)
    
    return maleWithoutOutlier, femaleWithoutOutlier

def convertAndRemoveLDBySex(fileNameMale, fileNameFemale, folder, name, plink2, run=True):
    withoutLDMale = convertAndRemoveLD(fileNameMale, folder, name, 'Male', plink2, run)
    withoutLDFemale = convertAndRemoveLD(fileNameFemale, folder, name, 'Female', plink2, run)
    return withoutLDMale, withoutLDFemale

def getIndFromMalesAndFemales(male, female, vcfFile, folder, name, run= True):
    
    fileOut = open(f'{folder}/ToKeep_FemalesAndMales.txt', 'w')
    fileOutLog = open(f'{folder}/ToKeep_FemalesAndMales.log', 'w')
    
    file = gzip.open(male)
    for line in file:
        line = line.decode("utf-8")
        if "#CHROM" in line:
            split = line.strip().split()
            
            for i in range(9, len(split)):
                fileOut.write(split[i]+ "\n")
                fileOutLog.write(split[i]+ "\tM\n")
            break

    file = gzip.open(female)
    for line in file:
        line = line.decode("utf-8")
        if "#CHROM" in line:
            split = line.strip().split()
            
            for i in range(9, len(split)):
                fileOut.write(split[i]+ "\n")
                fileOutLog.write(split[i]+ "\tF\n")
            break
    fileOut.close()
    execute(f'bcftools view -S {folder}/ToKeep_FemalesAndMales.txt -Oz -o {folder}/{name}_Both_step1.vcf.gz {vcfFile} --force-samples', run)
    return f'{folder}/{name}_Both_step1.vcf.gz'
    
    
    
    

def convertAndRemoveLD(fileName, folder, name, sex, plink2, run = True):
    execute(f"{plink2} --vcf {fileName} --make-pgen --out {folder}/{name}_{sex}_toLD", run)
    
    execute(f"{plink2} --pfile {folder}/{name}_{sex}_toLD --out {folder}/{name}_{sex} --indep-pairwise 200 50 0.2", run)
    
    file = open(f'{folder}/{name}_{sex}.prune.in')
    SNPs = []
    for line in file:
        line = line.strip()
        SNPs.append(line)
    file.close()
    
    file = open(f'{folder}/{name}_{sex}_toLD.pvar')
    fileToKeep = open(f'{folder}/{name}_{sex}_withoutLD.txt', 'w')

    print("Creating the list of variants to keep")
    if run:
        header = True
        for line in file:
            if header:
                if "#CHROM" in line:
                    header = False
            else:
                split = line.split()
                if split[2] in SNPs:
                    fileToKeep.write(f'chr{split[0]}\t{split[1]}\n')
        fileToKeep.close()

    execute(f"bcftools view -T {folder}/{name}_{sex}_withoutLD.txt -Oz -o {folder}/{name}_{sex}_withoutLD.vcf.gz {fileName}", run)
    return f'{folder}/{name}_{sex}_withoutLD.vcf.gz'
    
def runPCARBySex(vcfFileMale, vcfFileFemale, folder, name, PCA, run = True):
    maleTSV = runPCA(vcfFileMale, folder, name, "Male", PCA, run)
    femaleTSV = runPCA(vcfFileFemale, folder, name, "Female", PCA, run)
    
    return maleTSV, femaleTSV
    
def runPCA(vcfFile, folder, name, sex, PCA, run=True):
    execute(f'Rscript {PCA} {vcfFile} {folder}/{name}_{sex}.gds {folder}/{name}_{sex}.tsv', run)
    
    return f'{folder}/{name}_{sex}.tsv'

def  convertToPLINK2AndRun(vcfFile, dictCovarLocal, vcfImputed, folder, name, sex, plink2, covar, cutoff, firth):
    
    
    #List of individuals to be extract from imputed data
    fileOut = open(f'{folder}/{name}_{sex}_ToExtractFromImputed', 'w')
    
    file = gzip.open(vcfFile)
    for line in file:
        line = line.decode("utf-8")
        if "#CHROM" in line:
            split = line.strip().split()
            
            for i in range(9, len(split)):
                fileOut.write(split[i]+ "\n")
            break
    fileOut.close()

    command = f'{plink2} --vcf {vcfImputed} --keep {folder}/{name}_{sex}_ToExtractFromImputed --make-pgen --out ' \
                 f'{folder}/{name}_{sex}_toRegression --extract-if-info \"R2 > {cutoff}\"'

    execute(command)
    execute(f'mkdir {folder}/backPSAM')
    execute(f'mkdir {folder}/result')
    execute(f'mv {folder}/{name}_{sex}_toRegression.psam {folder}/backPSAM/')
    
    filePSAMOriginal = open(f'{folder}/backPSAM/{name}_{sex}_toRegression.psam')
    filePSAMWithCovar = open(f'{folder}/{name}_{sex}_toRegression.psam', 'w')
    
    filePSAMWithCovar.write(f'#IID')

    header = True
    fieldsSet = False
    for line in filePSAMOriginal:
        if header:
            header = False
        else:
            split = line.split()
            ind = split[0]
            if not fieldsSet:
                headerOrder = []
                for field in dictCovarLocal[ind]:
                    headerOrder.append(field)
                    filePSAMWithCovar.write(f'\t{field}')
                filePSAMWithCovar.write('\n')
                fieldsSet = True
                
            
            filePSAMWithCovar.write(f'{ind}')
            print(f"{ind} -> {dictCovarLocal[ind]}")
            for field in headerOrder:
                print(f"{field} -> ", end="")
                if "PC" in field:
                    filePSAMWithCovar.write(f'\t{dictCovarLocal[ind][field][sex]}')
                else:
                     filePSAMWithCovar.write(f'\t{dictCovarLocal[ind][field]}')
            print(f"\n")
            filePSAMWithCovar.write(f'\n')
            
    filePSAMWithCovar.close()

    if not firth:
        command = f'{plink2} --pfile {folder}/{name}_{sex}_toRegression --pheno-name DISEASE --covar-variance-standardize ' \
              f'--glm hide-covar --out {folder}/result/{name}_{sex} --covar-name {covar} --ci 0.95'
    else:
        command = f'{plink2} --pfile {folder}/{name}_{sex}_toRegression --pheno-name DISEASE --covar-variance-standardize ' \
                  f'--glm hide-covar firth --out {folder}/result/{name}_{sex} --covar-name {covar} --ci 0.95'
    execute(command)

    return f"{folder}/{name}_{sex}_toRegression", f"{folder}/result/{name}_{sex}"
    
def buildCovarList(covarList, maxPC):
    sex = ""
    both = ""

    for covar in covarList:
        if covar.lower() != "sex":
            if sex == "":
                sex = covar
                both = covar
            else:
                sex = sex + " " + covar
                both = both + " " + covar

    for i in range(1, maxPC+1):
        if i == 1:
            sex = f"{sex} PC{i}"
            both = f"{both} PC{i}"
        else:
            sex = f"{sex} PC{i}"
            both = f"{both} PC{i}"

    return sex, both

def  removeHeterozygous(vcf, folder, name):
    import gzip

    inputFile = gzip.open(vcf)
    fileOut = open(f"{folder}/{name}_withoutHeterozygous.vcf", 'w')

    header = True
    for line in inputFile:
        line = line.decode('utf-8')
        if header:
            if line[0:6] == "#CHROM":
                split = line.strip().split()

                fileOut.write(f"{split[0]}")
                for i in range(1, len(split)):
                    fileOut.write(f"\t{split[i]}")
                fileOut.write("\n")
                header = False
            else:
                fileOut.write(line)
        else:
            split = line.strip().split()
            fileOut.write(f"{split[0]}")
            for i in range(1, len(split)):
                field = split[i]
                if field[0:3] == "0|1" or field[0:3] == "1|0":
                    data = field.split(":")
                    outString = ".|."
                    for j in range(1, len(data)):
                        outString = f"{outString}:."
                    fileOut.write(f"\t{outString}")
                else:
                    fileOut.write(f"\t{split[i]}")
            fileOut.write("\n")
    fileOut.close()

    return f"{folder}/{name}_withoutHeterozygous.vcf"


def runMetaMaleFemale(pfilesMale, regressionMale, pfilesFemale, regressionFemale, gwama, python, plink, folder, name, firth):
    newFolder = f"{folder}/{name}FemaleMale/"
    execute(f"mkdir {newFolder}")
    fileInput = open(f"{newFolder}/{name}toMetaAnalyse.txt", 'w')

    if not firth:
        regMale = f'{regressionMale}.DISEASE.glm.logistic.hybrid'
        regFemale = f'{regressionFemale}.DISEASE.glm.logistic.hybrid'
    else:
        regMale = f'{regressionMale}.DISEASE.glm.firth'
        regFemale = f'{regressionFemale}.DISEASE.glm.firth'

    fileInput.write(f"{name}_Female\t{regFemale}\t{pfilesFemale}\tF\n")
    fileInput.write(f"{name}_Male\t{regMale}\t{pfilesMale}\tM\n")
    fileInput.close()

    execute(f"{python} metaAnalysisGWAMA.py -l {newFolder}/{name}toMetaAnalyse.txt -n {name}_MetaFemaleMale "
            f"-f {newFolder} -G {gwama} -P {plink} -o -s")

def convertPLINK2VCF(autosomal, genotyped, folder, name, plink2, bcftools, run = True):
    execute(f"{bcftools} query {genotyped} -l > {folder}/toKeepChrX.txt", run)
    #time.sleep(2)


    fileInput = open(f"{folder}/toKeepChrX.txt")
    fileToUse = open(f"{folder}/toKeepChrX_2cols.txt", 'w')
    for line in fileInput:
        ID = line.strip()
        fileToUse.write(f"{ID}\t{ID}\n")
    fileToUse.close()
    fileInput.close()

    execute(f"{plink2} --bfile {autosomal} --keep {folder}/toKeepChrX_2cols.txt --recode vcf id-paste=iid --out {folder}/{name} --output-chr chr26", run)
    #input()
    return f"{folder}/{name}.vcf"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    data = parser.add_argument_group("Data arguments")
    data.add_argument('-g', '--genotyped', help='Genotyped file name', required=False)
    data.add_argument('-a', '--autosomal', help='Autosomal Genotyped file name', required=False)

    data.add_argument('-i', '--imputed', help='Imputed file name', required=False)
    data.add_argument('-t', '--tableCovar', help='File with covariatives to be added to the model', required=True)
    data.add_argument('-C', '--countryFile', help='File with relation Ind country', required=True)
    data.add_argument('-c', '--country',
                          help='Country to analyze (default: all). You can select more than one country',
                          required=False, default=["all"], nargs='+')

    parameters = parser.add_argument_group("Parameter arguments")
    parameters.add_argument('-l', '--covarList', help='List of covar to be used (do not include PCAs)', required=True, nargs="+")
    parameters.add_argument('-m', '--maxPC', help='max PC to be used on covar', required=True, type=int)
    parameters.add_argument('-r', '--r2', help='r2 cutoff', required=True, type=float)
    parameters.add_argument('-F', '--firth', help='Force all PLINK2 regressions use the firth', required=False, default = False, action="store_true")
    parameters.add_argument('-H', '--homozygousOnly', help='Remove heterozygous from imputed file', required=False,
                            default=False, action="store_true")

    output = parser.add_argument_group("Output arguments")
    output.add_argument('-n', '--name', help='Name to use', required=False)
    output.add_argument('-f', '--folder', help='Folder to output files', required=False)


    programs = parser.add_argument_group("Programs arguments")
    programs.add_argument('-G', '--gwama', help='GWAMA program (default = gwama)', required=False, default="gwama")
    programs.add_argument('-p', '--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")
    programs.add_argument('-P', '--python', help='Path of Python 3 (default = python)', required=False, default="python")
    programs.add_argument('-R', '--runPCA', help='Path of runPCA.R script (default = runPCA.R)', required=False,
                          default="runPCA.R")
    args = parser.parse_args()
    
    execute(f'mkdir {args.folder}')
    
    covarDict = readInformationAboutSamples(args.tableCovar, args.countryFile)
    fileToExtractMale, fileToExtractFemale = filterData(args.country, args.folder, args.name, covarDict)
    
    autosomalVCF = convertPLINK2VCF(args.autosomal, args.genotyped, args.folder, args.name, args.plink2, "bcftools")

    #==================================== Extract ================================================
    fileExtractedMaleA, fileExtractedFemaleA = bcftoolsExtractAutosomal(args.genotyped, autosomalVCF, fileToExtractMale,
                                                                        fileToExtractFemale, args.folder, args.name)


    # ==================================== Remove LD ================================================
    fileWithoutLDMale, fileWithoutLDFemale = convertAndRemoveLDBySex(fileExtractedMaleA, fileExtractedFemaleA, args.folder,
                                                                     args.name+"_Autosomal", args.plink2)

    # ==================================== Run PCA ================================================
    PCAFileMale, PCAFileFemale = runPCARBySex(fileWithoutLDMale, fileWithoutLDFemale, args.folder, args.name+"_Autosomal",
                                              args.runPCA)

    # ==================================== Outliers ================================================
    outlierFileMale, outlierFileFemale, dictCovar = createOutlierFileBySex(PCAFileMale, PCAFileFemale, args.folder,
                                                                           args.name, covarDict, 3)

    # ==================================== Remove outliers ================================================
    maleWithoutOutlierA, femaleWithoutOutlierA = removeOutliersMalesAndFemales(outlierFileMale, outlierFileFemale,
                                                                             fileExtractedMaleA, fileExtractedFemaleA,
                                                                             args.name+"_Autosomal", args.folder)

    # ==================================== Both ================================================
    bothBeginAutosomal = getIndFromMalesAndFemales(maleWithoutOutlierA, femaleWithoutOutlierA, autosomalVCF, args.folder, args.name+"_Autosomal")
    fileWithoutLDBoth = convertAndRemoveLD(bothBeginAutosomal, args.folder, args.name+"_Autosomal", 'Both', args.plink2)
    PCABoth = runPCA(fileWithoutLDBoth, args.folder, args.name, "Both", args.runPCA)

    outlierBoth, dictCovar = createOutlierFile(PCABoth, args.folder, args.name, covarDict, 3, "Both")
    bothWithoutOutlierA = removeOutiler(outlierBoth, bothBeginAutosomal, args.name, args.folder, "both_OutlierMaleAndFemale")

    print("Building covar list (covariates + PCs). In our tests PLINK2 add automatically the SEX, causing this error message:")
    print("Error: Cannot proceed with --glm regression on phenotype 'DISEASE', since correlation between covariates "
          "'SEX' and 'SEX' is too high (CORR_TOO_HIGH). You may want to remove redundant covariates and try again.")

    covarSex, covarBoth = buildCovarList(args.covarList, args.maxPC)

    if args.homozygousOnly:
        args.imputed = removeHeterozygous(args.imputed, args.folder, args.name)

    print(dictCovar)
    pfilesMale, regressionMale = convertToPLINK2AndRun(maleWithoutOutlierA, dictCovar, args.imputed, args.folder, args.name, "Male", args.plink2, covarSex, args.r2, args.firth)
    pfilesFemale, regressionFemale = convertToPLINK2AndRun(femaleWithoutOutlierA, dictCovar, args.imputed, args.folder, args.name, "Female", args.plink2, covarSex, args.r2, args.firth)
    pfilesBoth, regressionBoth = convertToPLINK2AndRun(bothWithoutOutlierA, dictCovar, args.imputed, args.folder, args.name, "Both", args.plink2, covarBoth, args.r2, args.firth)

    runMetaMaleFemale(pfilesMale, regressionMale, pfilesFemale, regressionFemale, args.gwama, args.python, args.plink2, args.folder, args.name, args.firth)