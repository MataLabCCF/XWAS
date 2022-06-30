# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:07:10 2022

@author: PEIXOTT
"""

import os
import argparse
import numpy as np
import gzip


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
    
    
        
    

def execute(line):
    print(" ==================== ======================== ===================")
    print(line)
    os.system(line)
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
    print('We are also checking if there is anyt covar field that is empty or NA')
    
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

def bcftoolsExtract(fileName, fileToExtractMale, fileToExtractFemale, folder, name):
    bcftoolsIndex(fileName)
    execute(f'bcftools view -S {fileToExtractMale} -Oz -o {folder}/{name}_Male_step1.vcf.gz {fileName} --force-samples')
    execute(f'bcftools view -S {fileToExtractFemale} -Oz -o {folder}/{name}_Female_step1.vcf.gz {fileName} --force-samples')
    
    return f'{folder}/{name}_Male_step1.vcf.gz', f'{folder}/{name}_Female_step1.vcf.gz'


def bcftoolsIndex(fileName):
    execute(f'bcftools index {fileName}')

def removeOutiler(listToRemove, vcfFile, name, folder, sex):
    bcftoolsIndex(vcfFile)
    execute(f'bcftools view -S ^{listToRemove} -Oz -o {folder}/{name}_{sex}_withoutOutlier.vcf.gz {vcfFile} --force-samples')
    
    return f'{folder}/{name}_{sex}_withoutOutlier.vcf.gz'
    

def removeOutliersMalesAndFemales(outlierFileMale, outlierFileFemale, fileNameMale, fileNameFemale, name, folder):
    maleWithoutOutlier = removeOutiler(outlierFileMale, fileNameMale, name, folder, "Male")
    femaleWithoutOutlier = removeOutiler(outlierFileFemale, fileNameFemale, name, folder, "Female")
    
    return maleWithoutOutlier, femaleWithoutOutlier

def convertAndRemoveLDBySex(fileNameMale, fileNameFemale, folder, name, plink2):
    withoutLDMale = convertAndRemoveLD(fileNameMale, folder, name, 'Male', plink2)
    withoutLDFemale = convertAndRemoveLD(fileNameFemale, folder, name, 'Female', plink2)
    return withoutLDMale, withoutLDFemale

def getIndFromMalesAndFemales(male, female, vcfFile, folder, name):
    
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
    execute(f'bcftools view -S {folder}/ToKeep_FemalesAndMales.txt -Oz -o {folder}/{name}_Both_step1.vcf.gz {vcfFile} --force-samples')
    return f'{folder}/{name}_Both_step1.vcf.gz'
    
    
    
    

def convertAndRemoveLD(fileName, folder, name, sex, plink2):
    execute(f"{plink2} --vcf {fileName} --make-pgen --out {folder}/{name}_{sex}_toLD")
    
    execute(f"{plink2} --pfile {folder}/{name}_{sex}_toLD --out {folder}/{name}_{sex} --indep-pairwise 200 50 0.2")
    
    file = open(f'{folder}/{name}_{sex}.prune.in')
    SNPs = []
    for line in file:
        line = line.strip()
        SNPs.append(line)
    file.close()
    
    file = open(f'{folder}/{name}_{sex}_toLD.pvar')
    fileToKeep = open(f'{folder}/{name}_{sex}_withoutLD.txt', 'w')
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

    execute(f"bcftools view -T {folder}/{name}_{sex}_withoutLD.txt -Oz -o {folder}/{name}_{sex}_withoutLD.vcf.gz {fileName}")
    return f'{folder}/{name}_{sex}_withoutLD.vcf.gz'
    
def runPCARBySex(vcfFileMale, vcfFileFemale, folder, name):
    maleTSV = runPCA(vcfFileMale, folder, name, "Male")
    femaleTSV = runPCA(vcfFileFemale, folder, name, "Female")
    
    return maleTSV, femaleTSV
    
def runPCA(vcfFile, folder, name, sex):
    execute(f'Rscript runPCA.R {vcfFile} {folder}/{name}_{sex}.gds {folder}/{name}_{sex}.tsv')
    
    return f'{folder}/{name}_{sex}.tsv'

def  convertToPLINK2AndRun(vcfFile, dictCovarLocal, vcfImputed, folder, name, sex, plink2):
    
    
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
    
    execute(f'{plink2} --vcf {vcfImputed} --keep {folder}/{name}_{sex}_ToExtractFromImputed --make-pgen --out {folder}/{name}_{sex}_toRegression --extract-if-info \"R2 > 0.8\"')
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
            for field in headerOrder:
                if "PC" in field:
                    filePSAMWithCovar.write(f'\t{dictCovarLocal[ind][field][sex]}')
                else:
                     filePSAMWithCovar.write(f'\t{dictCovarLocal[ind][field]}')
            filePSAMWithCovar.write(f'\n')
            
    filePSAMWithCovar.close()
    
    if sex.upper() != "BOTH":
        execute(f'{plink2} --pfile {folder}/{name}_{sex}_toRegression --pheno-name DISEASE --covar-variance-standardize --glm hide-covar --out {folder}/result/{name}_{sex} --covar-name PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10')
    else:    
        execute(f'{plink2} --pfile {folder}/{name}_{sex}_toRegression --pheno-name DISEASE --covar-variance-standardize --glm hide-covar --out {folder}/result/{name}_{sex} --covar-name PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10')
    
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-g', '--genotyped', help='Genotyped file name', required=False)
    required.add_argument('-i', '--imputed', help='Imputed file name', required=False)
    required.add_argument('-n', '--name', help='Name to use', required=False)
    required.add_argument('-f', '--folder', help='Folder to output files', required=False)
    required.add_argument('-t', '--tableCovar', help='File with covariatives to be added to the model', required=True)
    required.add_argument('-C', '--countryFile', help='File with relation Ind country', required=True)
    required.add_argument('-c', '--country', help='Country to analyze (default: all). You can select more than one country', required=False, default = ["all"], nargs='+')
    required.add_argument('-p', '--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")

    args = parser.parse_args()
    
    execute(f'mkdir {args.folder}')
    
    covarDict = readInformationAboutSamples(args.tableCovar, args.countryFile)
    fileToExtractMale, fileToExtractFemale = filterData(args.country, args.folder, args.name, covarDict)
    
    
    fileExtractedMale, fileExtractedFemale = bcftoolsExtract(args.genotyped, fileToExtractMale, fileToExtractFemale, args.folder, args.name)
    fileWithoutLDMale, fileWithoutLDFemale = convertAndRemoveLDBySex(fileExtractedMale, fileExtractedFemale, args.folder, args.name, args.plink2)
    
    PCAFileMale, PCAFileFemale = runPCARBySex(fileWithoutLDMale, fileWithoutLDFemale, args.folder, args.name)
    outlierFileMale, outlierFileFemale, dictCovar = createOutlierFileBySex(PCAFileMale, PCAFileFemale, args.folder, args.name, covarDict, 3)
    maleWithoutOutlier, femaleWithoutOutlier = removeOutliersMalesAndFemales(outlierFileMale, outlierFileFemale, fileExtractedMale, fileExtractedFemale, args.name, args.folder)
    
    bothBegin = getIndFromMalesAndFemales(maleWithoutOutlier, femaleWithoutOutlier, args.genotyped, args.folder, args.name)
    fileWithoutLDBoth = convertAndRemoveLD(bothBegin, args.folder, args.name, 'Both', args.plink2)
    PCABoth = runPCA(fileWithoutLDBoth, args.folder, args.name, "Both")
    outlierBoth, dictCovar = createOutlierFile(PCABoth, args.folder, args.name, covarDict, 3, "Both")
    bothWithoutOutlier = removeOutiler(outlierBoth, bothBegin, args.name, args.folder, "both_OutlierMaleAndFemale")
    
    convertToPLINK2AndRun(maleWithoutOutlier, dictCovar, args.imputed, args.folder, args.name, "Male", args.plink2)
    convertToPLINK2AndRun(femaleWithoutOutlier, dictCovar, args.imputed, args.folder, args.name, "Female", args.plink2)
    convertToPLINK2AndRun(bothWithoutOutlier, dictCovar, args.imputed, args.folder, args.name, "Both", args.plink2)
