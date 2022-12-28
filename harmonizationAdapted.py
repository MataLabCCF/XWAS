import os
import time
import gzip
import argparse
import numpy as np
from scipy.stats import fisher_exact

# Global variables
plink = ""
plink2 = ""
bcftools = ""
python = ""
eagle = ""
NAToRA = ""
king = ""
bgzip = ""

def execute(command):
    print(f'\n\n{command}\n\n')
    os.system(command)


def setFamToCovar(fileName, covar, phenoName, continuous, outputPrefix, folder):
    toRemoveMissing = open(f"{folder}/{output}_missingCovar", "w")

    file = open(covar)
    header = True

    binaryFix = False

    covarDict = {}
    phenoIndex = -1
    sexIndex = -1
    idIndex = -1
    for line in file:
        if header:
            headerSplit = line.strip().split()
            for i in range(len(headerSplit)):
                if headerSplit[i].lower() == phenoName.lower():
                    phenoIndex = i
                elif headerSplit[i].lower() == "sex":
                    sexIndex = i
                elif headerSplit[i].lower() == "id" or headerSplit[i].lower() == "iid":
                    idIndex = i
            header = False
        else:
            split = line.strip().split()
            id = split[idIndex]
            sex = split[sexIndex]
            pheno = split[phenoIndex]

            covarDict[id] = {}
            covarDict[id]['sex'] = sex
            covarDict[id]['pheno'] = pheno

            if not continuous:
                if pheno == '0' and not binaryFix:
                    print("We have pheno with value = 0. We will consider case as 2 and control as 1 (PLINK requires)")
                    binaryFix = True


    command = f"cp {fileName}.bed {folder}/{output}.bed"
    execute(command)
    command = f"cp {fileName}.bim {folder}/{output}.bim"
    execute(command)
    file = open(f"{fileName}.fam")
    fileOut = open(f"{folder}/{output}.fam", "w")
    for line in file:
        split = line.strip().split()
        if split[1] in covarDict:
            sex = covarDict[split[1]]['sex']
            pheno = covarDict[split[1]]['pheno']

            if invalid(sex) or invalid(pheno):
                toRemoveMissing.write(f'{split[1]} {split[1]}\n')
                fileOut.write(line)
            else:
                fileOut.write(f'{split[1]} {split[1]} {split[2]} {split[3]} ')
                if sex[0].lower() == 'm':
                    fileOut.write(f'1 ')
                elif sex[0].lower() == 'f':
                    fileOut.write(f'2 ')
                else:
                    fileOut.write(f'{sex} ')
                if binaryFix:
                    pheno = int(pheno) + 1
                fileOut.write(f'{pheno}\n')
        else:
            toRemoveMissing.write(f'{split[1]} {split[1]}\n')
            fileOut.write(f'{split[1]} {split[1]} {split[2]} {split[3]} {split[4]} {split[5]}\n')

    fileOut.close()
    toRemoveMissing.close()

    command = lineToRemoveIndividual(f"{folder}/{output}",f"{folder}/{output}_missingCovar", f"{folder}/{output}_nonNA")
    execute(command)

    return f"{folder}/{output}_nonNA", "bed"

def lineToRemoveIndividual(fileIn, fileToRemove, fileOut):
    return f"plink --bfile {fileIn} --remove {fileToRemove} --make-bed --out {fileOut}"

def lineToRemoveSNPs(fileIn, fileToExclude, fileOut):
    return f"plink --bfile {fileIn} --exclude {fileToExclude} --make-bed --out {fileOut}"


def invalid(field):
    NAs = ['', ' ', 'na', 'nan', '-1']
    if field in NAs:
        return True
    return False


def countInd(fileName, fileType):
    if fileType == "bed":
        count = 0
        file = open(fileName + ".fam")
        for line in file:
            count = count + 1
        return count


def countSNPs(fileName, fileType):
    if fileType == "bed":
        count = 0
        file = open(fileName + ".bim")
        for line in file:
            count = count + 1
        return count

def countFiles(fileName, type, file, typeData, step):
    print(f'Counting data from {fileName}')
    SNPs = countSNPs(fileName, type)
    ind = countInd(fileName, type)
    file.write(f'{typeData}\t{step}\t{SNPs}\t{ind}\n')

    return SNPs, ind

def removeStructural(fileName, structural, output, folder):
    structuralFile = open(structural)
    remove = open(f'{folder}/{output}_structural','w')

    for line in structuralFile:
        split = line.split()
        remove.write(f'{split[0]}\n')
    remove.close()

    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_structural', f'{folder}/{output}_NonStructural')
    execute(command)
    return f'{folder}/{output}_NonStructural', 'bed'

def readLine(line, fileType):
    split = line.strip().split()
    if fileType == "bed":
        return split[0], split[3], split[4], split[5], split[1]
    elif fileType == "vcf":
        return split[0], split[1], split[3], split[4], split[2], split[7]

def removeDuplicated(fileName, output, folder):
    file = open(f'{fileName}.bim')
    duplicate = open(f'{folder}/{output}_duplicateList', 'w')
    dictSNPs = {}
    toRemove = {}

    for line in file:
        chr, pos, a1, a2, ID = readLine(line, "bed")

        if chr not in dictSNPs:
            dictSNPs[chr] = {}

        if pos not in dictSNPs[chr]:
            dictSNPs[chr][pos] = []
            dictSNPs[chr][pos].append(ID)
        else:
            if ID not in dictSNPs[chr][pos]:
                dictSNPs[chr][pos].append(ID)

            if chr not in toRemove:
                toRemove[chr] = []

            if pos not in toRemove[chr]:
                toRemove[chr].append(pos)

    for chr in toRemove:
        for pos in toRemove[chr]:
            for ID in dictSNPs[chr][pos]:
                duplicate.write(f'{ID}\n')
    duplicate.close()

    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_duplicateList', f'{folder}/{output}_nonDuplicate')
    execute(command)

    return f'{folder}/{output}_nonDuplicate', 'bed'

def removePotentialProbe(fileName, output, folder, gnomad, chrX = False):
    file = open(f'{fileName}.bim')
    probe = open(f'{folder}/{output}_probe', 'w')
    reasonFile = open(f'{folder}/{output}_potentialProbeSite', 'w')

    toRemove = []
    reasonList = []

    chrList = []
    dictSNPs = {}
    for line in file:
        chrom, pos, a1, a2, ID = readLine(line, "bed")
        if chrX:
            chrom = "X"

        if chrom not in chrList:
            print(f'Reading the file {gnomad.replace("*", chrom)}')
            chrList.append(chrom)
            dictSNPs = putVCFInDict(gnomad.replace('*', chrom))
            print(f'Done')

        for i in range(-20, 21):
            if i != 0:
                position = int(pos) + i
                posStr = str(position)

                if posStr in dictSNPs:
                    probe.write(f'{ID}\n')
                    reasonFile.write(f'{ID} -> {chrom}:{posStr}\n')
                    break

    probe.close()
    reasonFile.close()
    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_probe', f'{folder}/{output}_probeSolved')
    execute(command)
    return f'{folder}/{output}_probeSolved', 'bed'



def putVCFInDict(vcfFile):
    vcfFile = gzip.open(vcfFile)
    dictSNPs = {}
    header = True
    for line in vcfFile:
        line = line.decode('utf-8')
        if header:
            if line[0:6] == "#CHROM":
                headerLine = line.strip().split()
                header = False
        else:
            split = line.strip().split()
            pos = split[1]
            data = split[-1].split("=")
            MAF = float(data[-1])

            if MAF > 0.5:
                MAF = 1 - MAF
            if MAF > 0.01:
                dictSNPs[pos] = MAF
    return dictSNPs





# def removePotentialProbe(fileName, output, folder, gnomad, chrX = False):
#     file = open(f'{fileName}.bim')
#     probe = open(f'{folder}/{output}_probe', 'w')
#     reasonFile = open(f'{folder}/{output}_potentialProbeSite', 'w')
#
#     toRemove = []
#     reasonList = []
#
#     dictSNPs = {}
#     for line in file:
#         chr, pos, a1, a2, ID = readLine(line, "bed")
#         if chrX:
#             chr = "X"
#
#         if chr not in dictSNPs:
#             dictSNPs[chr] = {}
#         dictSNPs[chr][pos] = ID
#
#     for chrom in dictSNPs:
#         start = time.time()
#         gnomadWithChr = gnomad.replace('*', chrom)
#         print(f'Looking gnomAD chr {chrom} => {gnomadWithChr}')
#         tb = tabix.open(gnomadWithChr)
#
#         for pos in dictSNPs[chrom]:
#             records = tb.query('chr'+chrom, int(pos)-20, int(pos)+20)
#             reason = getProbProblem(records)
#             if reason:
#                 toRemove.append(dictSNPs[chrom][pos])
#                 reasonList.append(reason)
#
#         end = time.time() - start
#         print(f'\tDone in {end} s')
#     for i in len(toRemove):
#         probe.write(f'{toRemove[i]}\n')
#         reasonFile.write(f'{toRemove[i]} -> {reasonList[i]}\n')
#     probe.close()
#     reasonFile.close()
#     command = lineToRemoveSNPs(fileName, f'{folder}/{output}_probe', f'{folder}/{output}_probeFix')
#     execute(command)
#     return f'{folder}/{output}_probeFix', "bed"
#
#
# def getProbProblem(records):
#     for rec in records:
#         info=rec[-1].strip().split('=')
#         print(info)
#         if float(info[1]) > (1/100):
#             return f'{rec[0]}:{rec[1]}'
#     return ""

def HWE(fileName, output, folder, group, pval, females=False):
    command = f"{plink} --bfile {fileName} --hardy --out {folder}/{output}_HWE"
    if females:
        command = command+" --filter-females"
    execute(command)

    file = open(f'{folder}/{output}_HWE.hwe')
    fileRemove = open(f'{folder}/{output}_FailedHWE', 'w')
    header = True
    for line in file:
        if header:
            header = False
        else:
            split = line.strip().split()
            ID = split[1]
            groupTested = split[2]
            if split[-1] != 'NA':
                pvalue = float(split[-1])
                if groupTested == group:
                    if pvalue < pval:
                        fileRemove.write(f'{ID}\n')
    fileRemove.close()
    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_FailedHWE', f'{folder}/{output}_HWE')
    print(f'{command}')
    execute(command)

    return f'{folder}/{output}_HWE', 'bed'

def genoMind(fileName, output, folder, cutoff):
    command = f"{plink} --bfile {fileName} --geno {cutoff} --mind {cutoff} --make-bed --out {folder}/{output}_genoMind"
    execute(command)

    return f'{folder}/{output}_genoMind', 'bed'

def pruneLD(fileName, output, folder):
    command = f"{plink} --bfile {fileName} --indep-pairwise 50 10 0.5 --out {folder}/{output}_LD"
    execute(command)
    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_LD.prune.out', f'{folder}/{output}_LDPruned')
    execute(command)
    return f'{folder}/{output}_LDPruned', 'bed'


def kinshipControl(fileName, output, folder, remove):
    command = f'{king} -b {fileName}.bed --kinship --prefix {folder}/{output}'
    execute(command)

    file = open(f'{folder}/{output}.kin0')
    fileOut = open(f'{folder}/{output}.NAToRA','w')

    header = True
    for line in file:
        if header :
            header = False
        else:
            split = line.strip().split()
            ID1 = split[1]
            ID2 = split[3]
            kinship = split[-1]
            fileOut.write(f'{ID1}\t{ID2}\t{kinship}\n')
    fileOut.close()

    command = f'{python} {NAToRA} -i {folder}/{output}.NAToRA -k -d {remove} -o {folder}/{output}_NAToRA'
    execute(command)


    file = open(f'{folder}/{output}_NAToRA_toRemove.txt')
    fileOut = open(f'{folder}/{output}_toRemovePLINK','w')

    for line in file:
        fileOut.write(f'{line.strip()} {line.strip()}\n')
    fileOut.close()

    command = lineToRemoveIndividual(fileName, f'{folder}/{output}_toRemovePLINK', f'{folder}/{output}_Relatedness')
    execute(command)

    return f'{folder}/{output}_Relatedness', "bed"

def removeMonomorphic(fileName, outputPrefix, folder):
    fileFam = open(f'{fileName}.fam')

    numInd = 0
    for line in fileFam:
        numInd = numInd +1

    freqNonMonomorphic = 1/((numInd*2)*10) #10x to have sure that just remove monomorphic

    command = f'{plink} --bfile {fileName} --maf {freqNonMonomorphic} --out {folder}/{outputPrefix}_nonMonomorphic --make-bed'
    execute(command)

    return f'{folder}/{outputPrefix}_nonMonomorphic', 'bed'

def sexCheck(fileName, output, folder, build):
    command = f'{plink} --bfile {fileName} --split-x {build} no-fail --make-bed --out {folder}/{output}_splitX '
    execute(command)
    command = f'{plink} --bfile {folder}/{output}_splitX --check-sex --out {folder}/{output}_sexCheck'
    execute(command)

    file = open(f'{folder}/{output}_sexCheck.sexcheck')
    fileRemove = open(f'{folder}/{output}_toRemoveCheck', 'w')
    for line in file:
        split = line.strip().split()
        if split[4] == "PROBLEM":
            if split[3] != "0":
              fileRemove.write(f'{split[0]} {split[1]}\n')
    fileRemove.close()

    command = lineToRemoveIndividual(f'{folder}/{output}_splitX', f'{folder}/{output}_toRemoveCheck', f'{folder}/{output}_sexChecked')
    execute(command)
    return f'{folder}/{output}_sexChecked', "bed"

def removeIndividualsOnChrX(fileName, autosomal, output, folder):
    autosomalFam = open(f'{autosomal}.fam')
    toKeep = open(f'{folder}/{output}_toKeep','w')

    for line in autosomalFam:
        split = line.split()
        toKeep.write(f'{split[0]} {split[1]}\n')
    toKeep.close()

    command = lineToRemoveIndividual(fileName, f'{folder}/{output}_toKeep', f'{folder}/{output}_Relationship_X')
    command = command.replace('--remove','--keep')
    execute(command)

    return f'{folder}/{output}_Relationship_X', "bed"

def hetInMales(fileName, output, folder):
    command = f'{plink} --bfile {fileName} --recode --out {folder}/{output}'
    execute(command)

    command = f'cp {folder}/{output}.map {folder}/{output}_Het.map'
    execute(command)

    pedFile = open(f'{folder}/{output}.ped')
    pedFileOut = open(f'{folder}/{output}_Het.ped', 'w')

    diff = 0
    for line in pedFile:
        split = line.strip().split()
        if split[4] == "1":
            pedFileOut.write(f'{split[0]} {split[1]} {split[2]} {split[3]} {split[4]} {split[5]}')
            for i in range(6,len(split),2):
                if split[i] != split[i+1]:
                    diff = diff + 1
                    pedFileOut.write(f' 0 0')
                else:
                    pedFileOut.write(f' {split[i]} {split[i+1]}')
            pedFileOut.write('\n')
        else:
            pedFileOut.write(line)

    pedFileOut.close()
    command = f'{plink} --file {folder}/{output}_Het --make-bed --out {folder}/{output}_HetMale'
    execute(command)

    print(f'We set 2*{diff} Alleles as missing due to heterozygous SNPs in males ')
    return f'{folder}/{output}_HetMale', 'bed'

def differentialMissingCaseAndControl(fileName, output, folder, cutoff):
    command = f'{plink} --bfile {fileName} --test-missing --out {folder}/{output}_differentialMissingness'
    execute(command)

    file = open(f'{folder}/{output}_differentialMissingness.missing')
    fileToRemove = open(f'{folder}/{output}_missingDiff', 'w')
    header = True
    for line in file:
        if header:
            header = False
        else:
            split = line.strip().split()
            if float(split[-1]) < cutoff:
                fileToRemove.write(split[1]+"\n")
    fileToRemove.close()
    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_missingDiff', f'{folder}/{output}_missingDifferential')
    execute(command)

    return f'{folder}/{output}_missingDifferential', "bed"

def differentialMissingMaleAndFemale(fileName, output, folder, cutoff):
    fileFam = open(f'{fileName}.fam')
    fileFamNew = open(f'{fileName}_DiffBetweenSex.fam','w')

    for line in fileFam:
        split = line.strip().split()
        fileFamNew.write(f'{split[0]} {split[1]} {split[2]} {split[3]} {split[4]} {split[4]}\n')
    fileFamNew.close()
    command = f'cp {fileName}.bed {fileName}_DiffBetweenSex.bed'
    execute(command)
    command = f'cp {fileName}.bim {fileName}_DiffBetweenSex.bim'
    execute(command)

    fileNameNew, type = differentialMissingCaseAndControl(f'{fileName}_DiffBetweenSex', output, folder, cutoff)
    command = f'cp {fileName}.fam {fileNameNew}.fam'
    execute(command)
    return fileNameNew, "bed"

def differentialMAFMaleAndFemale(fileName, output, folder, cutoff):
    command = f'{plink} --bfile {fileName} --freq counts --out {folder}/{output}_FreqMale --filter-males'
    execute(command)

    command = f'{plink} --bfile {fileName} --freq counts --out {folder}/{output}_FreqFemale --filter-females'
    execute(command)

    file = open(f'{folder}/{output}_FreqMale.frq.counts')
    maleDict = {}
    header = True
    for line in file:
        if header:
            header = False
        else:
            split = line.split()
            maleDict[split[1]] = line

    file = open(f'{folder}/{output}_FreqFemale.frq.counts')
    fileToRemove = open(f'{folder}/{output}_differentialMAF', 'w')
    header = True
    for line in file:
        if header:
            header = False
        else:
            splitFemale = line.strip().split()
            splitMale = maleDict[splitFemale[1]].strip().split()

            if splitFemale[3] == splitMale[3]:
                maleMajor = splitMale[5]
                maleMinor = splitMale[4]
                femaleMajor = splitFemale[5]
                femaleMinor = splitFemale[4]
            else:
                maleMajor = splitMale[4]
                maleMinor = splitMale[5]
                femaleMajor = splitFemale[5]
                femaleMinor = splitFemale[4]

            #Table
            #                       Female  Male
            # Major Allele Count
            # Minor Allele Count

            odds, pvalue = fisher_exact([[femaleMajor, maleMajor], [femaleMinor, maleMinor]], alternative='two-sided')

            if pvalue < cutoff:
                print(f'{splitMale[1]} F:{femaleMajor} {femaleMinor} -- M:{maleMajor} {maleMinor}')
                fileToRemove.write(f'{splitMale[1]}\n')

    fileToRemove.close()
    command = lineToRemoveSNPs(fileName, f'{folder}/{output}_differentialMAF', f'{folder}/{output}_diffMAF')
    execute(command)

    return f'{folder}/{output}_diffMAF', 'bed'

def flip(fileName, out, folder, oneThousand, threads, geneticMap, homozygous):
    sexDict = {}
    famFile = open(fileName+".fam")
    for line in famFile:
        split = line.split()
        sexDict[split[1]] = split[4]

    command = f'{plink} --bfile {fileName} --recode vcf-iid --out {folder}/{out}_toFixSex'
    execute(command)

    vcfFile = open(f'{folder}/{out}_toFixSex.vcf')
    fileOut = open(f'{folder}/{out}_sexFixed.vcf', 'w')
    header = True


    for line in vcfFile:
        if header:
            fileOut.write(line)
            if line[0:6] == "#CHROM":
                headerLine = line.strip().split()
                header = False
        else:
            split = line.strip().split()
            fileOut.write(f'chrX')
            for i in range(1, len(split)):
                if i < 9:
                    fileOut.write(f'\t{split[i]}')
                else:
                    sex = str(sexDict[headerLine[i]])
                    if sex == "1" or sex == "M":
                        if not homozygous:
                            fileOut.write(f'\t{split[i][0]}')
                        else:
                            fileOut.write(f'\t{split[i][0]}/{split[i][0]}')
                    else:
                        fileOut.write(f'\t{split[i]}')
            fileOut.write('\n')
    fileOut.close()


    command = f'{bgzip} {folder}/{out}_sexFixed.vcf'
    execute(command)
    command = f'{bcftools} index {folder}/{out}_sexFixed.vcf.gz'
    execute(command)

    vcfInput = f'{folder}/{out}_sexFixed.vcf.gz'

    command = f'{eagle} --vcfTarget {vcfInput} --vcfRef {oneThousand} --allowRefAltSwap --outPrefix {folder}/{out}_Phased ' \
              f'--numThreads {threads} --geneticMapFile {geneticMap} --keepMissingPloidyX'
    execute(command)


def chrXSteps(fileName, outputPrefix, folder, structural, logFile, covar, pheno, continuous, gnomad, autosomal, build):
    dataType = "ChrX"
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Begin')

    fileName, type = setFamToCovar(fileName, covar, pheno, continuous, outputPrefix, folder)
    fileName, fileType = removeIndividualsOnChrX(fileName, autosomal, outputPrefix, folder)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Covar and relatedness')

    # Looks like this is just chr X
    fileName, type = removeStructural(fileName, structural, outputPrefix, folder)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Structural')

    fileName, type = removeDuplicated(fileName, outputPrefix, folder)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Duplicated')

    fileName, type = removeMonomorphic(fileName, outputPrefix, folder)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Monomorphic')

    #Check
    if gnomad:
        fileName, type = removePotentialProbe(fileName, outputPrefix, folder, gnomad, True)

        SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Probe')

    fileName, type = HWE(fileName, outputPrefix, folder, "UNAFF", 0.00001)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'HWE')

    fileName, type = genoMind(fileName, outputPrefix, folder, 0.05)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Geno')

    fileName, type = sexCheck(fileName, outputPrefix, folder, build)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'SexCheck')

    fileName, type = hetInMales(fileName, outputPrefix, folder)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'HetMale')

    #fileName, type = differentialMissingCaseAndControl(fileName, outputPrefix, folder, 0.00001)
    #SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'DifferentialMissing')

    #fileName, type = HWE(fileName, outputPrefix, folder, "ALL", 0.00001, True)
    #SNPs, ind = countFiles(fileName, "bed",  logFile, dataType, 'HWE females')

    fileName, type = removeMonomorphic(fileName, outputPrefix, folder)
    SNPs, ind = countFiles(fileName, "bed",  logFile, dataType, 'Monomorphic')

    fileName, type = differentialMissingMaleAndFemale(fileName, outputPrefix, folder, 0.00001)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'DifferentialMissingMaleFemale')

    fileName, type = differentialMAFMaleAndFemale(fileName, outputPrefix, folder, 0.00001)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'DifferentialMAFMaleFemale')

    return fileName


def autosomalSteps(fileName, outputPrefix, folder, structural, logFile, covar, pheno, continuous, gnomad, remove):
    dataType = "Autosomal"
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Begin')

    fileName, type = setFamToCovar(fileName, covar, pheno, continuous, outputPrefix, folder)

    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Covar')

    #Looks like this is just chr X
    fileName, type = removeStructural(fileName, structural, outputPrefix, folder)

    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Structural')

    fileName, type = removeDuplicated(fileName, outputPrefix, folder)

    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Duplicated')

    fileName, type = removeMonomorphic(fileName, outputPrefix, folder)

    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Monomorphic')

    if gnomad:
        fileName, type = removePotentialProbe(fileName, outputPrefix, folder, gnomad)

        SNPs, ind = countFiles(fileName, "bed", logFile,"Autosomal", "Probe")

    fileName, type = HWE(fileName, outputPrefix, folder, "UNAFF", 0.00001)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'HWE')

    fileName, type = genoMind(fileName, outputPrefix, folder, 0.05)
    SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Geno')

    if remove:
        fileName, type = kinshipControl(fileName, outputPrefix, folder, remove)
        SNPs, ind = countFiles(fileName, "bed", logFile, dataType, 'Relationship')
    return fileName

def setPrograms(args):
    global plink
    global plink2
    global bcftools
    global eagle
    global NAToRA
    global king
    global python
    global bgzip


    if args.plink:
        plink = args.plink
    else:
        plink = "plink"

    if args.plink2:
        plink2 = args.plink2
    else:
        plink2 = "plink2"

    if args.bcftools:
        bcftools = args.bcftools
    else:
        bcftools = "bcftools"

    if args.NAToRA:
        NAToRA = args.NAToRA
    else:
        NAToRA = "NAToRA"

    if args.eagle:
        eagle = args.eagle
    else:
        eagle = "eagle"

    if args.king:
        king = args.king
    else:
        king = "king"

    if args.interpreter:
        python = args.interpreter
    else:
        python = "python"

    if args.bgzip:
        bgzip = args.bgzip
    else:
        bgzip = "bgzip"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='XWAS cleaning')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-a', '--autosomal', help='Autosomal input file', required=True)
    required.add_argument('-x', '--xchromosome', help='Chromosome X input file', required=True)
    required.add_argument('-s', '--structural', help='Variants in strutural region', required=True)
    required.add_argument('-o', '--output', help='Output file prefix', required=True)
    required.add_argument('-O', '--oneThousand', help='Path for One Thousand Genomes vcf file', required=True)
    required.add_argument('-f', '--folder', help='Folder to store output files', required=True)
    required.add_argument('-c', '--covar', help='File with the covariatives', required=True)
    required.add_argument('-g', '--geneticMap', help='Path to genetic Map from eagle', required=True)
    required.add_argument('-p', '--phenotype', help='Phenotype column on covar file', required=True)
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-I', '--interpreter', help='Path to python3 (default: python)', required=False, default="")
    optional.add_argument('-P', '--plink', help='Plink path (default: plink)', required=False, default="")
    optional.add_argument('-T', '--plink2', help='Plink2 path (default: plink2)', required=False, default="")
    optional.add_argument('-B', '--bcftools', help='BCFTools path (default: bcftools)', required=False, default="")
    optional.add_argument('-b', '--bgzip', help='Bgzip path (default: bgzip)', required=False, default="")
    optional.add_argument('-E', '--eagle', help='Eagle path (default: eagle)', required=False, default="")
    optional.add_argument('-N', '--NAToRA', help='NAToRA path (default: NAToRA)', required=False, default="")
    optional.add_argument('-K', '--king', help='KING path (defaul: king)', required=False, default="")
    optional.add_argument('-d', '--degree', help='Relationship degree to be removed. (default = "", that means'
                                                 'to not remove relationship) ',
                          required=False, default="")
    optional.add_argument('-C', '--continuous', help='Set the phenotype as continuos', required=False, default=False,
                          action="store_true")
    optional.add_argument('-G', '--gnomAD', help='Path to gnomAD files with the chr replaced by * (default = "", '
                                                 'not perform this step)', required=False, default="")
    optional.add_argument('-r', '--referenceBuild', help='Human Genome build (default : hg38)', required=False, default="hg38")
    optional.add_argument('-t', '--threads', help='Number of threads to Eagle (default = 1)', required=False, default="1")
    optional.add_argument('-H', '--homozygousMale', help='The males has homozygous diploid chr X (default False)', required=False,
                          default=False, action = "store_true")


    args = parser.parse_args()

    autosomalPrefix = args.autosomal
    chrXPrefix = args.xchromosome
    structural = args.structural
    continuous = args.continuous
    pheno = args.phenotype
    remove = args.degree
    output = args.output
    folder = args.folder
    gnomad = args.gnomAD
    covar = args.covar



    setPrograms(args)


    # Pre-cleaning
    command = f"mkdir {folder}"
    execute(command)
    outLog = open(f'{folder}/{output}.log', "w")
    outLog.write('Data\tStep\tN SNPs\tN Individuals\n')


    autosomalFiles = autosomalSteps(autosomalPrefix, output, folder, structural, outLog, covar, pheno, continuous, gnomad, remove)
    chrX = chrXSteps(chrXPrefix, output+"_chrX", folder, structural, outLog, covar, pheno, continuous, gnomad,
                     autosomalFiles, args.referenceBuild)

    flip(chrX, output+"_chrX", folder, args.oneThousand, args.threads, args.geneticMap, args.homozygousMale)
    outLog.close()
