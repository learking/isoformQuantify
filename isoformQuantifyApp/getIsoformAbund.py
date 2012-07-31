from django.http import HttpResponse
from django.utils import simplejson
import json
import subprocess
import os

from scipy.optimize import nnls
from models import *

debug_flag = False

def getReads(gene):
    geneID = gene['id']
    geneChr = gene['chr']
    geneStart = gene['start']
    geneEnd = gene['end']

    samtoolsFile = '/home/s-kwang/software/samtools-0.1.18/samtools'
    samtoolsOption = 'view'
    samtoolsOutput = '-o'
    samFile = '/home/s-kwang/webApp/isoformQuantify/isoformQuantifyApp/data/sam/' + geneID + '.sam'
    seqFile = '/home/s-kwang/webApp/isoformQuantify/isoformQuantifyApp/data/bam/cno.bam'
    chrRegion = 'chr' + geneChr + ':' + str(geneStart) + '-' + str(geneEnd)

    reads = Reads()
    if not os.path.exists(samFile):
        process = subprocess.Popen([samtoolsFile, samtoolsOption, samtoolsOutput, samFile, seqFile, chrRegion])
        #for line in open('./isoformQuantifyApp/data/Apoe.sam'):
        if process.wait() == 0:
            for line in open(samFile):
                line = line.rstrip()
                parts = line.split('\t')
                reads.addRead(ShortRead(parts))
    else:
        for line in open(samFile):
            line = line.rstrip()
            parts = line.split('\t')
            reads.addRead(ShortRead(parts))

    reads.seperateReads()
    reads.initializeLength()    
    return reads

def getGeneExons(geneExons):
    exons = []
    for exon in geneExons:
        exons.append(geneExon(exon['start'], exon['end']))
    return exons

def getGeneTranscripts(geneTranscripts):
    transcripts = []
    for transcript in geneTranscripts:
        transcripts.append(geneTranscript(transcript))
    return transcripts

def getExonBounds(exons):
    exonBounds = []
    for exon in exons:
        if exon.start not in exonBounds:
            exonBounds.append(exon.start)
        if exon.end not in exonBounds:
            exonBounds.append(exon.end)
    exonBounds.sort()
    return exonBounds


def getSubExons(exons):
    subExons = []
    exonBounds = getExonBounds(exons)
    for bound_index in range( len(exonBounds) - 1 ):
        isSubExon = False
        nextBound_index = bound_index + 1
        currentSubExon = SubExon(exonBounds[bound_index],exonBounds[nextBound_index])
        for exon in exons:
            if currentSubExon.start >=  exon.start and currentSubExon.end <= exon.end:
                isSubExon = True
                exon.subExons.append(len(subExons))
        if isSubExon:
            subExons.append(currentSubExon)
    return subExons

'''
def initialize_eachTranscript(readLength):
    for transcript in targetGene.transcripts:
        currentTranscript = targetGene.transcripts[transcript]
        initialize_subexonForTranscript(currentTranscript)
        currentTranscript.initialize_length()
        currentTranscript.calculate_subExonPercent()
        currentTranscript.initialize_spliceJunctions(readLength)
        currentTranscript.initialize_sudo_spliceJunctions(readLength)    
'''
def iniTranscripts(transcripts, exons, subExons):
    for transcript in transcripts:
        for exon in transcript.exons:
            transcript.length += exons[exon].length
            for subExon in exons[exon].subExons:
                transcript.subExons.append(subExon)
        transcript.subExons.sort()

        for currentIndex, subExon in enumerate(transcript.subExons):
            transcript.subExonPercents.append(float(subExons[subExon].length) / float(transcript.length))

            nextIndex = currentIndex + 1
            if nextIndex < len(transcript.subExons):
                currentJunction = (transcript.subExons[currentIndex], transcript.subExons[nextIndex])
                if subExons[transcript.subExons[currentIndex]].end != subExons[transcript.subExons[nextIndex]].start :
                    transcript.spliceJunctions[currentJunction] = spliceJunction()
                else:
                    transcript.sudoSpliceJunctions[currentJunction] = sudoSpliceJunction()

def iniJunctions(transcripts, exons, subExons, readLength):
    for transcript in transcripts:
        if transcript.sudoSpliceJunctions:
            for key in transcript.sudoSpliceJunctions:
                sudoSpliceJunction = transcript.sudoSpliceJunctions[key]
                for exon in transcript.exons:
                    if key[0] in exons[exon].subExons:
                        sudoSpliceJunction.middle = subExons[key[0]].end
                        sudoSpliceJunction.left = max(exons[exon].start, sudoSpliceJunction.middle - readLength + 2)
                        sudoSpliceJunction.right = min(exons[exon].end, sudoSpliceJunction.middle + readLength - 2)
                        sudoSpliceJunction.length = sudoSpliceJunction.right - sudoSpliceJunction.left + 1
                        sudoSpliceJunction.ratio = float(sudoSpliceJunction.length)/float(transcript.length)
        if transcript.spliceJunctions:
            for key in transcript.spliceJunctions:
                spliceJunction = transcript.spliceJunctions[key]
                spliceJunction.leftEnd = subExons[key[0]].end
                spliceJunction.rightStart = subExons[key[1]].start
                leftLength = 0
                for exon in transcript.exons:
                    if exons[exon].end <= spliceJunction.leftEnd:
                        leftLength += exons[exon].length
                rightLength = transcript.length - leftLength
                spliceJunction.leftLength = min(readLength - 3, leftLength)
                spliceJunction.rightLength = min(readLength - 3, rightLength)
                spliceJunction.length = spliceJunction.leftLength + spliceJunction.rightLength
                spliceJunction.ratio = float(spliceJunction.length) / float(transcript.length)

def count_readsWithinSubexons(subExons, reads):
    readLength = reads.length
    for subExon in subExons:
        for read in reads.wholeReads:
            if int(read.start) >= int(subExon.start) and int(read.start) + (readLength - 1) <= int(subExon.end):
                subExon.readCounts += 1

def getSpliceJunctions(transcripts):
    spliceJunctions = []
    for transcript in transcripts:
        for spliceJunction in transcript.spliceJunctions:
            if not (spliceJunction in spliceJunctions):
                spliceJunctions.append(spliceJunction)
    return spliceJunctions

def getSudoSpliceJunctions(transcripts):
    sudoSpliceJunctions = []
    for transcript in transcripts:
        for sudoSpliceJunction in transcript.sudoSpliceJunctions:
            if not (sudoSpliceJunction in sudoSpliceJunctions):
                sudoSpliceJunctions.append(sudoSpliceJunction)
    return sudoSpliceJunctions

def count_spliceJunction_reads(spliceJunctions, subExons, reads):
    spliceJunctionReadCounts = []
    for spliceJunction in spliceJunctions:
        readCount = 0
        splicePosition = (subExons[spliceJunction[0]].end , subExons[spliceJunction[1]].start)
        for spliceRead in reads.splicedReads:
            if spliceRead.splicePositions[0] == splicePosition:
                readCount += 1
            if len(spliceRead.splicePositions) == 2:
                if spliceRead.splicePositions[1] == splicePosition:
                    readCount += 1
        spliceJunctionReadCounts.append(readCount)
    return spliceJunctionReadCounts

def count_sudoSpliceJunction_allTranscripts(transcripts, reads):
    readLength = reads.length
    for transcript in transcripts:
        for sudoSpliceJunction in transcript.sudoSpliceJunctions:
            currentSudoSpliceJunction = transcript.sudoSpliceJunctions[sudoSpliceJunction]
            for read in reads.wholeReads:
                if int(read.start) >= int(currentSudoSpliceJunction.left) and int(read.start) + (readLength - 1) <= int(currentSudoSpliceJunction.right):
                    currentSudoSpliceJunction.readCount += 1                

def get_sudoSpliceJunction_MaxReads(sudoSpliceJunctions, transcripts):
    sudoSpliceJunctionReadCounts = []
    for sudoSpliceJunction in sudoSpliceJunctions:
        maxReadCount = 0
        readCount = 0
        for transcript in transcripts:
            if sudoSpliceJunction in transcript.sudoSpliceJunctions.keys():
                readCount = transcript.sudoSpliceJunctions[sudoSpliceJunction].readCount
                if readCount > maxReadCount:
                    maxReadCount = readCount
        sudoSpliceJunctionReadCounts.append(maxReadCount)
    return sudoSpliceJunctionReadCounts

def returnIsoformAbund(request):
    if(debug_flag):
        f = open('./isoformQuantifyApp/data/log.txt', 'w')

    gene = simplejson.loads(request.body)
    exons = getGeneExons(gene['exons'])
    subExons = getSubExons(exons)
    transcripts = getGeneTranscripts(gene['transcripts'])

    reads = getReads(gene)
    readLength = reads.length

    count_readsWithinSubexons(subExons, reads)
    iniTranscripts(transcripts, exons, subExons)
    iniJunctions(transcripts, exons, subExons, readLength)
    spliceJunctions = getSpliceJunctions(transcripts)
    sudoSpliceJunctions = getSudoSpliceJunctions(transcripts)

    spliceJunctionReadCounts = count_spliceJunction_reads(spliceJunctions, subExons, reads)
    count_sudoSpliceJunction_allTranscripts(transcripts, reads)
    sudoSpliceJunctionReadCounts = get_sudoSpliceJunction_MaxReads(sudoSpliceJunctions, transcripts)
    
    designMatrix = []
    response = []
    for subExon_index, subExon in enumerate(subExons):
        if subExon.length >= readLength:
            currentFeature = []
            for transcript in transcripts:
                if subExon_index in transcript.subExons:
                    subExonPosition = transcript.subExons.index(subExon_index)
                    currentFeature.append(transcript.subExonPercents[subExonPosition])
                else:
                    currentFeature.append(0)
            if(debug_flag):
                f.write(str(subExon_index) + ":length:" + str(subExon.length) + ":" + str(currentFeature) + " readCount:" + str(subExon.readCounts) + "\n")
            designMatrix.append(currentFeature)
            response.append(subExon.readCounts)

    for spliceJunction_index, spliceJunction in enumerate(spliceJunctions):
        if spliceJunctionReadCounts[spliceJunction_index] > 0:
            currentFeature = []
            for transcript in transcripts:
                if spliceJunction in transcript.spliceJunctions.keys():
                    currentFeature.append(transcript.spliceJunctions[spliceJunction].ratio)
                else:
                    currentFeature.append(0)
            if(debug_flag):
                f.write(str(currentFeature) + " readCount:" + str(spliceJunctionReadCounts[spliceJunction_index]) + "\n")
            designMatrix.append(currentFeature)
            response.append(spliceJunctionReadCounts[spliceJunction_index])

    for sudoSpliceJunction_index, sudoSpliceJunction in enumerate(sudoSpliceJunctions):
        if sudoSpliceJunctionReadCounts[sudoSpliceJunction_index] > 0:
            currentFeature = []
            for transcript in transcripts:
                if sudoSpliceJunction in transcript.sudoSpliceJunctions.keys():
                    currentFeature.append(transcript.sudoSpliceJunctions[sudoSpliceJunction].ratio)
                else:
                    currentFeature.append(0)
            if(debug_flag):
                f.write(str(currentFeature) + " readCount:" + str(sudoSpliceJunctionReadCounts[sudoSpliceJunction_index]) + "\n")
            designMatrix.append(currentFeature)
            response.append(sudoSpliceJunctionReadCounts[sudoSpliceJunction_index])        

    result = nnls(designMatrix, response)[0]
    normalized_result = [x/sum(result) for x in result ]
    if(debug_flag):
        f.write(str(result) + "\n")

    '''
    for transcript in transcripts:
        f.write(transcript.id + " length:"+ str(transcript.length) +" exons: " + str(transcript.exons) + " subexons: " + str(transcript.subExons) + "\n")
        #f.write(str(transcript.subExonPercents) + " sum:" + str(sum(transcript.subExonPercents)) + "\n")
        f.write("spliceJunction: ")
        for key in transcript.spliceJunctions:
            f.write("||" + str(key) + " length:" + str(transcript.spliceJunctions[key].length)
 + " ratio:" + str(transcript.spliceJunctions[key].ratio))
        f.write("\n")
        f.write("sudoSpliceJunction: ")
        for key in transcript.sudoSpliceJunctions:
            f.write("||" + str(key) + " length:" + str(transcript.sudoSpliceJunctions[key].length) + " ratio:"+ str(transcript.sudoSpliceJunctions[key].ratio)  + " reads:" + str(transcript.sudoSpliceJunctions[key].readCount))
        f.write("\n")

    for exon_index, exon in enumerate(exons):
        f.write("exon " + str(exon_index) + " length:"+ str(exon.length) + " :" + str(exon)+" ")
        f.write(str(exon.subExons))
        f.write('\n')

    for subexon_index, subexon in enumerate(subExons):
        f.write("subexon " + str(subexon_index) + " :" + str(subexon))
        f.write('\n')

    f.write(str(spliceJunctions) + "\n")
    f.write(str(sudoSpliceJunctions) + "\n")
    f.write(str(spliceJunctionReadCounts) + "\n")
    f.write(str(sudoSpliceJunctionReadCounts) + "\n")
    '''
    if(debug_flag):
        f.close()

    resultMapping = {}
    for transcript_index, transcript in enumerate(transcripts):
        resultMapping['transcript_'+str(transcript_index)] = normalized_result[transcript_index]

    return HttpResponse(json.dumps(resultMapping))
