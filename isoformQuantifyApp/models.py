from django.db import models
import re

# Create your models here.
class Gene:

    def __init__(self, targetGeneID):
        self.id = targetGeneID
        self.chr = ""
        self.strand = ""
        self.exons = []
        self.transcripts = []
        self.start = 0
        self.end = 0

    def isTranscriptExist(self, transcriptID, transcripts):
        existFlag = False
        for transcript in transcripts:
            if transcript.id == transcriptID:
                existFlag = True
        return existFlag

    def getTranscripts(self, inputData):
        transcriptP = re.compile('transcript_id\s"(\w+)"')
        transcripts = []
        for line in inputData:
            transcriptID =  transcriptP.search(line).group(1)
            if not self.isTranscriptExist(transcriptID, transcripts):
                transcripts.append(Transcript(transcriptID))
        return transcripts

    def isExonExist(self, exon, exons):
        existFlag = False
        for currentExon in exons:
            if currentExon == exon:
                existFlag = True
        return existFlag

    def getExons(self, inputData):
        exons = []
        for line in inputData:
            parts = line.split('\t')
            exon = Exon(parts)
            if not self.isExonExist(exon, exons):
                exons.append(exon)
        return exons
    
    def getTranscriptByID(self, transcriptID):
        for transcript in self.transcripts:
            if transcript.id == transcriptID:
                return transcript

    def initializeTranscripts(self, inputData):
        transcriptP = re.compile('transcript_id\s"(\w+)"')
        for line in inputData:
            transcriptID = transcriptP.search(line).group(1)
            tmpTranscript = self.getTranscriptByID(transcriptID)
            parts = line.split('\t')
            tmpExon = Exon(parts)
            for exon_index, exon in enumerate(self.exons):
                if exon == tmpExon:
                    tmpTranscript.exons.append(exon_index)

    def getStrand(self, inputData):
        exonNumberP = re.compile('exon_number\s"(\d+)"')
        firstExonNumber = exonNumberP.search(inputData[0]).group(1)
        secondExonNumber = exonNumberP.search(inputData[1]).group(1)
        if firstExonNumber == "1" and secondExonNumber == "2":
            firstExonStart = inputData[0].split()[3]
            secondExonStart = inputData[1].split()[3]
            if firstExonStart < secondExonStart:
                return "+"
            else:
                return "-"
        else:
            print "ERROR: the first transcript only has one exon!"

    def initializeGeneStartEnd(self):
        for exon_index, exon in enumerate(self.exons):
            if exon_index == 0:
                self.start = int(exon.start)
                self.end = int(exon.end)
            if self.start > int(exon.start):
                self.start = int(exon.start)
            if self.end < int(exon.end):
                self.end = int(exon.end)

    def initializeGene(self, inputData):
        if len(inputData) == 0:
            print "ERROR: input file empty!"
        self.chr = inputData[0].split()[0]
        self.strand = self.getStrand(inputData)
        self.transcripts = self.getTranscripts(inputData)
        self.exons = self.getExons(inputData)
        self.initializeTranscripts(inputData)
        self.initializeGeneStartEnd()
        '''
        for transcript in self.transcripts:
            print transcript.id
            print transcript.exons
            for exon_index in transcript.exons:
                print str(self.exons[exon_index])
        '''

    def exons2Dict(self):
        exons = []
        for exon in self.exons:
            exons.append(exon.toDict())
        return exons

    def transcripts2Dict(self):
        transcripts = []
        for transcript in self.transcripts:
            transcripts.append(transcript.toDict())
        return transcripts

    def toDict(self):
        return {
            'id' : self.id,
            'chr': self.chr,
            'strand': self.strand,
            'start': self.start,
            'end' : self.end,
            'exons': self.exons2Dict(),
            'transcripts': self.transcripts2Dict(),
        }

class Transcript:
    def __init__(self, transcriptID):
        self.id = transcriptID
        self.exons = []

    def toDict(self):
        return {
            'id' : self.id,
            'exons' : self.exons,
        }

class Exon:
    def __init__(self, parts):
        self.start = int(parts[3])
        self.end = int(parts[4])

    def toDict(self):
        return {
            'start': self.start,
            'end': self.end,
        }

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end:
            return True
        else:
            return False

    def __str__(self):
        return '[%d %d]' % (self.start, self.end)

class ShortRead:
    
    def __init__(self, parts):
        self.start = parts[3]
        self.cigar = parts[5]

    def isWholeRead(self):
	return re.search("^\d+[a-zA-Z]+$", self.cigar)

    def __str__(self):
        return '[%s, %s]' % (self.start, self.cigar) 

class splicedRead(ShortRead):

    def __init__(self, shortRead):
	self.maxLength = 5
	self.start = shortRead.start
	self.cigar = shortRead.cigar
	self.cigarNumbers = re.findall(r"\d+" , self.cigar)
	self.readPositions = [int(x) for x in self.cigarNumbers]

	if len(self.readPositions) > self.maxLength:
		print "ERROR: exceed max length ", str(self)

	self.splicePositions = []
	self.initialize_splicePositions()
	#print self.splicePositions

    def initialize_splicePositions(self):
	if len(self.readPositions) < 3:
		print "ERROR: this is not a spliced read!"

	if len(self.readPositions) == 3:
		self.splicePositions = [(int(self.start) + self.readPositions[0] - 1,  int(self.start) + sum(self.readPositions[0:2]) )]

	if len(self.readPositions) == 5:
		self.splicePositions = [(int(self.start) + self.readPositions[0] - 1,  int(self.start) + sum(self.readPositions[0:2]) ), (int(self.start) + sum(self.readPositions[0:3]) - 1,  int(self.start) + sum(self.readPositions[0:4]) )]


class Reads:

    def __init__(self):
        self.reads = []
	self.wholeReads = []
	self.splicedReads = []
	self.length = 0

    def addRead(self, ShortRead):
        self.reads.append(ShortRead)

    def __str__(self):
        allReads = ''
	allReads += "whole reads:\n"
        for read in self.wholeReads:
            allReads += str(read)
	allReads += "spliced reads:\n"
        for read in self.splicedReads:
            allReads += str(read)
        return allReads

    def seperateReads(self):
	for read in self.reads:
		if read.isWholeRead():
			self.wholeReads.append(read)
		else:
			self.splicedReads.append(splicedRead(read))

    def initializeLength(self):
	if len(self.wholeReads) == 0:
		print "there is not whole read in .sam file!"
	else:
	   self.length = int(re.search("^\d+",self.wholeReads[0].cigar).group(0))

class geneExon:

    def __init__(self, startPosition, endPosition):
        self.start = startPosition
        self.end = endPosition
        self.length = self.end - self.start + 1
        self.subExons = []
        
    def __str__(self):
        return '[%d %d]' % (self.start, self.end)

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end:
            return True
        else:
            return False

class geneTranscript:

    def __init__(self, transcript):
        self.id = transcript['id']
        self.exons = transcript['exons']
        self.length = 0
        self.subExons = []
        self.subExonPercents = []
        self.spliceJunctions = {}
        self.sudoSpliceJunctions = {}

class SubExon():

    def __init__(self, startPosition, endPosition):
        self.start = startPosition
        self.end = endPosition
        self.length = self.end - self.start + 1
        self.readCounts = 0

    def __str__(self):
        return '[%d %d]' % (self.start, self.end)

class spliceJunction:

    def __init__(self):
        self.leftLength = 0
        self.leftEnd = 0
        self.rightStart = 0
        self.rightLength = 0
        self.length = 0
        self.ratio = 0

class sudoSpliceJunction:

    def __init__(self):
        self.left = 0
        self.middle = 0
        self.right = 0
        self.length = 0
        self.ratio = 0
        self.readCount = 0
