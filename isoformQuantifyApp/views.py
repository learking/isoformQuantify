# Create your views here.
from django.shortcuts import render_to_response
from django.http import HttpResponse
from models import *
import json
import subprocess
import re
import os

json_content_type = 'application/json'

def index(request):
    params = {
        'STATIC_URL' : '/static/',
    }
    return render_to_response('index.html', params)

def writeGTF(fileName, geneID):
    command_process = subprocess.Popen(['grep', geneID, './isoformQuantifyApp/data/gtf/Mus_musculus.NCBIM37.67.exons.gtf'], stdout=subprocess.PIPE)
    command_output = command_process.communicate()[0]

    if command_output != "":
        FILE = open(fileName,"w")
        FILE.write(command_output)
        FILE.close()

def initializeGeneJSON(request, geneID):
    fileName = "./isoformQuantifyApp/data/gtf/" + geneID + ".gtf"

    if not os.path.exists(fileName):
        writeGTF(fileName, geneID)

    try:
        inputFile = open(fileName,"r")
        inputData = inputFile.readlines()
    except IOError:
        print "File open failed!"

    targetGene = Gene(geneID)
    targetGene.initializeGene(inputData)

    if os.path.exists(fileName):
        return HttpResponse(json.dumps(targetGene.toDict()), content_type = json_content_type)
    else:
        return HttpResponse(json.dumps({'fileExist': 0}), content_type = json_content_type)
