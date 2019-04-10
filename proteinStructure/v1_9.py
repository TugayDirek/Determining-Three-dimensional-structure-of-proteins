import numpy as np
import os,glob
#import vpython
import math
#from numba import jit
import time
# Import pairwise2 module
from Bio import pairwise2
# Import format_alignment method
from Bio.pairwise2 import format_alignment
from scipy import spatial
from scipy.cluster.hierarchy import average, fcluster
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import pandas as pd





start_time = time.time()

data = []
i=0
counter = 0
j=-1
#counter_3 = False


#read data from pdb files
def readFromFile(data):
    i = 0
    counter = 0
    j = -1
    for filename in glob.glob(os.path.join('*.pdb')):
        j = j + 1
        data.append([filename])
        i = 0
        counter = 0
        #counter_3 = False
        with open(filename,'r') as f:
            for line in f:

                a = i
                counter = 0

                if line.__contains__("ENDMDL"):
                    break

                if line.__contains__("ATOM") and line.__contains__(" CA "):
                    data[j].append(["CA"])
                    i = i + 1
                    data[j][a + 1].append(line.split()[1])
                    data[j][a + 1].append(line.split()[3])
                    #x[j][a + 1].append(line.split()[6])
                    #x[j][a + 1].append(line.split()[7])
                    #x[j][a + 1].append(line.split()[8])
                    data[j][a + 1].append(line[30:38].strip(" "))
                    data[j][a + 1].append(line[38:46].strip(" "))
                    data[j][a + 1].append(line[46:54].strip(" "))


#print(x)

#delete element with less than three points
def deletePoorData(data):
    for point in range(0,len(data)):
        if(len(data[point])<=3):
            del(data[point])


angles = []
plane_angles = []
vectors = []
#angles_2 = []

i=-1

#calculate angles between vectors and planes
def calculateAngles(data):
    i = -1
    for element in data:
        i+=1
        plane_angles.append([element[0]])
        angles.append([element[0]])
        for point in range(0, len(element)):

            if point == len(element) - 3:
                break
            # print([element[0]],(element[point + 3][3]), (element[point + 3][4]), (element[point + 3][5]))
            #try:
            a = np.array([float(element[point + 1][3]), float(element[point + 1][4]), float(element[point + 1][5])])
            b = np.array([float(element[point + 2][3]), float(element[point + 2][4]), float(element[point + 2][5])])
            c = np.array([float(element[point + 3][3]), float(element[point + 3][4]), float(element[point + 3][5])])
            #except:
            #    continue

            ba = a - b
            bc = c - b

            cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            # print(cosine_angle)
            angle = np.arccos(cosine_angle)

            angles[i].append([np.degrees(angle)])

        for point in range(0,len(element)):

            if point == len(element)-4:
                break
            #print([element[0]],(element[point + 3][3]), (element[point + 3][4]), (element[point + 3][5]))
            #try:
            a = np.array([float(element[point+1][3]), float(element[point+1][4]), float(element[point+1][5])])
            b = np.array([float(element[point+2][3]), float(element[point+2][4]), float(element[point+2][5])])
            c = np.array([float(element[point+3][3]), float(element[point+3][4]), float(element[point+3][5])])
            d = np.array([float(element[point + 4][3]), float(element[point + 4][4]), float(element[point + 4][5])])
            #except:
            #    continue


            ba = a - b
            bc = c - b

            cb = b - c
            cd = d - c

            first_normal = np.cross(ba,bc)
            second_normal = np.cross(cb,cd)

            cosine_angle_2 = np.dot(first_normal, second_normal) / (np.linalg.norm(first_normal) * np.linalg.norm(second_normal))
            angle_2 = np.arccos(cosine_angle_2)


            #print(np.degrees(angle))
            plane_angles[i].append([np.degrees(angle_2)])



aadict = {
"GLU":"1 0 0 0 0","AGLU":"1 0 0 0 0","BGLU":"1 0 0 0 0","CGLU":"1 0 0 0 0","DGLU":"1 0 0 0 0",
"LEU":"0 1 0 0 0","ALEU":"0 1 0 0 0","BLEU":"0 1 0 0 0","CLEU":"0 1 0 0 0","DLEU":"0 1 0 0 0",
"THR":"0 0 1 0 0","ATHR":"0 0 1 0 0","BTHR":"0 0 1 0 0","CTHR":"0 0 1 0 0","DTHR":"0 0 1 0 0",
"PRO":"0 0 0 1 0","APRO":"0 0 0 1 0","BPRO":"0 0 0 1 0","CPRO":"0 0 0 1 0","DPRO":"0 0 0 1 0",
"ASP":"0 0 0 0 1","AASP":"0 0 0 0 1","BASP":"0 0 0 0 1","CASP":"0 0 0 0 1","DASP":"0 0 0 0 1",
"GLN":"1 1 0 0 0","AGLN":"1 1 0 0 0","BGLN":"1 1 0 0 0","CGLN":"1 1 0 0 0","DGLN":"1 1 0 0 0",
"HIS":"1 0 1 0 0","AHIS":"1 0 1 0 0","BHIS":"1 0 1 0 0","CHIS":"1 0 1 0 0","DHIS":"1 0 1 0 0",
"PHE":"1 0 0 1 0","APHE":"1 0 0 1 0","BPHE":"1 0 0 1 0","CPHE":"1 0 0 1 0","DPHE":"1 0 0 1 0",
"ILE":"1 0 0 0 1","AILE":"1 0 0 0 1","BILE":"1 0 0 0 1","CILE":"1 0 0 0 1","DILE":"1 0 0 0 1",
"MET":"1 1 0 1 0","AMET":"1 1 0 1 0","BMET":"1 1 0 1 0","CMET":"1 1 0 1 0","DMET":"1 1 0 1 0",
"SER":"1 1 0 0 1","ASER":"1 1 0 0 1","BSER":"1 1 0 0 1","CSER":"1 1 0 0 1","DSER":"1 1 0 0 1",
"TYR":"1 1 1 0 0","ATYR":"1 1 1 0 0","BTYR":"1 1 1 0 0","CTYR":"1 1 1 0 0","DTYR":"1 1 1 0 0",
"ASN":"1 1 1 1 0","AASN":"1 1 1 1 0","BASN":"1 1 1 1 0","CASN":"1 1 1 1 0","DASN":"1 1 1 1 0",
"LYS":"1 1 1 0 1","ALYS":"1 1 1 0 1","BLYS":"1 1 1 0 1","CLYS":"1 1 1 0 1","DLYS":"1 1 1 0 1",
"ARG":"1 1 1 1 1","AARG":"1 1 1 1 1","BARG":"1 1 1 1 1","CARG":"1 1 1 1 1","DARG":"1 1 1 1 1",
"ALA":"0 1 1 0 0","AALA":"0 1 1 0 0","BALA":"0 1 1 0 0","CALA":"0 1 1 0 0","DALA":"0 1 1 0 0",
"VAL":"0 1 1 1 0","AVAL":"0 1 1 1 0","BVAL":"0 1 1 1 0","CVAL":"0 1 1 1 0","DVAL":"0 1 1 1 0",
"GLY":"0 1 1 0 1","AGLY":"0 1 1 0 1","BGLY":"0 1 1 0 1","CGLY":"0 1 1 0 1","DGLY":"0 1 1 0 1",
"TRP":"0 1 1 1 1","ATRP":"0 1 1 1 1","BTRP":"0 1 1 1 1","CTRP":"0 1 1 1 1","DTRP":"0 1 1 1 1",
"CYS":"0 0 1 1 0","ACYS":"0 0 1 1 0","BCYS":"0 0 1 1 0","CCYS":"0 0 1 1 0","DCYS":"0 0 1 1 0"
}


pwdict = {
"GLU":"E","AGLU":"E","BGLU":"E","CGLU":"E","DGLU":"E",
"LEU":"L","ALEU":"L","BLEU":"L","CLEU":"L","DLEU":"L",
"THR":"T","ATHR":"T","BTHR":"T","CTHR":"T","DTHR":"T",
"PRO":"P","APRO":"P","BPRO":"P","CPRO":"P","DPRO":"P",
"ASP":"D","AASP":"D","BASP":"D","CASP":"D","DASP":"D",
"GLN":"Q","AGLN":"Q","BGLN":"Q","CGLN":"Q","DGLN":"Q",
"HIS":"H","AHIS":"H","BHIS":"H","CHIS":"H","DHIS":"H",
"PHE":"F","APHE":"F","BPHE":"F","CPHE":"F","DPHE":"F",
"ILE":"I","AILE":"I","BILE":"I","CILE":"I","DILE":"I",
"MET":"M","AMET":"M","BMET":"M","CMET":"M","DMET":"M",
"SER":"S","ASER":"S","BSER":"S","CSER":"S","DSER":"S",
"TYR":"Y","ATYR":"Y","BTYR":"Y","CTYR":"Y","DTYR":"Y",
"ASN":"N","AASN":"N","BASN":"N","CASN":"N","DASN":"N",
"LYS":"K","ALYS":"K","BLYS":"K","CLYS":"K","DLYS":"K",
"ARG":"R","AARG":"R","BARG":"R","CARG":"R","DARG":"R",
"ALA":"A","AALA":"A","BALA":"A","CALA":"A","DALA":"A",
"VAL":"V","AVAL":"V","BVAL":"V","CVAL":"V","DVAL":"V",
"GLY":"G","AGLY":"G","BGLY":"G","CGLY":"G","DGLY":"G",
"TRP":"W","ATRP":"W","BTRP":"W","CTRP":"W","DTRP":"W",
"CYS":"C","ACYS":"C","BCYS":"C","CCYS":"C","DCYS":"C"
}

encodedAaNames = []

i=-1



#insert encoded aminoacid codes
def insertEncodedNames(data):
    i = -1
    for element in data:
        i+=1
        encodedAaNames.append([element[0]])
        for point in range(0,len(element)):

            if point == len(element)-4:
                break

            #print(element)

            a = ((element[point+1][2]))
            b = ((element[point+2][2]))
            c = ((element[point+3][2]))
            d = ((element[point + 4][2]))

            #encodedAaNames[i].append([aadict.get(a),aadict.get(b),aadict.get(c),a_4.mag,b_4.mag,c_4.mag])
            encodedAaNames[i].append([aadict.get(a),aadict.get(b),aadict.get(c),aadict.get(d)])
        #print(encodedAaNames[i][0])


pairwiseAaNames = []
i = -1
#encode 3 letters aminoacid names to 1 letter
def pairwiseAlgnCode(data):
    i = -1
    for element in data:
        i+=1
        sequence = ""
        pairwiseAaNames.append([element[0]])
        for point in range(0,len(element)):

            if point == len(element)-1:
                break

            a = ((element[point+1][2]))
            sequence += pwdict.get(a)


        pairwiseAaNames[i].append(sequence)

pairwiseResultMatrix =  []
i = -1
j=-1

#do pairwise alignment and construct similarity matrix
def pairwiseAlignment(pairwiseAaNames):
    i = -1
    j = -1
    for element in range(0, len(pairwiseAaNames)):
        i += 1

        pairwiseResultMatrix.append([pairwiseAaNames[element][0]])
        X = pairwiseAaNames[element][1]

        for element_2 in range(0, len(pairwiseAaNames)):

            Y = pairwiseAaNames[element_2][1]

            alignments = pairwise2.align.globalms(X, Y, 2, -1, -0.5, -0.1,one_alignment_only=True)

            pairwiseResultMatrix[i].append([pairwiseAaNames[element_2][0],alignments])

        print(pairwiseResultMatrix[i][0])

pairwiseDistanceMatrix =  []
i = -1
j=-1

#construct distance matrix from similarity matrix
def constructDistMatrix(pairwiseResultMatrix):
    i = -1
    for element in range(0, len(pairwiseResultMatrix)):
        i += 1

        X = pairwiseResultMatrix[element][element+1][1][0][2]

        pairwiseDistanceMatrix.append([])

        for element_2 in range(0, len(pairwiseResultMatrix)):
            #i+=1


            Y = pairwiseResultMatrix[element][element_2+1][1][0][2]
            X_2 = pairwiseResultMatrix[element_2][element_2+1][1][0][2]
            X_3 = (X + X_2)/2

            distance = round(1 - round(Y/X_3,4),4)

            pairwiseDistanceMatrix[i].append(distance)


#insert vector angles
def insertVectorAngles(encodedAaNames):
    i = -1
    j = -1
    for element in encodedAaNames:
        i+=1
        j=-1
        #encodedAaNames.append([element[0]])
        for point in range(0,len(element)):

            j+=1

            if point == len(element)-1:
                break

            encodedAaNames[i][point+1].append(round(angles[i][point+1][0],2))
            encodedAaNames[i][point + 1].append(round(angles[i][point + 2][0],2))

i=-1
j=-1

#insert plane angles
def insertPlaneAngles(encodedAaNames):
    i = -1
    j = -1
    for element in encodedAaNames:
        i+=1
        j=-1
        for point in range(0,len(element)):

            j+=1

            if point == len(element)-1:
                break

            encodedAaNames[i][point+1].append(round(plane_angles[i][point+1][0],2))

orderedCluster = []
#arrange clusters to write seperatly to file
def arrangeClusters(pairwiseDistanceMatrix):
    y = spatial.distance.squareform(pairwiseDistanceMatrix)
    # print(y)
    Z = average(y)

    clusters = fcluster(Z, 0.83, criterion='distance')
    numberOfClusters = max(clusters)


    for point in range(0, (numberOfClusters)):
        orderedCluster.append([])

    for point in range(0, len(clusters)):
        orderedCluster[clusters[point] - 1].append(point)

    for element in orderedCluster:
        print(element)

    # Z = linkage(y, 'average')
    # print(Z)
    #fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z)
    plt.show()
    #print(Z)
    #print(clusters)


#write to file
def writeToFile(encodedAaNames,orderedCluster):
    mydir="C:\projeler\pdb files\pdb file\/folder\/folder/results"
    filelist = glob.glob(os.path.join(mydir, "*.txt"))
    for f in filelist:
        os.remove(f)

    i=0
    for element_2 in orderedCluster:
        f_2 = open("results/text"+str(i)+".txt", "w")
        i+=1
        text = ""
        for element_3 in element_2:
            #for element in encodedAaNames:
                for point in encodedAaNames[element_3]:
                    text = ""
                    for point_2 in point:
                        text += str(point_2) + "  "
                    #print(text)
                    if(text.__contains__("p  d  b")):
                        text = ""
                    f_2.write("%s \r \n" % text)

        f_2.close()

readFromFile(data)
print("reading data done")
deletePoorData(data)
calculateAngles(data)
print("calculating data done")

insertEncodedNames(data)
print("inserting encoded names done")

pairwiseAlgnCode(data)
pairwiseAlignment(pairwiseAaNames)
print("pairwise alignment done")

constructDistMatrix(pairwiseResultMatrix)

arrangeClusters(pairwiseDistanceMatrix)
insertVectorAngles(encodedAaNames)
insertPlaneAngles(encodedAaNames)
print("inserting angles done")

writeToFile(encodedAaNames,orderedCluster)


print("TIME: ")
print(time.time() - start_time)