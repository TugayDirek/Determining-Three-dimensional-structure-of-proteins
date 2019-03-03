import numpy as np
import os,glob
import vpython
import math
x = []
i=0
counter = 0
j=-1
#counter_3 = False

for filename in glob.glob(os.path.join('*.pdb')):
    j = j + 1
    x.append([filename])
    i = 0
    counter = 0
    #counter_3 = False
    with open(filename,'r') as f:
        for line in f:
            #if (counter_3):
                #print("first")
                #print(counter_3)
                #break

            #print(line)
            a = i
            counter = 0

            if line.__contains__("ENDMDL"):
                break

            if(line.__contains__("ATOM") and line.__contains__("CA")):
               x[j].append(["CA"])
               i = i + 1

            for word in line.split():
               if line.__contains__("ATOM") and line.__contains__("CA"):

                #print(x[0][0])
                counter = counter + 1
                if (counter == 2):
                 x[j][a + 1].append(word)
                elif(counter == 4):
                 x[j][a+1].append(word)
                elif (counter == 7):
                 x[j][a+1].append(word)
                elif (counter == 8):
                 x[j][a+1].append(word)
                elif (counter == 9):
                 x[j][a+1].append(word)


               #x[i].append(word)
               #print(word)


#print(x)
print("Alpha Carbons")
for element in x:
    for point in element:
        print(point)

''''                
                 for counter_2 in range(1,len(x[j])):

                     if (x[j][counter_2].__contains__(word)):
                      #print(x[j][counter_2])
                      counter_3 =  True
                      break
                 if(counter_3):
                     break
'''

angles = []
vectors = []
#angles_2 = []
i=-1
for element in x:
    i+=1
    angles.append([element[0]])
    for point in range(0,len(element)):

        if point == len(element)-3:
            break

        a = np.array([float(element[point+1][3]), float(element[point+1][4]), float(element[point+1][5])])
        b = np.array([float(element[point+2][3]), float(element[point+2][4]), float(element[point+2][5])])
        c = np.array([float(element[point+3][3]), float(element[point+3][4]), float(element[point+3][5])])


        ba = a - b
        bc = c - b

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        #print(cosine_angle)
        angle = np.arccos(cosine_angle)

        #ann = vpython.vector(-2.05 , -1.247 , 2.976).diff_angle(vpython.vector(-0.191 , 3.68 ,  1.014))
        #print(ba)
        #print(bc)
        #angles_2.append(vpython.degrees(ann))


        #print(np.degrees(angle))
        angles[i].append([np.degrees(angle)])

i=-1
for element in x:
    i += 1
    vectors.append([element[0]])
    for point in range(0,len(element)):
        
        if point == len(element)-2:
            break

        a = np.array([float(element[point + 1][3]), float(element[point + 1][4]), float(element[point + 1][5])])
        b = np.array([float(element[point + 2][3]), float(element[point + 2][4]), float(element[point + 2][5])])

        a_2 = vpython.vector(b[0]-a[0], b[1]-a[1], b[2]-a[2])

        vectors[i].append(a_2)
       
print("Angles")
for element in angles:
    for point in element:
        print(point)
print(angles)
#print(angles_2)
#print(len(x[0])+len(x[1])+len(x[2]))
#print(len(vectors))
#v1 = vpython.vector( 3.042 , 0.594, -2.146)
#gg = vpython.rotate(v1,angle=vpython.radians(angles[0]),axis=vpython.vector(1,1,1))
#print(v1)
#print(gg)
#print(vectors)
'''''
print("Vectors")
for element in vectors:
    for point in element:
        print(point)
'''
vectors_mag = []


i=-1
for element in vectors:
    i += 1
    vectors_mag.append([element[0]])
    for point in range(0, len(element)):

        if point == len(element)-1:
            break

        vectors_mag[i].append(vpython.mag(element[point + 1]))

'''''
print("Vector Magnitudes")
for element in vectors_mag:
    for point in element:
        print(point)
'''

i=-1
trans_vectors = []

for element in x:
    i += 1
    trans_vectors.append([element[0]])
    for point in range(0, len(element)):


        if point == len(element)-1:
            break

        a = np.array([float(element[point+1][3]), float(element[point+1][4]), float(element[point+1][5])])
        b = np.array([float(element[1][3]), float(element[1][4]), float(element[1][5])])



        trans_vectors[i].append(a - b)
'''''
print("Transformed Points")
for element in trans_vectors:
    for point in element:
        print(point)
'''

aadict = {
"GLU":"1 0 0 0 0",
"LEU":"0 1 0 0 0",
"THR":"0 0 1 0 0",
"PRO":"0 0 0 1 0",
"ASP":"0 0 0 0 1",
"GLN":"1 1 0 0 0",
"HIS":"1 0 1 0 0",
"PHE":"1 0 0 1 0",
"ILE":"1 0 0 0 1",
"MET":"1 1 0 1 0",
"SER":"1 1 0 0 1",
"TYR":"1 1 1 0 0",
"ASN":"1 1 1 1 0",
"LYS":"1 1 1 0 1",
"ARG":"1 1 1 1 1",
"ALA":"0 1 1 0 0",
"VAL":"0 1 1 1 0",
"GLY":"0 1 1 0 1",
"TRP":"0 1 1 1 1",
"CYS":"0 0 1 1 0"
}

encodedAaNames = []

i=-1
for element in x:
    i+=1
    encodedAaNames.append([element[0]])
    for point in range(0,len(element)):

        if point == len(element)-3:
            break

        a = ((element[point+1][2]))
        b = ((element[point+2][2]))
        c = ((element[point+3][2]))

        a_distance = 0
        b_distance = 0
        c_distance = 0

        a_4 = [0,0,0]#vpython.vector(0,0,0)
        b_4 = [0,0,0]#vpython.vector(0,0,0)
        c_4 = [0,0,0]#vpython.vector(0,0,0)

        a_vectors = []
        b_vectors = []
        c_vectors = []

        a_2 = np.array([float(element[point + 1][3]), float(element[point + 1][4]), float(element[point + 1][5])])
        b_2 = np.array([float(element[point + 2][3]), float(element[point + 2][4]), float(element[point + 2][5])])
        c_2 = np.array([float(element[point + 3][3]), float(element[point + 3][4]), float(element[point + 3][5])])

        j = -1
        for element_2 in x:
            j += 1
            #encodedAaNames.append([element[0]])
            for point_2 in range(0, len(element_2)):

                if point_2 == len(element_2) - 1:
                    break

                a_3 = np.array([float(element_2[point_2 + 1][3]), float(element_2[point_2 + 1][4]), float(element_2[point_2 + 1][5])])


                a_distance = math.sqrt(math.pow(a_2[0]- a_3[0],2) + math.pow(a_2[1]- a_3[1],2) +math.pow(a_2[1]- a_3[1],2))
                b_distance = math.sqrt(math.pow(b_2[0] - a_3[0], 2) + math.pow(b_2[1] - a_3[1], 2) + math.pow(b_2[1] - a_3[1], 2))
                c_distance = math.sqrt(math.pow(c_2[0] - a_3[0], 2) + math.pow(c_2[1] - a_3[1], 2) + math.pow(c_2[1] - a_3[1], 2))

                #print(c_distance)

                if(a_distance<10):
                    #a_4 = a_4 + #vpython.vector(a_2[0]-a_3[0],a_2[1]-a_3[1],a_2[2]-a_3[2])
                    a_4[0] += a_2[0]-a_3[0]
                    a_4[1] += a_2[1] - a_3[1]
                    a_4[2] += a_2[2] - a_3[2]

                if (b_distance < 10):
                    #b_4 = b_4 + vpython.vector(b_2[0]-a_3[0],b_2[1]- a_3[1],b_2[2]- a_3[2])
                    b_4[0] += b_2[0]-a_3[0]
                    b_4[1] += b_2[1] - a_3[1]
                    b_4[2] += b_2[2] - a_3[2]

                if (c_distance < 10):
                    #c_4 = c_4 + vpython.vector(c_2[0]-a_3[0],c_2[1]- a_3[1],c_2[2]- a_3[2])
                    c_4[0] += c_2[0]-a_3[0]
                    c_4[1] += c_2[1] - a_3[1]
                    c_4[2] += c_2[2] - a_3[2]

                #print(a_3[0])

        #encodedAaNames[i].append([aadict.get(a),aadict.get(b),aadict.get(c),a_4.mag,b_4.mag,c_4.mag])
        encodedAaNames[i].append([aadict.get(a),aadict.get(b),aadict.get(c),a_4[0],a_4[1],a_4[2],b_4[0],b_4[1],b_4[2],c_4[0],c_4[1],c_4[2]])


print("encodedNames")
print(encodedAaNames)


i=-1
j=-1
for element in encodedAaNames:
    i+=1
    j=-1
    #encodedAaNames.append([element[0]])
    for point in range(0,len(element)):

        j+=1

        if point == len(element)-1:
            break

        encodedAaNames[i][point+1].append(angles[i][point+1][0])


print(encodedAaNames)

f_2 = open("text.txt","w")
text = ""
print("encodedAaNames Points")
for element in encodedAaNames:
    for point in element:
        text = ""
        for point_2 in point:
            text += str(point_2) + "  "
        print(text)
        if(text.__contains__("p  d  b")):
            text = ""
        f_2.write("%s\r\n" % text)

f_2.close()