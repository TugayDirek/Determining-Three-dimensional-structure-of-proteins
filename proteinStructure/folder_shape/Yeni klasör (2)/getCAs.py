import numpy as np
import os,glob


data = []
i=0
counter = 0
j=-1
#counter_3 = False
filenames =""

def readFromFile(data):
    i = 0
    counter = 0
    j = -1
    filenames = ""
    for filename in glob.glob(os.path.join('*.pdb')):
        j = j + 1
        data.append([filename])

        i = 0
        #counter = 0
        #counter_3 = False
        with open(filename,'r') as f:
            counter += 1
            filenames += filename.replace('.pdb', ', ')
            #print(filenames)
            #print(counter)
            for line in f:

                a = i
                #counter = 0

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



input_Nn = []
output_Nn = []

def deriveNeuralNetworkCodes(input_Nn,output_Nn):
    i = -1
    j = -1
    path = "results"
    for filename in glob.glob(os.path.join('text1.txt')):
        i = i + 1
        j=-1

        input_Nn.append([])
        output_Nn.append([])


        with open(filename, 'r') as f:

            if str(f).__contains__("text1"):
                #print(f)
                for line in f:

                    if len(line.strip())==0:
                        continue

                    input_Nn[i].append([])
                    output_Nn[i].append([])
                    j = j + 1


                    line=line.rstrip() #delete last blank spaces

                    input_1 = line[0:40].strip(" ")
                    input_2 = line[41:85].strip(" ")
                    input_3 = line[86:130].strip(" ")
                    input_4 = line[131:175].strip(" ")

                    output_1 = line[180:186].strip(" ")
                    output_2 = line[187:198].strip(" ")
                    output_3 = line[198:210].strip(" ")

                    #input_final = input_1+input_2+input_3+input_4
                    #output_final = output_1+output_2+output_3

                    for element in input_1.rsplit():
                        input_Nn[i][j].append(element)

                    for element in input_2.rsplit():
                        input_Nn[i][j].append(element)

                    for element in input_3.rsplit():
                        input_Nn[i][j].append(element)

                    for element in input_4.rsplit():
                        input_Nn[i][j].append(element)

                    #input_Nn[i][j].append(input_2.replace(" ", ""))

                    output_Nn[i][j].append(output_1)
                    output_Nn[i][j].append(output_2)
                    output_Nn[i][j].append(output_3)




#readFromFile(data)
#print(data)

deriveNeuralNetworkCodes(input_Nn,output_Nn)
#print(input_Nn)
#print(output_Nn)

seq = "["

for point in range(0, len(input_Nn[0])):
    code = "["
    for point2 in range(0,len(input_Nn[0][point])):
        code_2 = input_Nn[0][point][point2]
        code = code + code_2 + ","

    code = code[:-1]
    code = code +"]"
    seq = seq + code + ",\n"
    print(code)#("("+data[0][point][3]+", "+data[0][point][4]+", "+data[0][point][5]+"),")

seq = seq[:-1]
seq = seq[:-1]
seq = seq + "]"
print(seq)

#for point in range(1, len(data[0])):
#    print("("+data[0][point][3]+", "+data[0][point][4]+", "+data[0][point][5]+"),")