# This script reads files given as 'input' format in the HEP data website.
# e.g. see ATLAS susy analyses.

"""
This function reads the X,Y and values from the input file.
It returns the three corresponding lists of objects.

First it creates a general list, for which the entries are the lines that are read, in form of a list:
i.e. ListOfLines = [   [x1,y1,z1] , [x2,y2,z2] , ... , [xn,yn,zn]   ]
Then it extracts the entries number 0,1 and 2 for each list, and fills the vectors XArrays, YArrays ... ZArrays, that are then returned.
Note that you have to adapt the numbers of verctor that are returned and read according to the numbers of columns present in the file.
num_col is the total number of entrie you need (x,y,z for efficiencies - x,y for exlusion lines)
column is the colomun that you need
"""
def Reading_Values(input,num_col,column):
    print 'Reading the values from the input file: ',input,' .The column containing the chosen values is number: ',column , ' . \n'
    ListOfLines = []
    ListOfSRLines = []
    inputFile = open(input,'r')
    midpoint = False
    for line in inputFile:
        if (line[0] == ' ' or line[0] == '\n') and not midpoint: midpoint = True
        if line[0] != '#' and line[0] != '*' and line[0] != '\n' and line[0] != '$' and line[0] != ' ' and line[0] != '\t' and line[0] != 'M':
            #            print line
            lineElements = line.split(',') #;
            #print lineElements;
            elementsList = []
            for element in lineElements:
                #print element
                if element and element != '\n':
                    #fElement=float(element)
                   element = element.replace(' ','')
                   element = element.replace('\n','')
                   element = element.replace('$\oplus$','+')
                   #print element
                   elementsList.append(element)
            if midpoint: ListOfSRLines.append(elementsList)
            else: ListOfLines.append(elementsList)
    inputFile.close()
    
    XArray = []
    YArray = []
    ZArray = []
    #print ListOfLines
    #print ListOfSRLines

    # If saves in the third list the values contained in the column number you specified in the parameter 'column'.
    #if(num_col ==3):
    for list in ListOfLines:
       XArray.append(list[0])
       YArray.append(list[1])
       ZArray.append(list[column])
    print ListOfSRLines
    return XArray, YArray, ZArray, ListOfSRLines

    #elif(num_col ==2):
    #    if(column == num_col):
    #        for list in ListOfLines:
    #            XArray.append(list[0])
    #            YArray.append(list[1])
    #        return XArray, YArray



"""
This function produces the efficiency maps: it multiplies the values for acceptance and efficiency and creates the .txt files for each region
The input parameters are the two name of the Acc and Eff files;
topo and SR are used to create the name of the output files.
BE CAREFUL if you want to divide or not by 10.000   ( i.e. if the values given are in percentage or absolute ): you can state this option
by inputting a normalization value in Norm
"""
def Map_Multiplier(topo, SR, accFile, effFile,num_col,column, Norm):
    X1,Y1,Acc,SRLines = Reading_Values(accFile,num_col,column)
    X2,Y2,Eff,unused  = Reading_Values(effFile,num_col,column)

    SRName = []
    #print SRLines
    for line in SRLines:
        sr = line[2]
        print sr
        found = False
        for s in SRName:
            if s == sr: found = True
        if not found: SRName.append(sr)

    #print SRName
    SRList = []
    for name in SRName: SRList.append([])

    for x1,y1,acc in zip(X1,Y1,Acc):
        for x2,y2,eff in zip(X2,Y2,Eff):
            if x1==x2 and y1==y2:
                for i in range(len(SRLines)):
                    if SRLines[i][0] == x1 and SRLines[i][1] == y1:
                        for j in range(len(SRName)):
                            if SRName[j] == SRLines[i][2]:     
                                SRList[j].append(x1 + ' ' + y1 + ' ' + str(float(acc)*float(eff)/Norm) + '\n')

    #print SRList[0]
    for i in range(len(SRName)):
        name = 'EffMap_'+topo+"_"+SRName[i]+".txt"
        outputMap = open(name,'w')
        outputMap.write('# MassX , MassY , Eff*Acc  '+'\n')
        for line in SRList[i]:
            outputMap.write(line)
        outputMap.close()
        print name + ' written!'

"""
    outputMap = open('EffMap_'+topo+"_"+SR+".txt",'w')
    outputMap.write('# MassX , MassY , Eff*Acc  '+'\n')
    for x1,y1,acc in zip(X1,Y1,Acc):
        for x2,y2,eff in zip(X2,Y2,Eff):
            if x1==x2 and y1==y2:
#   print x1 + ' ' + x2 + ' ' + y1 + ' ' + y2 + ' \n'  # just to check if the selected values from the two files matches
                outputMap.write(x1 + ' ' + y1 + ' ' +  str(float(acc)*float(eff)/Norm) + '\n')
    print "Map ",'EffMap_'+topo+"_"+SR+".txt", ' written!'
"""


"""
This function simply rewrite in a file .dat that you want to plot, in a SModelS friendly format. It takes the values of the arrays from the Reading_Values function.
Give as parameters the two array you want to plot, and the name of the output file.
With 'type of data' you specify what kind of values are you extracting. 
"""
def Simple_Map_Producer(X,Y,Z,type_of_data,outputName):
    output = open(outputName+'.dat','w')
    output.write('# MassX , MassY ' + type_of_data+'\n')
    for x,y,z in zip(X,Y,Z):
        output.write(x+' '+y+' '+z +'\n')

def Simple_Exclusion_Producer(X,Y,type_of_data,outputName):
    output = open(outputName+'.dat','w')
    output.write('# MassX , MassY ' + type_of_data+'\n')
    for x,y in zip(X,Y):
        output.write(x+' '+y+'\n')

#files = [25,27,29,31,33,35]
#SR	  = ['SF-loose','SF-tight','DF-100','DF-150','DF-200','DF-300']

files = [93]#,89,93,95,97]
topo  = ['T6bbWW']#,'T2tt','T6bbWW','T6bbWW','T6bbWW']
SR	  = ['']#,'slep-b','slep-c','slep-d','slep-e']

for i in range(len(files)):
    ACC = 'Table'+str(files[i])+'.csv'
    EFF = 'Table'+str(files[i]+1)+'.csv'
    Map_Multiplier(topo[i],SR[i],ACC,EFF,2,2,1) #norm=10000


#X,Y = Reading_Values("T2cc_Obs_Excl.csv",2,2) #Obs_Line.dat
#Simple_Exclusion_Producer(X,Y,"Obs_Excl","T2cc_Obs_Excl.dat")
#X,Y = Reading_Values("T2cc_Exp_Excl.csv",2,2) #Exp_Line.dat
#Simple_Exclusion_Producer(X,Y,"Exp_Excl","T2cc_Exp_Excl.dat")

#X,Y,Z = Reading_Values("T2cc_Obs_UL.csv",3,3)
#Simple_Map_Producer(X,Y,Z,"Obs_UL","T2cc_Obs_UL.dat")









