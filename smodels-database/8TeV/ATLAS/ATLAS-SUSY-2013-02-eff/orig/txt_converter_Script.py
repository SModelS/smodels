import os, sys

thisDirectory = os.getcwd()
mapsInFolder = os.listdir(thisDirectory)

mapsList = []
for file in mapsInFolder:
    if '.txt' in file and '~' not in file:   # skipping the weird ~ objects
        mapsList.append(file)
mapsList.sort()                          # Sort Alphabetically

#print 'List of maps:',mapsList



regionNameList = []



for Map in mapsList :
    motherMass= []
    daughterMass = []
    efficiency =[]
    #mappa = open(Map,"r").read().split('/n')
    mappa = open(Map,"r")
    mappaOut = open(Map[:-3]+'dat',"w")
    finalList =[]
    
    for line in mappa:
        
            lineElement = line.split('\t')
            for element in lineElement:
                #print element
                if '\n' in element:
                    element = element.replace('\n','')
            
            convertedLineElement = []
            for element in lineElement:
                if element and element !='\n':
                    if 'x' not in element and 'y' not in element:
                        convertedLineElement.append(element)
            finalList.append(convertedLineElement)
#            print finalList



    for line in finalList[2:]:
        if (len(line)>=6):
            
            motherMass.append(line[0])
            daughterMass.append(line[3])
            efficiency.append(line[6])

    for i in range(0,len(motherMass)):
        mappaOut.write(motherMass[i] +'  '+ daughterMass[i] + '  ' +efficiency[i])
    mappaOut.close()


'''

helper_to_Convert = open('helper_to_Convert.txt','w')

listOfTopo = ['T1tttt']

for dataFile in mapsList:
    dataset = dataFile.replace(".txt","").replace("HEPdata.","").replace(".","_")
    for topo in listOfTopo:
        helper_to_Convert.write( topo + '.efficiencyMap.setSource( ' + dataFile + ', "txt" , objectName =  "' + dataFile + '", index = None , dataset = "' + dataset + '" ) ' + '\n')
        helper_to_Convert.write( topo + '.efficiencyMap.setStatistics( observedN=  ,  expectedBG=   , bgError=   )'  + '\n' + '\n'+ '\n')

helper_to_Convert.close()


'''



'''

sourceFile = ''
sourceType = ''
obName = ''
dSet = ''

for topo in listOfTopo:



 T5WZ.efficiencyMap.setSource( "orig/SUS13012_XsecLimits_T5VV.root", "root", objectName = "h_EffAcc_3NJet6_1250HT1500_450MHTinf", index = None, dataset="3NJet6_1250HT1500_450MHTinf" )
 T5WZ.efficiencyMap.setStatistics ( observedN=23, expectedBG=17.6, bgError=4.1 )



#    region.replace("/n", " ")
#    print region



'''












