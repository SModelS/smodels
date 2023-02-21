#!/usr/bin/env python

import csv

# def clean_row (row):
#     i = 0
#     while i < len(row):
#         if row[i] in ['',';','}']:
#             del row[i]
#         else :
#             i+=1
#     return row

def replace (element):
    element = element.replace('\\','')
    element = element.replace('WZi{','WZ_')
    element = element.replace('i{','')
    element = element.replace('{','')
    element = element.replace('}','')
    element = element.replace('z','_0')
    element = element.replace('n','_n')
    element = element.replace('j','J')
    element = element.replace('tcc','')
    element = element.replace('$','')
    element = element.replace('~','')
    element = element.replace('SF','_SF')
    element = element.replace('DF','_DF')
    return element

######## tables = {topology_1 : {SR_1 : {observed : x, expected : y, bgError : z}, SR_2 : {...}, ... }, topology_2 : {...}, ... } ########
tables={}
with open("ANA-SUSY-2019-09-PAPER.csv",newline='') as csvfile :
    reader = csv.reader(csvfile, delimiter=" ") #Every space is a seperator between two elements of the row
    for row in reader : #We go rw by row
        if row != [] : #If the row is not empty

            #The following can be used to retrieve the name in the label of the tex table and stored in the variable topo
            # if row[0][0:7] == '\\label{' :
            #     if row[0][11:18] == 'results':
            #         topo = row[0][19:-1]
            #     else :
            #         topo = row[0][13:-1]


            if row[0] in ['{Regions}','{Regions', 'Region'] : #If the row contains the SR names
                list_SR = [] #We store the SR order in order to write the following yields in the corresponding dictionary
                for element in row :
                    if element != '' and '\SR' in element : #If the element of the row really is a SR name
                        element = replace(element) #We reshape the SR name to match the name in the corresponding .csv file title
                        topo = element[0:element.index('_')] #The name 'topo' is used to easily loop over a given topology (not mandatory). We only keep the topology part of the name
                        list_SR.append(element) #Store the SR in the list of SR names, at the following of the others
                        if topo not in tables : #If it is the first time we encounter this topology
                            tables[topo] = {element : {'observed' : None, 'expected' : None, 'bgError' : None}}
                        else : #If the dictionary for this topoogy already exists
                            tables[topo][element] = {'observed' : None, 'expected' : None, 'bgError' : None}

            if row[0] == 'Observed' : #If the row corresponds to the 'observed' yields
                index_SR = 0 #We want to write the yields in the order given by the previous 'Region' row
                # row = clean_row(row)
                for element in row :
                    if '$' in element : #If the element is a number
                        tables[topo][list_SR[index_SR]]['observed'] = float(replace(element)) #Add the number to the dictionary following the order given by list_SR (the order of the SR obtained according to the previous 'Region' row)
                        index_SR += 1  #If we just wrote a yield, the next one will be for the next SR according to the order of the previous 'Region' row

            if row[0] == 'Fitted': #If the row corresponds to the 'expected' and 'bgError' yields
                index_SR = 0 #We want to write the yields in the order given by the previous 'Region' row
                bgError = False #In the row, the 'expected' yield always comes before the corresponding 'bgError'. The boolean is to know if the element of the 'for' loop is the 'expected' or not
                # row = clean_row(row)
                for element in row :
                    if '$' in element : #If the element is a number
                        if bgError : #If the 'expected' yield has already been written, the next number is the corresponding 'bgError'
                            tables[topo][list_SR[index_SR]]['bgError'] = float(replace(element)) #Same as for the 'Observed' row
                            bgError = False #The boolean goes to False because the next number in the row must be an 'expected' yield
                            index_SR += 1 #The next number will not be for the same SR
                        else : #If it is the first number or the previous one was a 'bgError', this one must now be an 'expected' yield
                            tables[topo][list_SR[index_SR]]['expected'] = float(replace(element)) #Same as for the 'Observed' row
                            bgError = True #The next numer will be a 'bgError' but for the same SR


for topo in tables : #To go through the results
    for SR in tables[topo] :
        print ('{' + topo + ' : {' + SR + ' : ' + str(tables[topo][SR]) + '}}') #Print the content returned if you want to check
        for event in tables[topo][SR] :
            if tables[topo][SR][event] == None : #Check that every entry of the dictionary has been rewritten
                print ("ERROR : Not all None rewritten")
