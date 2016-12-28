#!/usr/bin/env python

"""
.. module:: tools.BrowserGUI.dataScreensGUI
   :synopsis:  Contains all the definitions for the main screens/layouts for the browser GUI
                        (Kivy-based)

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import sys,os
from kivy.logger import Logger
import webbrowser  #To open URL links
from kivy.uix.image import Image
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.core.window import Window
from kivy.uix.listview import ListItemButton, ListView
from kivy.adapters.listadapter import ListAdapter
from kivy.adapters.simplelistadapter import SimpleListAdapter
from kivy.graphics import Color,Rectangle
from kivy.metrics import sp
Window.clearcolor = [0.08627450980392157, 0.10196078431372549, 0.10588235294117647, 1.0]
from dataObjectsGUI import ValItem, ExpResItem, PropItem, datasetItem, txnameItem
sys.path.append(os.path.join(os.getenv("HOME"),'smodels'))
from smodels.tools.physicsUnits import GeV


class MainScreen(BoxLayout):
    """
    The Main screen. Holds all the sub-screens as children.
    
    :ivar browser: databaseBrowser.Browser object holding the experimental results
    """
    
    def __init__(self,browser, **kwargs):
        super(MainScreen, self).__init__(**kwargs)
        
        #Set up the main screen options
        self.orientation = 'horizontal'
        self.padding = [sp(10),sp(10),sp(10),sp(10)]
        self.spacing = sp(15)
        self.browser = browser
        self.background_color = (0.21, 0.7, 1., 1.)
        self.size_hint_x = None
        
        #Builds the four sub-screens:        
        #The properties screen (contains all Browser atributes):
        self.propScreen = PropertiesScreen(browser)
        #The experimental results screen (contains all selected ExpResult objects):
        self.resScreen = ResultsScreen(browser) 
        #Information screen (contains the globalInfo object and the list of datasets for the selected ExpResult)
        self.infoScreen = InfoScreen(browser._selectedExpResults[0])
        #Dataset information screen (contains the dataInfo object and the list of Txnames for the selected DataSet) 
        self.dataInfoScreen = DataInfoScreen(browser._selectedExpResults[0].datasets[0])
        
        #Set up the relative size of the sub-screens:        
        self.propScreen.size_hint = (0.22,1.)
        self.resScreen.size_hint = (0.22,1.)
        self.infoScreen.size_hint  = (0.3,1.)
        self.dataInfoScreen.size_hint  = (0.25,1.)
        
        #Add the screens to the main screen:
        self.add_widget(self.propScreen)
        self.add_widget(self.resScreen)
        self.add_widget(self.infoScreen)
        self.add_widget(self.dataInfoScreen)
        
        #Set up the main screen absolute size:
        self.width = sp(1200)  #Use sp to build screen-independent pixels
        Window.size = (self.width,list(Window.size)[1])
        
        
        

class PropertiesScreen(BoxLayout):
    """
    The Properties screen. Contains a selectable list of all
    the attributes of the Browser.
    
    :ivar browser: databaseBrowser.Browser object holding the selected experimental results
    """
    
    def __init__(self, browser,**kwargs):
        super(PropertiesScreen, self).__init__(**kwargs)

        #Set up global options                        
        self.orientation = 'vertical'
        self.spacing = sp(10)
        
        #Build list of attributes (property items)
        self.propItems = []        
        for attr in sorted(browser.getAttributes()):
            if  attr[0] == '_': continue #Skip private attributes
            propItem = PropItem(attr,browser)  #Create the item based on the attribute
            propItem.valScreen = ValuesScreen(attr,browser) #Each item contains the values screen
            self.propItems.append(propItem) 
         
        #Create Button list of the properties:       
        args_converter = lambda row_index, obj: {'text': obj.name,
                                         'size_hint_y': None,
                                         'height': sp(27),
                                          'deselected_color' : [0.35, 0.35, 0.35, 1.0],
                                          'selected_color' : [0.6,0.6,0.6,1.],                                          
                                          'font_size' : sp(20),                                        
                                         'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.propItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='multiple')
        
        #Label for the screen
        l = Label(text="Properties",font_size=sp(30))
        l.size_hint_y = 0.1  #Relative label size        
        #List of buttons:     
        propList = ListView(adapter=list_adapter)
        propList.children[0].bar_width = 5  #Increases the size of the scroll bar
        propList.size_hint_y = 0.75 #Relative list size
        #Button to select the results (apply the filter)
        filterButton = Button(text="Select", font_size =25)
        filterButton.size_hint_y = 0.08 #Button relative size
        filterButton.bind(on_press=selectResults) #When pressed, call selectResults
        #Button to clear the selection:
        clearButton = Button(text="Clear Selection", font_size =25)
        clearButton.size_hint_y = 0.06
        clearButton.bind(on_press=clearSelection)#When pressed, call clearSelection
        
        #Add widgets:
        self.add_widget(l)
        self.add_widget(propList)
        self.add_widget(filterButton)
        self.add_widget(clearButton)


class ResultsScreen(BoxLayout):
    """
    The Selected Results screen. Contains a selectable list of all
    the experimental results, after selection.
    
    :ivar browser: databaseBrowser.Browser object holding the selected experimental results    
    """
    
    def __init__(self,browser, **kwargs):
        super(ResultsScreen, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = sp(12)        
        self.browser = browser
        self.resItems = []  #List of experimental result items
        
    def clear(self):
        """
        Remove all widgets 
        """
        self.clear_widgets()
        
    def update(self):
        """
        Update the list of selected results
        """
        self.clear()
        self.resItems = []
        
        #Re-build the list of selected results using the active browser:
        for expRes in self.browser:
            self.resItems.append(ExpResItem(expRes,self.browser))
        args_converter = lambda row_index, obj: {'text': obj.name,
                                         'size_hint': (1,None),
                                         'height': sp(40),
                                          'deselected_color' : [96./255.,384./255.,503./255.,1.],
                                          'selected_color' : [147./255.,473./255.,600./255.,1.],
                                          'font_size' : sp(20),
                                         'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.resItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='single',propagate_selection_to_data=True)
        resList = ListView(adapter=list_adapter)
        resList.children[0].bar_width = sp(5)
        resList.size_hint_y = 0.9
        
        l = Label(text="Selected:",font_size=sp(30))
        l.size_hint_y = 0.1
         
        #Add widgets                    
        self.add_widget(l)
        self.add_widget(resList)
            

class InfoScreen(BoxLayout):
    """
    The Info screen. Contains the information about a given experimental result.
    Includes both the meta data in expRes.globalInfo and the list of datasets in expRes.datasets.
    
    :ivar expRes: smodels.experiment.databaseObjects.ExpResult object    
    """
    
    def __init__(self,expRes, **kwargs):
        super(InfoScreen, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = sp(12)
        self.dataItems = []
            
    def clear(self):
        """
        Clear the screen (remove all widgets)
        """
        self.clear_widgets()
        
    def updateWith(self,expRes):
        """
        Update the screen with information about expRes.
        :param expRes: smodels.experiment.databaseObjects.ExpResult object
        """                
        self.clear()
        self.dataItems = []
        
        #First block is global Info
        info_items = []
        for name,val in expRes.globalInfo.__dict__.items():            
            if 'http' in str(val):
                iline = '[b] [size=17]'+str(name) + '[/size] [/b] : [ref='+str(val)+'] [color=ff9999] '
                iline  += str(val) + '[/color] [/ref]'
            else:
                iline = '[b] [size=17]'+str(name) + '[/size] [/b] : ' + str(val)
            info_items.append(iline)
        list_adapter = SimpleListAdapter(data=sorted(info_items),
                                     cls=myLabel, selection_mode='single')        
        globalInfoList = ListView(adapter=list_adapter)
        globalInfoList.children[0].bar_width = sp(5)
        
        
        #Second block is dataset list:
        for dataset in expRes.datasets:
            self.dataItems.append(datasetItem(dataset))        
        args_converter = lambda row_index, obj: {'text': str(obj.name),
                                 'size_hint': (0.9,None),
                                 'height': sp(40),
                                 'deselected_color' : [96./255.,384./255.,503./255.,1.],
                                 'selected_color' : [147./255.,473./255.,600./255.,1.],
                                  'font_size' : sp(20),
                                   'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.dataItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='single',propagate_selection_to_data=True)
        datasetList = ListView(adapter=list_adapter)
        datasetList.children[0].bar_width = sp(5)
                
        #Define relative sizes
        globalInfoList.size_hint = (1.,0.5)
        datasetList.size_hint = (1.,0.4)
        
        infoLabel = Label(text="Experimental Result Info:",font_size=30,size_hint_y=0.07)
        infoLabel.bind(size=infoLabel.setter('text_size')) 
        self.add_widget(infoLabel)
        self.add_widget(globalInfoList)
        dataLabel = Label(text="Data Sets:",font_size=30,size_hint_y=0.08, halign='left')
        dataLabel.bind(size=dataLabel.setter('text_size')) 
        self.add_widget(dataLabel)
        self.add_widget(datasetList)


class DataInfoScreen(BoxLayout):
    """
    The DataSet Info screen. Contains the information about a given Data Set.
    Includes both the meta data in dataset.dataInfo and the list of TxNames in dataset.txnameList.
    
    :ivar dataset: smodels.experiment.datasetObject.DataSet object    
    """
    
    def __init__(self,dataset, **kwargs):
        super(DataInfoScreen, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = sp(12)
            
    def clear(self):
        """
        Clear the screen (remove all widgets)
        """        
        self.clear_widgets()
        
    def updateWith(self,dataset):
        """
        Update the screen with information about dataset.
        :param dataset: smodels.experiment.datasetObject.DataSet object
        """              
        self.clear()
        self.txnameItems = []
        for txname in dataset.txnameList:
            txItem = txnameItem(txname)
            txItem.txnameScreen = TxnameScreen(txname)
            self.txnameItems.append(txItem)
        
                     
        #First block is dataset Info
        info_items = []
        for name,val in dataset.dataInfo.__dict__.items():
            info_items.append('[b] [size=17]'+str(name) + '[/size] [/b] : ' + str(val))
        list_adapter = SimpleListAdapter(data=sorted(info_items),
                                     cls=myLabel, selection_mode='single')        
        dataInfoList = ListView(adapter=list_adapter)
        dataInfoList.children[0].bar_width = sp(5)
        
        
        #Second block is txnames list:
        args_converter = lambda row_index, obj: {'text': str(obj.name),
                                 'size_hint': (1,None),
                                 'height': sp(40),
                                 'deselected_color' : [96./255.,384./255.,503./255.,1.],
                                 'selected_color' : [147./255.,473./255.,600./255.,1.],
                                  'font_size' : sp(20),
                                   'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.txnameItems,args_converter=args_converter,                                   
                                    cls=ListItemButton, selection_mode='single',propagate_selection_to_data=True)
        txnameList = ListView(adapter=list_adapter)
        txnameList.children[0].bar_width = sp(5)
        
        #Third block is Upper Limit getter
        self.ulGetter = ULgetter(dataset)
        
        dataInfoList.size_hint = (1.,0.25)
        txnameList.size_hint = (None,0.3)
        txnameList.width = sp(max([len(str(txitem.name)) for txitem in self.txnameItems])*15.)
        txnameList.width = sp(min(txnameList.width,self.width))
        self.ulGetter.size_hint = (1.,0.4) 
        
        infoLabel = Label(text="Data Set Info:",font_size=sp(30),size_hint_y=0.07)
        infoLabel.bind(size=infoLabel.setter('text_size')) 
        self.add_widget(infoLabel)
        self.add_widget(dataInfoList)
        txnameLabel = Label(text="TxNames:",font_size=sp(30),size_hint_y=0.08, halign='left')
        txnameLabel.bind(size=txnameLabel.setter('text_size')) 
        self.add_widget(txnameLabel)
        self.add_widget(txnameList)
        self.add_widget(self.ulGetter)

        

class ValuesScreen(BoxLayout):
    """
    The Values screen. Contains the possible values for a given browser attribute.
    It is associated to a given PropItem object.
    
    :ivar text: name of the associate attribute/property
    :ivar browser: databaseBrowser.Browser object holding the experimental results
    """    
    def __init__(self,text,browser, **kwargs):
        super(ValuesScreen, self).__init__(**kwargs)

        self.orientation = 'vertical'
        self.spacing = 12
        self.padding = 8
        self.name = text
        self.browser = browser
        self.valItems = []

        #Set up the background color (draw a retancle in the background)
        with self.canvas.before:
            Color(0.6,0.6,0.6,1.0)            
            self.rect = Rectangle(size=self.size,pos=self.pos)
            self.bind(pos=update_rect, size=update_rect)
        
    def create(self):
        """
        Create the screen, if not yet created (for improved time performance).
        Checks the browser for all possible values of the attribute/property and
        creates a selectable list with the set of all values
        """
        if not self.children:   
            #Get all possible values for the attributes using the full database.
            #To make sure that the values are taken from all exp results, temporarily restore the
            #browser to the pre-selection state
            tmpRes = self.browser._selectedExpResults[:]
            self.browser.loadAllResults()            
            values = dict([[str(v),v] for v in self.browser.getValuesFor(attribute=self.name)])
            self.browser._selectedExpResults = tmpRes                       
            for val_str in sorted(values.keys()):
                self.valItems.append(ValItem(self.name,values[val_str]))        
            
            args_converter = lambda row_index, obj: {'text': str(obj.value),
                                                'deselected_color' : [0.5, 0.5, 0.5, 1.0],
                                                'selected_color' : [0.95,0.95,0.95,1.],                                                     
                                              'size_hint_y' : None,
                                              'height' : sp(30),
                                             'on_press' : obj.changeSel}
                
            list_adapter = ListAdapter(data=self.valItems,args_converter=args_converter,
                                        cls=ListItemButton, selection_mode='multiple')            

            valList = ListView(adapter=list_adapter)
            valList.children[0].bar_width = sp(5)
            valList.size_hint_y = 0.9
            
            valLabel = Label(text='Values for %s' %self.name, font_size = sp(20))
            valLabel.size_hint_y = 0.1
            
            self.add_widget(valLabel)
            self.add_widget(valList)
            
            
class TxnameScreen(BoxLayout):
    """
    The TxName screen. Contains the meta data about a given Txname.
        
    :ivar txname: smodels.experiment.txnameObject.TxName object
    """        
    
    def __init__(self,txname, **kwargs):
        super(TxnameScreen, self).__init__(**kwargs)

        self.orientation = 'vertical'
        self.txname = txname
        self.padding = 10
        self.infoItems = []
        with self.canvas.before:
            Color(0.10390763765541741, 0.1238898756660746, 0.12788632326820604, 1.0)
            self.rect = Rectangle(size=self.size,pos=self.pos)
            self.bind(pos=update_rect, size=update_rect)
        
    def create(self):
        """
        Create the screen, if not yet created (for improved time performance).
        Loads all the meta data information (skipping unwanted fields) and
        builds a simple scrollable liste with the info. 
        """        
        
        if not self.children:            
            for name,val in self.txname.__dict__.items():
                if name.lower() == 'efficiencymap' or name.lower() == 'upperlimits': continue
                if name.lower() == 'txnamedata': continue
                if name.lower() == 'globalinfo': continue
                if name[0] == '_': continue
                if 'http' in str(val):
                    iline = '[b] [size=17]'+str(name) + '[/size] [/b] : [ref='+str(val)+'] [color=ff9999] '
                    iline  += str(val) + '[/color] [/ref]'
                else:
                    iline = '[b] [size=17]'+str(name) + '[/size] [/b] : ' + str(val)                
                self.infoItems.append(iline)
            list_adapter = SimpleListAdapter(data=sorted(self.infoItems),
                                         cls=myLabel, selection_mode='single')        
            txnameInfoList = ListView(adapter=list_adapter)
            txnameInfoList.children[0].bar_width = sp(5)
            txnameInfoList.size_hint_y = 0.6
            txnameLabel = Label(text='%s Info:' %self.txname.getInfo('txName'), font_size = sp(20))
            txnameLabel.size_hint_y = 0.1
            
            if os.path.isfile('./feyn/'+self.txname.getInfo('txName')+'_feyn.png'):
                txImage = Image(source='./feyn/'+self.txname.getInfo('txName')+'_feyn.png')
                txImage.size_hint = (0.2,0.5)
                with txImage.canvas.before:
                    Color(1.,1.,1.,1.)            
                    txImage.rect = Rectangle(size=txImage.size,pos=txImage.pos)
                    txImage.bind(pos=update_rect, size=update_rect)
            else:
                txImage = myLabel(text="[i] [size=20] (TxName image not found) [/size] [/i]")
                txImage.markup = True
                txImage.size_hint_y = 0.5
            
            self.add_widget(txnameLabel)
            self.add_widget(txnameInfoList)
            self.add_widget(txImage)
            

class ULgetter(FloatLayout):
    """
    The screen for the upper limit calculator.
    In the case of an EM result/dataset, only contains the fixed upper limit for the
    corresponding dataset (signal region). For a UL-type result/dataset, contains a dynamic
    upper limit calculator.
    
    :ivar dataset:  smodels.experiment.datasetObject.DataSet object  
    """
    
    def __init__(self,dataset,**kwargs):
        super(ULgetter, self).__init__(**kwargs)
        
        self.cols = 2        
        self.spacing = 5
        self.padding = 2
        self.dataset = dataset
        self.txname = None
        
        #Defines the background color
        with self.canvas.before:
            Color(0.10761166226227092, 0.15942468483299393, 0.16341030195381884, 1.0)
            self.rect = Rectangle(size=self.size,pos=self.pos)
            self.bind(pos=update_rect, size=update_rect)
        
        #If UL-type, build the calculator
        if self.dataset.getValuesFor('dataType')[0] == 'upperLimit':
            self.createULcalc()
        else: #If EM-type, just show the fixed value
            fixedUL = Label(text = "Upper Limit for %s:\n   %s"
                            %(dataset.getValuesFor('dataId')[0],str(dataset.getValuesFor('upperLimit')[0])),font_size=20)
            fixedUL.size_hint = (0.5,0.5)
            fixedUL.pos_hint = {'top' : 1., 'right' : 0.7}
            self.add_widget(fixedUL)
    
    def createULcalc(self):
        """
        Builds the UL calculator for UL-type results/datasets.
        Contains the mass input box and other labels.
        """
        ulLabel = Label(text="[i] (select TxName) [/i]",font_size=23,size_hint=(1,0.2), halign='left', 
                        pos_hint = {'top' : 1.,'right' : 1.})
        ulLabel.bind(size=ulLabel.setter('text_size'))
        ulLabel.shorten = True
        ulLabel.markup = True
        
        massALabel  = Label(text="Masses 1 (GeV):",font_size=15,size_hint=(0.5,0.2), 
                            pos_hint = {'top' : 0.8,'right' :0.45})        
        massBLabel  = Label(text="Masses 2 (GeV):",font_size=15,size_hint=(0.5,0.2), 
                            pos_hint = {'top' : 0.6,'right' : 0.45})
        massAInput = TextInput(text="",size_hint = (0.6,0.14), 
                               pos_hint = {'top' : 0.75,'right' : 1.})
        massAInput.multiline = False
        massAInput.foreground_color = [0.44, 0.44, 0.44, 1.0]
        massBInput = TextInput(text="",size_hint = (0.6,0.14), 
                               pos_hint = {'top' : 0.6,'right' : 1.})
        massBInput.multiline = False
        massBInput.foreground_color = [0.44, 0.44, 0.44, 1.0]
        
        ulbutton = Button(text="Get UL",size_hint = (0.2,0.2), pos_hint = {'top' : 0.35,'right' : 0.2})
        ulbutton.bind(on_press=self.getUL)
        resLabel  = Label(text="",font_size=20, size_hint = (0.2,0.7), 
                          pos_hint = {'top' : 0.6,'right' : 0.65})
        resLabel.markup = True
        
        self.ulLabel = ulLabel        
        self.button = ulbutton
        self.massInput = [massAInput,massBInput]
        self.resLabel =resLabel
        self.add_widget(ulLabel)
        self.add_widget(massALabel)
        self.add_widget(massAInput)
        self.add_widget(massBLabel)
        self.add_widget(massBInput)        
        self.add_widget(resLabel)
        self.add_widget(ulbutton)

    def updateWith(self,txname):
        """
        For UL-results: Updates the calculator with a new TxName
        (It has no effect for  EM-results)
         
        :param txname: smodels.experiment.txnameObject.TxName object
        """
        
        if not self.dataset.getValuesFor('dataType')[0] == 'upperLimit':
            return
        self.txname = txname
        self.ulLabel.text = "UL for %s :" %txname.txName
        #if not txname.txnameData._data:
        #    txname.txnameData.loadData()
        # massarray =  txname.txnameData._data[0][0]
        massarray =  "FIXME"
        for ib in range(len(self.massInput)):
            if ib == 0:
                mLabel = "m"
            else:
                mLabel = "M"
            self.massInput[ib].text = str([mLabel+str(i) for i in range(1,len(massarray[ib])+1)])
            self.massInput[ib].text = self.massInput[ib].text.replace("'","")
     
    def clear(self):
        """
        For UL-results: clear the calculator. Removes the TxName associated to it. 
        (It has no effect for  EM-results)        
        """
        if not self.dataset.getValuesFor('dataType')[0] == 'upperLimit':
            return
        
        self.ulLabel.text = "[i] (select TxName) [/i]"
        self.txname = None
        for massInput in self.massInput:
            massInput.text = ""
        self.resLabel.text = ""
        
    def getUL(self,button=None):
        """
        For UL-results: tries to compute the UL for the mass input given.
        Updates the result to the upper limit calculator screen.
        (It has no effect for  EM-results)        
        """
        if not self.dataset.getValuesFor('dataType')[0] == 'upperLimit':
            return
        
        try:
            massarray = []
            for massStr in self.massInput:
                masses = [m*GeV for m in list(eval(massStr.text))]
                massarray.append(masses)
        except:
            self.resLabel.text = 'Wrong mass input'
            return
        
        ul = self.txname.txnameData.getValueFor(massarray)
        
        if  ul is False:
            self.resLabel.text = 'Error computing limit'
        else:
            self.resLabel.text = "UL  = " + str(ul)
        
    


class myLabel(Label):
    """
    Auxiliary class to define labels aligned to the left
    """    
    def __init__(self,  *args, **kwargs):
        super(myLabel, self).__init__(*args, **kwargs)

        self.halign = 'left'
        self.markup = True
        self.bind(size=self.setter('text_size')) 
        self.shorten = True
        self.bind(on_ref_press=self.openlink)
        
    def openlink(self,mlabel,link):
        """
        Open link in the webbrowser
        :param mlabel: myLabel object
        :param link: url address
        """
        try:            
            webbrowser.open(link)
        except:
            Logger.error("Failed to open link: %s" %link)

def selectResults(button):
    """
    Use the values selected in the Properties screen to filter the experimental results
    """
    
    propScreen = button.parent
    propItems = propScreen.propItems
    restrictDict = {}
    for pitem in propItems:
        if not pitem.valScreen: continue
        for valitem in pitem.valScreen.valItems:
            if valitem.is_selected:
                if not pitem.name in restrictDict: restrictDict[pitem.name] = []
                restrictDict[pitem.name].append(valitem.value)
    
    mainScreen = propScreen.parent
    resScreen = mainScreen.resScreen
    resScreen.browser.selectExpResultsWith(**restrictDict)
    resScreen.update()
    
    infoScreen = mainScreen.infoScreen
    infoScreen.clear()

def clearSelection(button):
    """
    Restore all values in the properties screen to unchecked and restore the browser
    to all experimental results
    """
    
    propScreen = button.parent
    propItems = propScreen.propItems
    for pitem in propItems:
        if not pitem.is_selected: continue
        if not pitem.valScreen: continue
        for valitem in pitem.valScreen.valItems:
            if not valitem.is_selected:
                continue 
            else:
                valitem.is_selected = False
                valitem.button.deselect()
        pitem.is_selected = False
        pitem.button.deselect()
    
    mainScreen = propScreen.parent    
    mainScreen.browser.loadAllResults()
    mainScreen.resScreen.clear()
    mainScreen.infoScreen.clear()
    mainScreen.dataInfoScreen.clear()

def update_rect(instance,value):
    """
    Auxiliary function to draw background boxes.
    """
    instance.rect.pos = instance.pos
    instance.rect.size = instance.size
