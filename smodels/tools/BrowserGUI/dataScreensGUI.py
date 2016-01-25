#!/usr/bin/env python

from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.button import Button
from kivy.core.window import Window
from kivy.uix.listview import ListItemButton, ListView
from kivy.adapters.listadapter import ListAdapter
from kivy.adapters.simplelistadapter import SimpleListAdapter
Window.clearcolor = [0.08627450980392157, 0.10196078431372549, 0.10588235294117647, 1.0]
from dataObjectsGUI import ValItem, ExpResItem, PropItem, datasetItem, txnameItem

class MainScreen(BoxLayout):
    
    def __init__(self,browser, **kwargs):
        super(MainScreen, self).__init__(**kwargs)
        
        self.orientation = 'horizontal'
        self.padding = [10,10,10,10]
        self.spacing = 15
        self.background_color = (0.21, 0.7, 1.0, 1.)
        self.propScreen = PropScreen(browser)
        self.resScreen = ResScreen(browser)
        self.infoScreen = InfoScreen(browser._selectedExpResults[0])
        self.dataInfoScreen = DataInfoScreen(browser._selectedExpResults[0].datasets[0])
        
        self.size_hint_x = None
        self.propScreen.size_hint = (0.22,1.)
        self.resScreen.size_hint = (0.22,1.)
        self.infoScreen.size_hint  = (0.3,1.)
        self.dataInfoScreen.size_hint  = (0.25,1.)
        self.add_widget(self.propScreen)
        self.add_widget(self.resScreen)
        self.add_widget(self.infoScreen)
        self.add_widget(self.dataInfoScreen)
        
        self.width = 1200
        Window.size = (self.width,list(Window.size)[1])
        
        
        

class PropScreen(BoxLayout):
    
    def __init__(self, browser,**kwargs):
        super(PropScreen, self).__init__(**kwargs)
                        
        self.orientation = 'vertical'
        self.spacing = 10
        self.propItems = []
        max_length = 0
        for attr in sorted(browser.getAttributes()):
            propItem = PropItem(attr,browser)
            propItem.valScreen = ValsScreen(attr,browser)
            self.propItems.append(propItem)
            max_length = max([max_length,len(attr)])
        
        l = Label(text="Properties",font_size=30)
        l.size_hint_y = 0.1
        args_converter = lambda row_index, obj: {'text': obj.name,
                                         'size_hint_y': None,
                                         'height': 25,
                                          'deselected_color' : [0.10980392156862745, 0.12549019607843137, 0.12941176470588237, 1.0],
                                          'selected_color' : [0.34509803921568627, 0.34509803921568627, 0.34509803921568627, 1.0],
                                          'font_size' : 20,                                        
                                         'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.propItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='multiple')
        propList = ListView(adapter=list_adapter)
        propList.children[0].bar_width = 5
        propList.size_hint_y = 0.75
        
        filterButton = Button(text="Select", font_size =25)
        filterButton.size_hint_y = 0.08
        filterButton.bind(on_press=selectResults)
        clearButton = Button(text="Clear Selection", font_size =25)
        clearButton.size_hint_y = 0.06
        clearButton.bind(on_press=clearSelection)
        
        self.add_widget(l)
        self.add_widget(propList)
        self.add_widget(filterButton)
        self.add_widget(clearButton)


class ResScreen(BoxLayout):
    
    def __init__(self,browser, **kwargs):
        super(ResScreen, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = 12
        self.resItems = []
        self.browser = browser
        
    def clear(self):
        self.clear_widgets()
        
    def update(self):
        self.clear()
        self.resItems = []
        for expRes in self.browser:
            self.resItems.append(ExpResItem(expRes,self.browser))
        
                   
        args_converter = lambda row_index, obj: {'text': obj.name,
                                         'size_hint': (1,None),
                                         'height': 40,
                                          'deselected_color' : [0.2, 0.6431372549019608, 0.8156862745098039, 1.0],
                                          'selected_color' : [0.21176470588235297, 0.9307282415630551, 1.0, 1.0],
                                          'font_size' : 20,
                                         'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.resItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='single',propagate_selection_to_data=True)
        resList = ListView(adapter=list_adapter)
        resList.children[0].bar_width = 5
        resList.size_hint_y = 0.9
        
        l = Label(text="Selected:",font_size=30)
        l.size_hint_y = 0.1
                             
        self.add_widget(l)
        self.add_widget(resList)
            

class InfoScreen(BoxLayout):
    
    def __init__(self,expRes, **kwargs):
        super(InfoScreen, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = 12
        self.dataItems = []
            
    def clear(self):
        self.clear_widgets()
        
    def updateWith(self,expRes):                   
        self.clear()
        self.dataItems = []
        
        #First block is global Info
        info_items = []
        for name,val in expRes.globalInfo.__dict__.items():
            info_items.append('[b] [size=17]'+str(name) + '[/size] [/b] : ' + str(val))
        list_adapter = SimpleListAdapter(data=sorted(info_items),
                                     cls=myLabel, selection_mode='single')        
        globalInfoList = ListView(adapter=list_adapter)
        globalInfoList.children[0].bar_width = 5
        
        
        #Second block is dataset list:
        max_length = 0
        for dataset in expRes.datasets:
            self.dataItems.append(datasetItem(dataset))
            max_length = max(max_length,len(str(dataset.dataInfo.dataId))*15.)        
        args_converter = lambda row_index, obj: {'text': str(obj.name),
                                 'size_hint': (None,None),
                                 'height': 40,                                 
                                 'width' : max_length,
                                  'deselected_color' : [0.2, 0.6431372549019608, 0.8156862745098039, 1.0],
                                  'selected_color' : [0.21176470588235297, 0.9307282415630551, 1.0, 1.0],
                                  'font_size' : 20,
                                   'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.dataItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='single')
        datasetList = ListView(adapter=list_adapter)
        datasetList.children[0].bar_width = 5
                
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
    
    def __init__(self,dataset, **kwargs):
        super(DataInfoScreen, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = 12
            
    def clear(self):
        self.clear_widgets()
        
    def updateWith(self,dataset):      
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
        dataInfoList.children[0].bar_width = 5
        
        
        #Second block is txnames list:
        args_converter = lambda row_index, obj: {'text': str(obj.name),
                                 'size_hint': (None,None),
                                 'height': 40,
                                  'deselected_color' : [0.2, 0.6431372549019608, 0.8156862745098039, 1.0],
                                  'selected_color' : [0.21176470588235297, 0.9307282415630551, 1.0, 1.0],
                                  'font_size' : 20,
                                   'on_press' : obj.changeSel}
            
        list_adapter = ListAdapter(data=self.txnameItems,args_converter=args_converter,
                                    cls=ListItemButton, selection_mode='single')
        txnameList = ListView(adapter=list_adapter)
        txnameList.children[0].bar_width = 5
        
        #Third block is Upper Limit getter
        ULgetter = Button(text="Get UL")
        
        dataInfoList.size_hint = (1.,0.3)
        txnameList.size_hint = (None,0.3)
        ULgetter.size_hint = (1.,0.2) 
        
        infoLabel = Label(text="Data Set Info:",font_size=30,size_hint_y=0.07)
        infoLabel.bind(size=infoLabel.setter('text_size')) 
        self.add_widget(infoLabel)
        self.add_widget(dataInfoList)
        txnameLabel = Label(text="TxNames:",font_size=30,size_hint_y=0.08, halign='left')
        txnameLabel.bind(size=txnameLabel.setter('text_size')) 
        self.add_widget(txnameLabel)
        self.add_widget(txnameList)
        ulLabel = Label(text="Compute UL:",font_size=27,size_hint_y=0.05, halign='left')
        ulLabel.bind(size=ulLabel.setter('text_size'))
        self.add_widget(ulLabel)
        self.add_widget(ULgetter)

        

class ValsScreen(BoxLayout):
    
    def __init__(self,text,browser, **kwargs):
        super(ValsScreen, self).__init__(**kwargs)

        self.orientation = 'vertical'
        self.name = text
        self.browser = browser
        self.valItems = []
        
    def create(self):
        if not self.children:   
            #values:
            values = dict([[str(v),v] for v in self.browser.getValuesFor(attribute=self.name)])                       
            for val_str in sorted(values.keys()):
                self.valItems.append(ValItem(self.name,values[val_str]))        
            
            args_converter = lambda row_index, obj: {'text': str(obj.value),
                                              'deselected_color' : [0.2746480740171258, 0.5532196919487817, 0.7886323268206039, 1.0],
                                              'selected_color' : [0.26706971305009763, 0.4634445020575226, 0.6323268206039077, 1.0],
                                              'height' : 30,
                                             'on_press' : obj.changeSel}
                
            list_adapter = ListAdapter(data=self.valItems,args_converter=args_converter,
                                        cls=ListItemButton, selection_mode='multiple')
            valList = ListView(adapter=list_adapter)
            valList.children[0].bar_width = 5
            valList.size_hint_y = 0.9
            
            valLabel = Label(text='Values for %s' %self.name, font_size = 20)
            valLabel.size_hint_y = 0.1
            
            self.add_widget(valLabel)
            self.add_widget(valList)
            
            
class TxnameScreen(BoxLayout):
    
    def __init__(self,txname, **kwargs):
        super(TxnameScreen, self).__init__(**kwargs)

        self.orientation = 'vertical'
        self.txname = txname
        self.infoItems = []
        
    def create(self):
        if not self.children:            
            for name,val in self.txname.__dict__.items():
                if name.lower() == 'efficiencymap' or name.lower() == 'upperlimits': continue
                if name.lower() == 'txnamedata': continue
                if name.lower() == 'globalinfo': continue
                if name[0] == '_': continue 
                self.infoItems.append('[b] [size=17]'+str(name) + '[/size] [/b] : ' + str(val))
            list_adapter = SimpleListAdapter(data=sorted(self.infoItems),
                                         cls=myLabel, selection_mode='single')        
            txnameInfoList = ListView(adapter=list_adapter)
            txnameInfoList.children[0].bar_width = 5
            txnameInfoList.size_hint_y = 0.9
            txnameLabel = Label(text='%s Info:' %self.txname.getInfo('txName'), font_size = 20)
            txnameLabel.size_hint_y = 0.1
            
            self.add_widget(txnameLabel)
            self.add_widget(txnameInfoList)
            
            



class myLabel(Label):    
        def __init__(self,  *args, **kwargs):
            super(myLabel, self).__init__(*args, **kwargs)

            self.halign = 'left'
            self.markup = True
            self.bind(size=self.setter('text_size')) 
            self.shorten = True

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
    resScreen = mainScreen.resScreen
    resScreen.browser.loadAllResults()
    resScreen.clear()    
    infoScreen = mainScreen.infoScreen
    infoScreen.clear()    
