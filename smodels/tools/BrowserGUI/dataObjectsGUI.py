#!/usr/bin/env python

from kivy.uix.label import Label
from kivy.adapters.models import SelectableDataItem
from kivy.uix.modalview import  ModalView
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.uix.boxlayout import BoxLayout

class PropItem(SelectableDataItem):
    
    def __init__(self, text, browser,**kwargs):
        super(PropItem, self).__init__(**kwargs)
        self.name = text
        self.is_selected = False
        self.broswer = browser        
        self.popup = None
        self.valScreen = None
        self.button = None
        
    def changeSel(self,button):
        if not self.popup:
            self.button = button
            self.valScreen.create()           
            max_length = 20
            for val in self.valScreen.valItems:
                max_length = max(max_length,len(str(val.value)))
            self.popup = ModalView(auto_dismiss=True)            
            self.popup.add_widget(self.valScreen)            
            self.popup.size_hint=(None,None)         
            self.popup.width = min(max_length*12.,0.9*button.parent.parent.parent.parent.parent.width)
            self.popup.height = min(len(self.valScreen.valItems)*30.+40.
                                    ,0.9*button.parent.parent.parent.parent.parent.height)
            

            self.popup.bind(on_dismiss=self.checkSelection)                
        self.popup.open()
    
    def checkSelection(self,button):
        
        self.is_selected = False
        for val in self.valScreen.valItems:
            if val.is_selected:
                self.is_selected = True
                break
                
        if self.is_selected:
            self.button.select()
        else:
            self.button.deselect()
        
            
class ExpResItem(SelectableDataItem):
    
    def __init__(self, expRes, browser,**kwargs):
        super(ExpResItem, self).__init__(**kwargs)
         
        self.expRes = expRes
        self.name = expRes.globalInfo.id
        self.is_selected = False
        self.broswer = browser        
        self.popup = None
        
    def changeSel(self,button):   
        resScreen = button.parent.parent.parent.parent
        mainScreen = resScreen.parent
        if not self.is_selected:
            mainScreen.infoScreen.clear()
            mainScreen.dataInfoScreen.clear()
            mainScreen.infoScreen.updateWith(self.expRes)
        else:
            mainScreen.infoScreen.clear()
            mainScreen.dataInfoScreen.clear()

class ResItem(SelectableDataItem):
    
    def __init__(self, expRes,**kwargs):
        super(ResItem, self).__init__(**kwargs)
        
        self.expRes = expRes
        self.name = expRes.globalInfo.id
        self.is_selected = False        
        self.popup = None
        
    def changeSel(self,button):   
        resScreen = button.parent.parent.parent.parent
        mainScreen = resScreen.parent
        if not self.is_selected:
            mainScreen.infoScreen.clear()
            mainScreen.infoScreen.updateWith(self.expRes)
        else:
            mainScreen.infoScreen.clear()


class datasetItem(SelectableDataItem):
    
    def __init__(self, dataset,**kwargs): 
        super(datasetItem, self).__init__(**kwargs)
               
        self.dataset = dataset
        self.name = dataset.dataInfo.dataId
        self.is_selected = False        
        self.popup = None
        
    def changeSel(self,button):   
        resScreen = button.parent.parent.parent.parent
        mainScreen = resScreen.parent
        if not self.is_selected:
            mainScreen.dataInfoScreen.clear()
            mainScreen.dataInfoScreen.updateWith(self.dataset)
        else:
            mainScreen.dataInfoScreen.clear()
            

class txnameItem(SelectableDataItem):
    
    def __init__(self, txname,**kwargs):
        super(txnameItem, self).__init__(**kwargs)
                
        self.txname = txname
        self.name = txname.getInfo('txName')
        self.is_selected = False 
        self.popup = None
        self.txnameScreen = None
        self.ulGetter = None
    
    def changeSel(self,button):
        mainScreen = button.parent
        while not 'MainScreen' in str(type(mainScreen)):
            mainScreen = mainScreen.parent
        if not self.popup:
            self.button = button
            self.txnameScreen.create()           
            max_length = 10
            for val in self.txnameScreen.infoItems:
                max_length = max(max_length,len(str(val)))
            self.popup = ModalView(auto_dismiss=True)            
            self.popup.add_widget(self.txnameScreen)            
            self.popup.size_hint=(None,0.8)         
            self.popup.width = min(max_length*8.,0.9*mainScreen.width)            
        
        if not self.is_selected:            
            mainScreen.dataInfoScreen.ulGetter.updateWith(self.txname)
            self.popup.open()
        else:
            mainScreen.dataInfoScreen.ulGetter.clear()


class ValItem(SelectableDataItem):
    
    def __init__(self, text, val, **kwargs):
        super(ValItem, self).__init__(**kwargs)
        
        self.name = text
        self.value = val
        self.is_selected = False
        self.button = None
        
    def changeSel(self,button):
        self.button = button
        self.is_selected = not self.is_selected
             
class ULgetter(BoxLayout):
    """
    Defines the box for the upper limit calculator
    """
    
    def __init__(self,txname,**kwargs):
        super(ULgetter, self).__init__(**kwargs)
        
        self.orientation = 'vertical'
        self.spacing = 5
        self.padding = 2
        self.txname = txname
        ulbutton = Button(text="Get UL")
        ulbutton.bind(on_press=self.getUL)
        massLabel  = Label(text="Mass (GeV):",font_size=20)        
        massInput = TextInput(text=str(txname.txnameData.data[0][0]),focus=True)
        massInput.multiline = False
        massInput.bind(on_text_validate=self.getUL)
        resLabel  = Label(text="Upper Limit (fb):",font_size=20)        
        self.button = ulbutton
        self.massInput = massInput
        self.resLabel =resLabel
        self.add_widget(massLabel)
        self.add_widget(massInput)
        self.add_widget(resLabel)
        self.add_widget(ulbutton)
        
    def getUL(self,button=None):        
        try:
            ul = self.txname.txnameData.getValueFor(eval(self.massInput.text))
        except:
            ul = False                
        if not ul:
            self.resLabel.text = 'Error computing limit'
        else:
            self.resLabel.text = "Upper Limit  = " + str(ul)
