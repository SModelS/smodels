#!/usr/bin/env python

from kivy.uix.label import Label
from kivy.adapters.models import SelectableDataItem
from kivy.uix.modalview import  ModalView

class PropItem(SelectableDataItem):
    
    def __init__(self, text, browser):
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
            max_length = 10
            for val in self.valScreen.valItems:
                max_length = max(max_length,len(str(val.value)))
            self.popup = ModalView(auto_dismiss=True)            
            self.popup.add_widget(self.valScreen)            
            self.popup.size_hint=(None,0.8)         
            self.popup.width = min(max_length*12.,0.9*button.parent.parent.parent.parent.parent.width)
            self.popup.background =  "atlas://data/images/defaulttheme/sliderh_background"

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
    
    def __init__(self, expRes, browser):        
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
            mainScreen.infoScreen.updateWith(self.expRes)
        else:
            mainScreen.infoScreen.clear()

class ResItem(SelectableDataItem):
    
    def __init__(self, expRes):        
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
    
    def __init__(self, dataset):        
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
    
    def __init__(self, txname):        
        self.dataset = txname
        self.name = txname.getInfo('txName')
        self.is_selected = False        
        self.popup = None
        self.txnameScreen = None        
    
    def changeSel(self,button):
        if not self.popup:
            self.button = button
            self.txnameScreen.create()           
            max_length = 10
            for val in self.txnameScreen.infoItems:
                max_length = max(max_length,len(str(val)))
            self.popup = ModalView(auto_dismiss=True)            
            self.popup.add_widget(self.txnameScreen)            
            self.popup.size_hint=(None,0.8)         
            self.popup.width = min(max_length*8.,0.9*button.parent.parent.parent.parent.parent.width)
            print type(button.parent.parent.parent.parent.parent),button.parent.parent.parent.parent.parent.width
            print self.popup.width
            self.popup.background =  "atlas://data/images/defaulttheme/sliderh_background"
           
        self.popup.open()


class ValItem(SelectableDataItem):
    
    def __init__(self, text, val):
        self.name = text
        self.value = val
        self.is_selected = False
        self.button = None
        
    def changeSel(self,button):
        self.button = button           
        self.is_selected = not self.is_selected
