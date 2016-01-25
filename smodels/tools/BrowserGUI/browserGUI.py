#!/usr/bin/env python

import sys,os
sys.path.append('../')
sys.path.append(os.path.join(os.getenv("HOME"),'smodels'))
from databaseBrowser import Browser
from dataScreensGUI import MainScreen
from kivy.app import App
    
class browserApp(App):
    
    def __init__(self,browser,**kwargs):
        super(browserApp, self).__init__(**kwargs)
        self.browserObj = browser
    
    def build(self):        
        return MainScreen(self.browserObj)    
    
if __name__ == "__main__":
    
#First load the browser:
#     browser = Browser(os.path.join(os.getenv("HOME"),'smodels-database'))
    browser = Browser('/home/lessa/smodels/test/database')
    app = browserApp(browser)
    app.run()    
        

