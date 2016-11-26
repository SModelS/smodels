#!/usr/bin/env python

"""
.. module:: tools.BrowserGUI.browserGUI
   :synopsis:  Central GUI module for the database browser (Kivy-based)

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from kivy.logger import Logger
Logger.setLevel(logging.WARNING)

import sys,os
sys.path.append('../')
sys.path.append(os.path.join(os.getenv("HOME"),'smodels'))
from databaseBrowser import Browser
from dataScreensGUI import MainScreen
from kivy.app import App
    
class databaseBrowserGUIApp(App):
    """
    Main Kivy app class. Builds the browser GUI
    """
    
    def __init__(self,browser,**kwargs):
        super(databaseBrowserGUIApp, self).__init__(**kwargs)
        self.browserObj = browser
    
    def build(self):        
        return MainScreen(self.browserObj)    
    
if __name__ == "__main__":
    
    #First check if the database folder exists:
    databasePath = os.path.join(os.getenv("HOME"),'smodels-database/database_light.pcl')
    if  len(sys.argv) > 1:
        databasePath = sys.argv[1]
    if not os.path.isdir(databasePath) and not os.path.isfile(databasePath):
        Logger.error("Database %s not found!" %databasePath)
        sys.exit()
    
    #Load the Browser and check if it is a valid database   
    browser = Browser(databasePath)
    if not browser:
        Logger.error("Error loading database in %s" %databasePath)
        sys.exit()
    
    #Launches the browser GUI    
    app = databaseBrowserGUIApp(browser)
    app.run()    
        

