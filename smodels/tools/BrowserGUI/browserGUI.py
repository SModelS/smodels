#!/usr/bin/env python

"""
.. module:: tools.BrowserGUI.browserGUI
   :synopsis:  Central GUI module for the database browser (Kivy-based)

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import sys,os
sys.path.append('../')
sys.path.append(os.path.join(os.getenv("HOME"),'smodels'))
from databaseBrowser import Browser
from dataScreensGUI import MainScreen
from kivy.app import App
import logging
logger = logging.getLogger(__name__)
    
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
    databasePath = os.path.join(os.getenv("HOME"),'smodels-database')
    if  len(sys.argv) > 1:
        databasePath = sys.argv[1]
    if not os.path.isdir(databasePath):
        logger.error("Database folder %s does not exist!" %databasePath)
        sys.exit()
     
    #Load the Browser and check if it is a valid database   
    browser = Browser(databasePath)
    if not browser:
        logger.error("Error loading database in %s" %databasePath)
        sys.exit()
    
    #Launches the browser GUI    
    app = databaseBrowserGUIApp(browser)
    app.run()    
        

