#!/usr/bin/env python

import sys,os

from kivy.uix.colorpicker import ColorPicker
from kivy.app import App
from kivy.uix.floatlayout import FloatLayout


# To monitor changes, we can bind to color property changes
def on_color(instance, value):
    print "RGBA = ", str(value)  #  or instance.color
    print "HSV = ", str(instance.hsv)
    print "HEX = ", str(instance.hex_color)



class clrpickerApp(App):
    
    def __init__(self,**kwargs):
        super(clrpickerApp, self).__init__(**kwargs)
    
    def build(self):
	f = FloatLayout()
        clr_picker = ColorPicker() 
	clr_picker.bind(color=on_color)
        f.add_widget(clr_picker)
        return f
    
    
if __name__ == "__main__":
    
    app = clrpickerApp()
    app.run()    
