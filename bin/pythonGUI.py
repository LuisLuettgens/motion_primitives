#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Tkinter import *
import subprocess 
import datetime

guiFont = "DejaVu Sans Mono"
buttonYPad = 10

# Colors
fontColor = '#bca3a3'
backgroundColor = "#333333"
buttonBackgroundColor = "#3f3f3f"
highlightFontColor = "#6f6f6f"
labelFontColor = '#6f6f6f'
buttonActiveFontColor = '#be7474'

# Options
buttonOptions = {'font' :(guiFont, 36), 'fg' : fontColor, 'background' : backgroundColor,
                 'bd' : 0, 'highlightthickness' : 0, 'activebackground' : backgroundColor, 'activeforeground' : buttonActiveFontColor}
frameOptions = {'background' : backgroundColor}
labelOptions = {'font' : (guiFont, 25), 'background' : backgroundColor, 'fg' : labelFontColor}
listBoxOptions = {'bd' : 0, 'highlightthickness' : 0,'background' : buttonBackgroundColor,
                  'font' :(guiFont, 22), 'setgrid' : 5, 'justify' : CENTER}
class AutonomousDrivingGUI:
    def __init__(self, master):
        self.master = master
        master.title("Motion Primitives")

        # Headline
        Label(master, text = 'Motion Primitives', **labelOptions).grid(columnspan = 4, row = 0, column = 0, pady=30)

        # Explanation for the GUI
        explanation = Label(master, text = 'Choose a set of maneuvers and trims, then click on ➤', **labelOptions)
        explanation.configure(font = (guiFont, 12))
        explanation.grid(row = 1, columnspan = 3, sticky = W, padx = 20)
        
        # Knöpfe für das Ausparken
        Label(master, text="Ausparken", **labelOptions).grid(columnspan = 3, sticky = W, padx = 20)
        t1 = '↫'
        Button(master, text=t1, command=lambda: self.addToList(t1), **buttonOptions).grid(row=3, column=1, padx=10, sticky=W)
        t2 = '↬'
        Button(master, text=t2, command=lambda: self.addToList(t2), **buttonOptions).grid(row=3, column=1, padx=10, sticky=E)
        
        # Knöpfe für das Kreuzungsverhalten
        Label(master, text="Kreuzungen", **labelOptions).grid(columnspan = 3, sticky = W, padx = 20)
        t3 = '↰'
        Button(master, text=t3, command=lambda: self.addToList(t3), **buttonOptions).grid(row=5,column=0, padx=20,sticky=E)
        t4 = '↑'                                                                      
        Button(master, text=t4,command=lambda: self.addToList(t4), **buttonOptions).grid(row=5,column=1, padx=20)
        t5 = '↱'
        Button(master, text=t5, command=lambda: self.addToList(t5), **buttonOptions).grid(row=5,column=2, padx=20, sticky=W)
        

        # Knöpfe für das Einparken
        Label(master, text="Einparken", **labelOptions).grid(columnspan = 3, sticky = W, padx = 20)
        t6 = 'P'
        Button(master, text=t6, command=lambda: self.addToList(t6), **buttonOptions).grid(row = 7,column = 1, padx = 20)
        
        
        # List box
        self.listBox = Listbox(master, **listBoxOptions)
        self.listBox.configure(state=DISABLED)
        self.listBox.grid(row = 1, column = 3, rowspan = 6, padx = 50)
        
            
        self.scrollbar = Scrollbar(master, bg = backgroundColor, bd = 0, highlightthickness = 0, troughcolor = backgroundColor, activebackground = buttonBackgroundColor)
        self.listBox.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.listBox.yview)
        self.scrollbar.grid(row = 1, column = 3, rowspan = 6, sticky = N+S+E, padx = 50)


        # Delete button
        delete_button = Button(master, text='⌫', command=self.deleteLastEntry, **buttonOptions)
        delete_button.grid(row = 7, column = 3)

        # Starting button
        buttonOptions['font'] = (guiFont, 50, 'bold')
        run_button = Button(master, text='➤', command= self.runADTF, **buttonOptions)
        run_button.grid(row = 8, columnspan = 4, pady = 60)
    
    def toggle_geom(self,event):
        geom=self.master.winfo_geometry()
        self.master.geometry(self._geom)
        self._geom=geom

    def addToList(self, text):
        self.listBox.configure(state=NORMAL)
        self.listBox.insert(END, text)
        self.listBox.configure(state=DISABLED)

    def deleteLastEntry(self):
        self.listBox.configure(state=NORMAL)
        self.listBox.delete(END, END)
        self.listBox.configure(state=DISABLED)

    def writeManeuverList(self):
        transformedEntries = {u'↫' : 'pull_out_left', u'↬' : 'pull_out_right',
                                  u'↰' : 'left', u'↑' : 'straight', u'↱' : 'right', 'P' : 'cross_parking'}
        now = datetime.datetime.now()
        filename = "maneuverfiles/autonomous-driving_{0}-{1}-{2}_{3}-{4}.xml".format(now.day, now.month, now.year, now.hour, now.minute)

        for name in [filename, "maneuver.xml"]:
            file = open(name,"w")  
            file.write("<?xml version=\"1.0\" encoding=\"iso-8859-1\" standalone=\"no\"?>\n")
            file.write("<AADC-Maneuver-List description=\"Teststrecke\">\n")
            file.write("	<AADC-Sector id=\"0\">\n")
            for (i, entry) in enumerate(self.listBox.get(0, END)):    
                file.write("		<AADC-Maneuver id=\"{0}\" action=\"{1}\" />\n".format(i,str(transformedEntries[entry]))) 
            file.write("	</AADC-Sector>\n")
            file.write("</AADC-Maneuver-List>")
            file.close()        

    def runADTF(self):
        # Write all entries from the list box into a maneuver list
        self.writeManeuverList()

        # Start ADTF
        # subprocess.Popen('/home/aadc/ADTF/config/start_LiveVisualizationConfig.sh', shell=True)        
        
        # End this GUI
        self.master.quit()

        

root = Tk()
root.configure(background = backgroundColor)        
my_gui = AutonomousDrivingGUI(root)
root.mainloop()



