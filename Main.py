
# Main module for DotLattice Experiment 
#
# Copyright (C) 2011 University of Virginia
# Supported by grants to the University of Virginia from the National Science Foundation 
# PI: Prof. Michael Kubovy <kubovy@virginia.edu>
# Author: Yuri Spitsyn <spitsyn@virginia.edu>
#
# version: 0.19.24
# last modified: 03.21.2011
#
# Distributed under the terms of the GNU Lesser General Public License (LGPL). 
#

from psychopy import visual, core, event #import some libraries from PsychoPy

from DLVisual import *
from DLSession import *

import pyglet.window

allScrs = pyglet.window.get_platform().get_default_display().get_screens()
defScr = allScrs[0]
winSize = [defScr.width, defScr.height]

#
# Create experiment window
#
expWin = visual.Window(size=winSize, units="pix") #create a window

#
# Create experiment session and configure it from the data file
#
ss = DLSession(expWin, "expdata.txt")
# ss.printExpdata()

#
# Create dot lattice
#
lattice = DotLattice(win=expWin, session=ss, dotVariantsList=[[ss.dotLum,1]])

if ss.isSessionValid():

    #
    # Create fixation
    #
    fixation = FixationScreen(win=expWin, session=ss)

    #
    # Create response screen
    #
    response = ResponseScreen(win=expWin, session=ss)

    #
    # Practice
    #
    if not ss.isDemo():
    
        ss.showInstructions()

        ss.startPractice()

        while not ss.isOver():

            ss.presentStimulus()
            ss.showResponseOptions()
            ss.finalizeTrial()

    #        
    # Start demo or experiment
    #
    ss.startSession()        

    while not ss.isOver():

        ss.presentStimulus()

        if ss.isDemo():
            ss.showDemoResponse()
        else:    
            ss.showResponseOptions()
            
        ss.finalizeTrial()

    ss.finalizeSession()
    
else:
    
    ss.showText("Session cannot proceed: not enough data!")
#
# Bye
#
core.quit()    
        

