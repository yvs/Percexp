

# DLSession  
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

import os
import csv
import time
import wx

from psychopy import gui, visual, core, event

from DLVisual import *



class DLSession:
    """
    Class representing an experiment session.
    Holds all the information about the session details
    and controls all session components: fixation, lattice, user response etc.
    After each trial it saves report to the output file and
    prepares the stimulus data for the next trial
    """
    def __init__(self, win, dataFileName, trialsCount='auto'):

        #
        # Session-wide parameter identifiers
        #
        self.expName='experimentName'
        self.workDir='workingDirectory'
        self.instrFile='instructions'
        self.sessId='sessionId'
        self.subjId='subjectId'
        self.asstId='assistantId'
        self.paid='paid'
        self.demo='demoMode'
        self.demoTrials='demoTrials'
        self.pracTrials='practiceTrials'
        self.sessTrials='sessionTrials'
        self.maskRad='maskingRadius'
        self.maskDur='maskingDuration'
        self.demoStimDur='demoStimulusDuration'
        self.fixDur='fixationDuration'
        self.pracStimDur='practiceStimulusDuration'
        self.trialStimDur='trialStimulusDuration'
        self.intertrialDur='intertrialDuration'
        self.respBgnd='responseBgnd'
        self.respSize='responseChoiceDiameter'
        self.respDist='responseChoicesDistance'
        self.winBgnd='windowBgnd'
        self.aperBgnd='apertureBgnd'
        self.aperSize='apertureSize'
        self.aperPos='aperturePos'
        self.aperShow='showAperture'
        self.aDist='aDistance'
        self.baAspect='baAspect'
        self.gamma='gamma'
        self.theta='theta'
        self.dotSize='dotSize'
        self.dotAngle='dotAngle'
        self.dotLum='dotLum'
        self.dotAlter='dotAlternation'
        self.dotTex='dotTex'
        self.dotMask='dotMask'

        #
        # Session's dictionary of parameters
        #
        self.expdata = {
           self.expName: [{ 'default': "DotLattice", 'reportAs':'experiment' }],
           self.workDir: [{ 'default': "/Users/Shared" }],
           self.instrFile: [{ 'default': "instructions.txt" }],
           self.sessId: [{ 'default': 1, 'reportAs':'session' }],
           self.subjId: [{ 'default': 1, 'reportAs':'subject' }],
           self.asstId: [{ 'default': 1, 'reportAs':'assistant' }],
           self.paid: [{ 'default': False }],
           self.demo: [{ 'default': False }],
           self.pracTrials: [{ 'default': 5 }],
           self.demoTrials: [{ 'default': 5 }],
           self.sessTrials: [{ 'default': 'auto' }],
           self.maskRad: [{ 'default': 30, 'reportAs':'maskRadius' }],
           self.fixDur: [{ 'default': 0.3 }],
           self.maskDur: [{ 'default': 0.21, 'reportAs':'maskDur' }],
           self.demoStimDur: [{ 'default': 5 }],
           self.pracStimDur: [{ 'default': 0.5, 'reportAs':'practiceStimDur' }],
           self.trialStimDur: [{ 'default': 0.3, 'reportAs':'trialStimDur' }],
           self.intertrialDur: [{ 'default': 0.58 }],
           self.respSize: [{ 'default': 160 }],
           self.respDist: [{ 'default': 180 }],
           self.respBgnd: [{ 'default': -1 }],
           self.winBgnd: [{ 'default': -1, 'reportAs':'winBgnd' }],
           self.aperBgnd: [{ 'default': 0 }],
           self.aperSize: [{ 'default': 600 }],
           self.aperPos: [{ 'default': (0,0) }],
           self.aperShow: [{ 'default': True }],
           self.aDist: [{ 'default': 30 }],
           self.baAspect: [{ 'default': 1.0 }],
           self.gamma: [{ 'default': 90 }],
           self.theta: [{ 'default': 0 }],
           self.dotAngle: [{ 'default': 0 }],
           self.dotLum: [{ 'default': 0.5, 'reportAs':'lum' }],
           self.dotSize: [{ 'default': 20 }],
           self.dotAlter: [{ 'default': 'a', 'reportAs':'alternation' }],
           self.dotTex: [{ 'default': 'flat' }],
           self.dotMask: [{ 'default': 'log2' }]
        }
        
        #
        # Parameter names to include in report
        #

        # 1 - general session parms
        self._cmhdr = [[
                self.pracTrials,
                self.sessTrials,
                self.sessId,
                self.subjId,
                self.asstId,
                self.paid,
                self.expName
                ]]
        
        # 2 - stimulus presentation related session parms
        self._sshdr = [[
                self.aperSize,
                self.winBgnd,
                self.aperBgnd,
                self.pracStimDur,
                self.trialStimDur,
                self.fixDur,
                self.maskDur,
                self.maskRad
                ]]
        
        # 3 - trial's common lattice header
        self._cthdr = [
               'trial',
               'response',
               'responseTime',
               'order',
               'theta', 
               'baAspect',
               'gamma',
               'lenA', 
               'angleA', 
               'lenB', 
               'angleB', 
               'lenC', 
               'angleC', 
               'lenD', 
               'angleD'
               ]
        
        # 4 - trial specific parameters
        self._sthdr = [[
               'dotSize',
               'dotLum1',
               'dotLum2',
               'dotAlternation'],
               [
               'dotSize',
               'lum1',
               'lum2',
               'alternation'
               ]]
        # lists to hold session related values for report (subset of self.expdata)
        self._cmval=[]
        self._ssval=[]

        #
        # Parameters to be set from the configuration file
        #
        self.blocks= {}         # Block sequences
        self.trialUpdates=[]    # Names of parameters to be updated per trial (by controller)
        # interpret data file
        self._parseSessionDataFile(dataFileName)
        
        #
        # Parameters to be set from the dialog window
        #
        self.fromGui = [
              { 'label':'' },
              { 'label':'Session Id' },
              { 'name':self.sessId, 'label':'Session Id' },
              { 'name':self.subjId, 'label':'Subject Id' },
              { 'name':self.asstId, 'label':'RA Id' },
              { 'label':'' },
              { 'name':self.demo, 'label':'Demo mode:' },
              { 'label':'' }
              ]
        self.runSetupDialog()

        #
        # Practice related data
        #
        self.pracdata={}
        self.practiceUpdates=[]
        # create some practice data
        self._generatePracticeData()
        
        self._practiceMode=False

        self._updates=None      # either trialUpdates or practiceUpdates
        
        #
        # Session window
        #
        self.win=win
        self.mouse= event.Mouse(win=win)
        
        #
        # Lattice (to be set later)
        #
        self.lattice=None
        self.latticeController=None

        #
        # Trials
        #
        if trialsCount != 'auto':
            if self.isDemo(): parm = self.demoTrials
            else:             parm = self.sessTrials
            self._setExpdata(parm, trialsCount)
        #internal trial counters
        self._trialsToGo=0
        self._curTrial=-1

        #
        # Various durations
        #
        self._stimWait=0 # wait period for stimulus (varies depending on self._practiceMode)
        self._maskWait=self.getExpdata(self.maskDur)
        self._interWait=self.getExpdata(self.intertrialDur)
        self.interTrialDur = self.getExpdata(self.intertrialDur)
        self.fixWait=self.getExpdata(self.fixDur)

        #
        # Fixation
        #
        self.fixation=None

        #
        # Response
        #
        self.responseObj=None
        self.responseBgnd=self.getExpdata(self.respBgnd)

        #
        # Instructions
        #
        self.instructions=''
        self._readInstructions()

        # text helper
        self.anykey = visual.TextStim(win, text="Press any key to continue...",
                                      pos=(0,-20), height=20, wrapWidth=400)

        # clocking
        self.clock = core.Clock()
        self.responseDur=-1

        #
        # Output related
        #
        self._resolveSessionHeaders()
        self._resolveSessionValues()
        self.sessionValues=None

        dataFile=self.getExpdata(self.expName)
        dataFileTemp = dataFile + "-" + time.strftime('%Y%m%d%H%M%S') + '.csv'
        self.dataPath=self.getExpdata(self.workDir)
        if self.dataPath[-1] != '/': self.dataPath += '/'
        self.dataPathTemp = self.dataPath + dataFileTemp
        self.dataPath += dataFile + '.csv'

        self.trialFile=None
        self.trialWriter=None

        


    """

    VISUAL

    """
    
    def runSetupDialog(self):
        t = self.getExpdata(self.expName) + '_experiment'
        dlg = gui.Dlg(title=t)
        for fld in self.fromGui:
            if not 'name' in fld.keys():
                dlg.addText(fld['label'])
            else:
                val=self.getExpdata(fld['name'])
                dlg.addField(fld['label'],val)
        dlg.show()
        if gui.OK:
            i=0
            for fld in self.fromGui:
                if 'name' in fld.keys():
                    self._setExpdata(fld['name'],dlg.data[i])
                    i += 1
        else:
            print "runSetupDialog(): user canceled."
                             

    def presentStimulus(self):
        
        self.mouse.setVisible(0)

        self.win.setColor(self.lattice.getBgndColor())
        self.win.flip()

        self.fixation.draw()
        self.win.flip()
        core.wait(self.fixWait)

        self.lattice.draw()
        self.win.flip()
        core.wait(self._stimWait)

        self.lattice.drawMask(self._maskWait)
        

    def showResponseOptions(self):
        
        self.responseObj.updateVariants()
        
        self.win.setColor(self.responseBgnd)
        self.win.flip()

        t=self.clock.getTime()
        self.responseDur=-1
        resp=None

        while resp==None:

            pos=self.mouse.getPos()
            choice=self.responseObj.variantForPoint(pos)
            
            self.mouse.setVisible(1)

            self.responseObj.draw()
            self.win.flip()

            btn=self.mouse.getPressed()
            if btn[0]==True and choice != -1:

                resp=True
                self.responseDur=self.clock.getTime()-t
                self.responseObj.blinkSelected()
                self.win.flip()
                
            event.clearEvents()        
        return choice


    def showDemoResponse(self):
        resp=None
        while resp==None:
            pos=self.mouse.getPos()
            self.lattice.basisVisible(self.lattice.pointInsideAperture(pos))
            self.lattice.draw()
            self.mouse.setVisible(1)
            self.win.flip()
            btn = self.mouse.getPressed()
            if btn[0] == True:
                resp=True
            allKeys=event.getKeys()
            for key in allKeys:
                if key in ['q', 'escape']:
                    core.quit()
            event.clearEvents()
        self.lattice.basisVisible(False)    
    
    def showInstructions(self):
        self.showText(self.instructions)

    def showText(self, text, autoPause=None):
        if self.win.winType == 'pyglet':
            winSize = self.win.winHandle.get_size()
        else: winSize = (1024,768) #not very clever; TO_DO: change when there is a time
        #print "winSize"
        #print winSize
        #print text
        texS = visual.TextStim(self.win, text=text, height=32, wrapWidth=winSize[0]-100)
        texS.setPos((0,texS.height))
        self.win.setColor(-1)
        self.win.flip()
        texS.draw()
        if autoPause==None:
            self.anykey.draw()
            self.win.flip()
            event.waitKeys()
        else:
            self.win.flip()
            core.wait(autoPause)
        self.win.flip()


    """

    SESSION OBJECTS ACCESS

    """
        
    def setFixation(self, fix):
        self.fixation=fix
        
    def setResponseObject(self, resp):
        self.responseObj=resp
        
    def getLattice(self):
        return self.lattice
    
    def setLattice(self, lattice):

        # set lattice and create controller for it
        #
        self.lattice=lattice
        self.latticeController=LatticeController(self.lattice)

        # create data for controller
        #

        #1 - block sequences
        for n, parms in self.blocks.iteritems():
            self.latticeController.createBlockSequence(n, parms['parnames'])
            self.latticeController.addBlocksToSequence(n, parms['blocks'])
            self.latticeController.finalizeSequence(n, parms['repeats'], parms['shuffled'])

    
    """

    SESSION FLOW CONTROL

    """
    
    def startPractice(self):
        self._practiceMode=True
        self._configureVaryingControl()
        self._curTrial=-1
        self._trialsToGo=int(self.getExpdata(self.pracTrials))
        self._stimWait=self.getExpdata(self.pracStimDur)
        self.latticeController.rewindAll()
        self.advanceTrialData()
        s = 'You will start with ' + str(self._trialsToGo) + ' practice trials'
        self.showText(s)

    def startSession(self):
        self._practiceMode=False
        self._configureVaryingControl()
        self._curTrial=-1
        self._trialsToGo=self.getSessionTrialsCount()
        if self.isDemo(): parm = self.demoStimDur
        else:             parm = self.trialStimDur
        self._stimWait=self.getExpdata(parm)
        self.latticeController.rewindAll()
        self.advanceTrialData()
        
        if not self.isDemo():
            if not os.path.exists(self.dataPath):
                self.writeReportFileHeader()
            #create data writer object
            self.trialFile = open(self.dataPathTemp, 'wb')
            self.trialWriter=csv.writer(self.trialFile,delimiter='\t')
            self.sessionValues = self.getSessionValues()
            self.showText('Experiment will now begin.')
        else:    
            s = 'You will start with ' + str(self._trialsToGo) + ' demo trials'
            self.showText(s)

    def finalizeTrial(self):

        if not (self.isDemo() or self.isPractice()):
            self.saveTrialData()
            
        self.advanceTrialData()
        
        core.wait(self._interWait)

    def finalizeSession(self):
        if not self.isDemo():
            self.trialFile.close() #close trial output file first
            #copy session data to the cumulative file
            fo = open(self.dataPath, 'a')
            fi = open(self.dataPathTemp, 'r')
            content = fi.readlines()
            fo.writelines(content)
            fi.close()
            fo.close()
            #say bye
            self.showText('Experiment complete. Thanks for your participation.',2)
        else:    
            self.showText('Demo complete.',2)
            
    def advanceTrialData(self):
        #update separately varying parameters
        self.latticeController.nextInControl(self._updates)
        if not self._practiceMode:
            #update block sequence(s)
            for block in self.blocks.iterkeys():
                self.latticeController.nextInSequence(block)
        self._curTrial += 1        
    
    def isOver(self):
       return not self._curTrial<self._trialsToGo


    """

    SESSION DATA ACCESS

    """

    def isSessionValid(self):
        return 0 < self.getSessionTrialsCount()

    def isDemo(self):
        return self.getExpdata(self.demo)

    def isPractice(self):
        return self._practiceMode

    def getSessionTrialsCount(self):
        if self.isDemo(): parm = self.demoTrials
        else:             parm = self.sessTrials
        tcnt = self.getExpdata(parm)
        if tcnt=='auto':
            tcnt = self.latticeController.getLongestSequenceLength()
        if tcnt==0:
            print "getSessionTrialsCount(): returned 0"
        return int(tcnt)    
    
    def _readInstructions(self):
        try:
            fname = self.getExpdata(self.instrFile)
            f = open(fname,'r')
            text = f.readlines()
            for line in text:
                self.instructions += line
        except:
            self.instructions= "Instructions file isn't found: %s" % fname
        
    def _generatePracticeData(self):
        self.pracdata[self.theta]={ 'stepCount':10 }
        self.pracdata[self.gamma]={ 'stepCount':10 }
        self.pracdata[self.baAspect]={ 'parMax':2.0, 'stepCount':10 }
        self.practiceUpdates=[ self.theta, self.gamma, self.baAspect]
        
    def _configureVaryingControl(self):
        self.latticeController.resetControl()
        if self._practiceMode==True:
            for p, v in self.pracdata.iteritems():
                self.latticeController.setControl(p, v)
            self._updates=self.practiceUpdates
        else:
            for n in self.trialUpdates:
                ndx, p = self._describeParm(n)
                try:
                    d = self.expdata[p][ndx]
                    try:
                        parms = d['varying']
                        self.latticeController.setControl(n, parms)
                    except: print "Varying parm %s has no associated data" % n
                except: print "Unknown parameter name: %s." % p
            self._updates=self.trialUpdates    
            
    def _setExpdata(self, dataName, dataVal, dataVariant=0):
        try:
            data = self.expdata[dataName]
            if dataVariant<len(data):
                data[dataVariant]['value']=dataVal
            else:
                print "_setExpdata(): dataVariant index is out of bounds: %d" % dataVariant
        except:
            print "_setExpdata(): dataName not found: %s" % dataName
            
    def getExpdata(self, dataName, ndx=0, ndxFailSafe=False):
        rv = None
        try:
            parm = self.expdata[dataName]
            try:
                d = parm[ndx]
                try:
                    rv = d['value']
                    if rv=='varying':
                        rv = d['default']
                except: rv = d['default']
            except:
                if ndxFailSafe == True:
                    self._expandData(dataName, ndx)
                    rv = self.getExpdata(dataName, ndx)
                else:
                    print "getExpdata(): wrong variant index (%d) for parameter %s" % (ndx, dataName)
        except:
            print "getExpdata(): unknown data name: %s" % dataName
        return rv

    def _expandData(self, dataName, ndx):
        parm = self.expdata[dataName]
        d = parm[0]
        v = {'default': d['default']}
        for i in range(1 + ndx - len(parm)):
            self.expdata[dataName].append(v)

    def getExpdataRaw(self, name):
        return self.expdata[name]

    def _describeParm(self, parmName):
        sn = parmName[-1]
        if sn.isdigit():
            n=int(sn)-1
            p = parmName[:-1]
        else:
            n=0
            p = parmName
        return n, p

    def _parseStringVal(self, sv):
        SV = sv.upper()
        if sv[0]=='(' or sv[0]=='[':
            sv = sv.strip("([])")
            splsv = sv.split(',')
            rv = []
            for s in splsv:
                v = self._parseStringVal(s)
                rv.append(v)
        elif SV == 'TRUE' or SV=='YES':
            rv = True
        elif SV == 'FALSE' or SV=='NO':
            rv = False
        elif SV == 'NONE':
            rv = None
        else:
            try:
                rv = float(sv)
                if -1 == sv.find('.'):
                    rv = int(sv)
            except: rv = sv
        return rv

    def _parseSessionDataFile(self, dataFileName):

        try:
            
            f = open(dataFileName,'r')
            text=f.readlines()
            blockBegin = False
            mode = None

            for i in range(len(text)):
                line = text[i].strip()
                if line != "":
                    if line[0]=='#': pass
                    else:
                        spl = line.split('#')
                        toparse = spl[0].strip()
                        if mode == None:
                            #opt.upper()
                            if -1 != toparse.find("<SESSION>"):
                                mode = 'session'
                            elif -1 != toparse.find("<FIXED>"):
                                mode = 'fixed'
                            elif -1 != toparse.find("<VARYING>"):
                                mode = 'varying'
                            elif -1 != toparse.find("<BLOCK"):
                                op = toparse[:-1]
                                mode = 'block'
                                seqhdr = op.split()
                                repeats=1
                                shuffled=True
                                for s in seqhdr:
                                    parm = s.split('=')
                                    if parm[0] == 'repeats':
                                        repeats = int(parm[1])
                                    elif parm[0] == 'shuffled':
                                        shuffled = self._parseStringVal(parm[1])
                                n = str(len(self.blocks)+1)
                                blockname = 'block'+n        
                                self.blocks[blockname]={ 'repeats': repeats, 'shuffled': shuffled }
                                self.blocks[blockname]['blocks']=[]
                                blockBegin=True
                        elif mode == 'block':
                            if -1 != toparse.find("</BLOCK>"):
                                mode = None
                            elif blockBegin == True:
                                p = toparse.split()
                                pnames = []
                                for name in p: 
                                    pnames.append(name)
                                self.blocks[blockname]['parnames']=pnames
                                blockBegin = False
                            else:
                                b = toparse.split()
                                blk = []
                                for v in b:
                                    blk.append(self._parseStringVal(v))
                                self.blocks[blockname]['blocks'].append(blk)
                        elif mode == 'varying':
                            if -1 != toparse.find("</VARYING>"): mode = None
                            else:
                                p = toparse.split()
                                ndx, parname = self._describeParm(p[0])
                                try:
                                    d = self.expdata[parname][ndx]
                                    d['value']='varying'
                                    v = p[1:] #remove parm name from further parsing
                                    dic = {}
                                    for field in v:
                                        try:
                                            ff = field.split('=')
                                            val = self._parseStringVal(ff[1])
                                            dic[ff[0]] = val
                                        except: print "Problem with varying parameter field: %s" % field
                                    d['varying']=dic    
                                    self.trialUpdates.append(p[0])
                                except: print "VARYING: Unknown parameter name: %s." % parname    
                        elif mode == 'fixed': 
                            if -1 != toparse.find("</FIXED>"): mode = None
                            else:
                                p = toparse.split()
                                ndx, parname = self._describeParm(p[0])
                                #make sure this parm name is known
                                try:
                                    parlist=self.expdata[parname]
                                    val = self._parseStringVal(p[1])
                                    ln = len(parlist)
                                    if ndx<ln:
                                        d = parlist[ndx]
                                        d['value']=val
                                    elif ndx==ln:
                                        d = { 'value': val }
                                        parlist.append(d)
                                    else:    
                                        print "FIXED: parameter index is too big: %d." % ndx
                                except: print "FIXED: Unknown parameter name: %s." % parname    
                        elif mode == 'session': 
                            if -1 != toparse.find("</SESSION>"): mode = None
                            else:
                                p = toparse.split()
                                parname = p[0]
                                #make sure this parm name is known
                                try:
                                    d = self.expdata[parname][0]
                                    d['value']= self._parseStringVal(p[1])
                                except:
                                    print "SESSION: Unknown parameter name: %s." % parname    
                        else:
                            pass
            f.close()
        except:
            print "_parseSessionDataFile(): Can't open the input file: %s" % dataFileName 


    def printExpdata(self):
        print "Printing content of the session's expdata dictionary:"
        for k, v in self.expdata.iteritems():
            print k
            print v

    """

    REPORTS

    """
    

    def _resolveSessionValues(self):
        # set session values for report
        warn="_resolveSessionValues(): Unknown parameter name: %s"    
        for name in self._cmhdr[0]:
            if name in self.expdata.keys():
                val = self.getExpdata(name)
                self._cmval.append(val)
            else: print warn % name
        for name in self._sshdr[0]:
            if name in self.expdata.keys():
                val = self.getExpdata(name)
                self._ssval.append(val)
            else: print warn % name
        
        
    def _resolveSessionHeaders(self):
        # set session header names for report
        for i in range(2):
            if i==0: lst=self._cmhdr
            else:    lst=self._sshdr
            hdrs=[]
            for h in lst[0]:
                name=h
                if h in self.expdata.keys():
                    d = self.expdata[h][0]
                    #print h
                    #print d
                    if 'reportAs' in d.keys():
                        name = d['reportAs']
                hdrs.append(name)
            lst.append(hdrs)
            #print lst

    def saveTrialData(self):
        timeStamp=str( time.ctime(time.time()))
        dataValues = self.getTrialValues() + self.sessionValues + [timeStamp]
        self.trialWriter.writerow(dataValues)
        

    def writeReportFileHeader(self):
        dataHeader = self.getTrialHeader() + self.getSessionHeader() + ['timeStamp']
        fo = open(self.dataPath,'wb')
        writer = csv.writer(fo, delimiter='\t')
        writer.writerow(dataHeader)
        fo.close()
        
    def getSessionHeader(self, hdr='all'):
        if hdr=='all': return self._sshdr[1] + self._cmhdr[1]
        elif hdr=='common': return self._cmhdr[1]
        elif hdr=='stimulus': return self._sshdr[1]

    def getSessionValues(self, val='all'):
        if val=='all': return self._ssval + self._cmval
        elif val=='common': return self._cmval
        elif val=='stimulus': return self._ssval

    def getTrialHeader(self, hdr='all'):
        if hdr=='all': return self._cthdr + self._sthdr[1]
        elif hdr=='common': return self._cthdr
        elif hdr=='specific': return self._sthdr[1]

    def getTrialValues(self, valmode='all'):
        valcmn=[]
        valspec=[]
        if valmode=='all' or valmode=='common':
            for n in self._cthdr:
                if n=='trial': val=self._curTrial+1 #since it's zero bazed
                elif n=='response': val=self.responseObj.getCurVariantName()
                elif n=='order': val=self.responseObj.getVariantsOrder()
                elif n=='responseTime': val=self.responseDur
                else: val=self.lattice.getLatticeParameter(n)
                valcmn.append(val)
        if valmode=='all' or valmode=='specific': 
            for n in self._sthdr[0]:
                parVariant, parName = self._describeParm(n)
                val = self.lattice.getDotParameter(parName, parVariant)
                valspec.append(val)
            #print "valspec"
            #print valspec
                
        return valcmn + valspec        




