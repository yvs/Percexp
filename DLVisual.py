
# DLVisual: Set of classes for DotLattice library visual display
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
import numpy as np
from numpy import pi, sqrt, sin, cos, arccos, arctan, tanh
from random import *

import OpenGL.GL as GL
import OpenGL.GLU as GLU
import OpenGL.GLUT as GLUT



class DotLattice(visual.ElementArrayStim):
    """
    Class for the dot lattice stimuli bound by a circular aperture
    """
    def __init__(self,
                 win,
                 session,
                 dotVariantsList=[],
                 tex_res = 128
                 ):
        """
        :Parameters:
        
            win :
                a :class:`~psychopy.visual.Window` object (required)
                                 
            session:
                instance of DLSession class providing the configuration data for lattice

            dotVariantsNameList: 
                contains names of those lattice parameters that are being expected to have
                variants
                
            tex_res : 
                the number of pixels in the dot texture (overridden if an array 
                or image is provided)                       
            
        """
        #
        # DLSession object is the main source of lattice data
        #
        self.ss=session
        session.setLattice(self) #make yourself known to session

        self.bgnd = session.getExpdata(session.winBgnd)
        
        #
        # Dot size - should be known early on
        #
        self.dotSize = session.getExpdata(session.dotSize)
        
        #
        # Quad definition
        #
        self.aDistance = session.getExpdata(session.aDist)
        self.baAspect = session.getExpdata(session.baAspect)
        self.gamma = session.getExpdata(session.gamma)
        self.theta = session.getExpdata(session.theta)

        #
        # Quad basis visibility (for demo mode)
        #
        self.drawBasis=False
        colors = ['red', 'green', 'blue']
        self.basis = [visual.ShapeStim(win=win, units="pix", vertices=((0,0),(1,1)), lineColor=colors[i], lineWidth=3, closeShape=False) 
                         for i in range(3)]

        #
        # Quad computed data
        #
        self.dirAngles = [0,0,0,0]
        self.dirMags = [0,0,0,0]
        self._basisQuad = [[0,0],[0,0],[0,0],[0,0]]
        self._updateBasisQuad()

        #
        # Possible dot variants
        #
        self.dotVariants=[]

        #
        # Dot variants alternation data
        #
        self.grid=[]
        
        #
        # Dot texture
        #
        dotTex = self._createDotTexture(tex_res)

        #
        # Dot mask
        #
        dotMask = self._createDotMask(tex_res)
            
        #
        # Create aperture
        #
        self.apertureSize = session.getExpdata(session.aperSize)
        self.aperturePos = session.getExpdata(session.aperPos)
        self.apertureColor = session.getExpdata(session.aperBgnd)
        self.apertureVisible = session.getExpdata(session.aperShow)
        self._abgnd = None
        self._createApertureBgnd(win)
        self.aperture = Aperture(win=win, size_px=self.apertureSize, pos_px=self.aperturePos)
        self.showAperture(session.getExpdata(session.aperShow))
        self.aperture.disable() #shouldn't be visible initially

        #
        # Create lattice layout
        #
        padding = 3
        self.dotSpacing = {}
        self._updateDotSpacing()
        self.dotsPerSide = {
            'a': int(self.apertureSize / (self.dotSize*2) + padding),
            'b': int(self.apertureSize / (self.dotSize*2) + padding)
        }
        self.numDots = self.dotsPerSide['a'] * self.dotsPerSide['b']
        self.dotXYs = self._calcLattice()
        # call base class init    
        visual.ElementArrayStim.__init__(self, win, units='pix', xys=self.dotXYs, elementTex=dotTex, elementMask=dotMask,
                                  texRes=tex_res, nElements = self.numDots, oris=0, sizes=self.dotSize, sfs=1.0)
        #
        # Dots luminance
        #
        self.dotLum=session.getExpdata(session.dotLum)

        #
        # Dots orientation
        #
        self._symbolicDotAngles=False
        self.dotAngle = session.getExpdata(session.dotAngle)
        self.setDotAngle(self.dotAngle)
        
        #
        # Create dots variants (if needed)
        #
        for variant in dotVariantsList:
            #parlist = session.getExpdataRaw(name)
            #ln = len(parlist)
            #if ln>1:
            if type(variant) in [list, tuple]:
                if len(variant)==2:
                    name = variant[0]
                    for i in range(variant[1]):
                        val = session.getExpdata(name, i+1, True)
                        self.setDotVariantParameters(i+1,{ name: val })
                else: print "dotVariantsList: the length of a variant should be 2"
            else: print "dotVariantsList: variant should be a list or tuple"        

        #
        # Compute dots alternation
        #
        self.dotAlternation=''
        self.setDotAlternation(session.getExpdata(session.dotAlter))

        #
        # Prepare dots permutations data for masking
        #
        d = float(2*session.getExpdata(session.maskRad)) #in pixels now
        self._pnum = 10
        self._plutx = range(self._pnum)
        self._perm = np.random.rand(self._pnum, self.numDots, 2)
        self._perm *= d
        self._perm -= d/2


    """

    LATTICE CONFIGURATION

    """

    def _createDotMask(self, res):
        dotMask=self.ss.getExpdata(self.ss.dotMask)    
        if type(dotMask) == str:    
            s = dotMask.lower()
            if s.find('log') != -1:
                n=1
                if len(s)==4:
                    sn=s[-1]
                    if sn.isdigit(): n=int(sn)
                pwr=1+(n-1)*2
                rad = visual.makeRadialMatrix(res)
                sq = np.sqrt(2)+0.0001
                mx=np.log(sq)**pwr
                scale=2/(1+mx)
                p = np.log(sq-rad)**pwr
                dotMask = (p>=-1)*((p+1)*scale)-1
        return dotMask        

    def _createDotTexture(self, res):    
        tex = self.ss.getExpdata(self.ss.dotTex)
        if type(tex) == str:
            s = tex.lower()
            if s.find('tanh') != -1:
                n=1
                if len(s)==5:
                    sn=s[-1]
                    if sn.isdigit(): n=int(sn)
                tanhA = 0.5*(2-int(n-1)/3)
                tanhC = 4+(0 == int(n)%3)
                tanhX = 5 + (0 != int(n-1)%3) +2*(0 == int(n)%3)
                tex = self.makeTanh(res, tanhA, tanhC, tanhX, -1, 1)
            elif s == 'flat':
                tex = np.ones((res,res), dtype=float)
        return tex
    
    def setPracticeMode(self, mode):
        self.practice=mode

    def basisVisible(self, bv):
        self.drawBasis=bv

    def setDotAlternation(self, mode=None):
        self.dotAlternation=mode
        if mode == 'a':
            self.setDotVariantsGrid([[0,1]])
        elif mode == 'b':    
            self.setDotVariantsGrid([0,1])
        else:
            self.setDotVariantsGrid([])

    def _createApertureBgnd(self, win):
        v = self.apertureSize/2
        p = self.aperturePos
        quad = [[-v,-v],[v,-v],[v,v],[-v,v]]
        self._abgnd = visual.ShapeStim(win=win, units="pix", pos=p, vertices=quad,
                               lineColor=self.apertureColor, fillColor=self.apertureColor)
        
    def setDotVariantParameters(self, variant, parms):
        """
        This method sets parameters for the dot variant.
        'parms' should be a dictionary of parameter names: values.
        """
        if variant>0:
            variant -= 1 #make it to be an index into self.dotVariants
            vlen = len(self.dotVariants)
            if variant <= vlen:
                if variant==vlen: #add new dot variant
                    self.dotVariants.append(parms)
                else: #update existing dot variant
                    vdict = self.dotVariants[variant]
                for k, v in parms.iteritems():
                    if k=="dotAngle":
                        self.setDotAngle(v, variant)
                    elif k=="dotLum":
                        self.setDotLum(v, variant)
                    else: vdict[k]=v
            else: print "addDotVariantParameters(): 'variant' parameter (%d) is larger than dotVariants list length (%d)" % (variant, vlen)
        else: print "addDotVariantParameters(): 'variant' parameter (%d) has to be positive" % variant

    def setDotVariantsGrid(self, grid):
        """
        Grid specifies the pattern of dot variants (if more than one is present in the lattice)
        'grid' parameter defines the pattern in a row-wise fashion. For instance, if two
        dot variants are present, then
        a) to create lattice with alternating columns grid=[[0,1]]
        b) to create lattice with alternating rows grid=[0,1]
        c) for checkerboard pattern grid=[[0,1],[1,0]] etc
        In examples 0 specifies the main dot variant that is initially created with the lattice,
        1 specifies the first dot variant in the dotVariants list.
        """
        self.grid=grid
        self._updateDotAngles()
        self._updateDotLums()

    def _updateBasisQuad(self):
        w = self.aDistance #1.0 #width of quad
        h = w * self.baAspect * sin(np.radians(self.gamma)) #height of quad
        bx = w * self.baAspect * cos(np.radians(self.gamma)) #x coord of point b
        dp = w + bx #x projection of vector d
        cp = np.abs(w - bx) #x projection of vector c
        
        # Compute directions
        self.dirAngles[0]= self.theta
        self.dirAngles[1]= self.theta + self.gamma
        self.dirAngles[2]= self.theta + 180 - np.degrees(arctan(h/cp))
        self.dirAngles[3]= self.theta + np.degrees(arctan(h/dp))

        # Compute vector magnitudes
        self.dirMags[0]= w
        self.dirMags[1]= w * self.baAspect
        self.dirMags[2]= sqrt(cp**2 + h**2)
        self.dirMags[3]= sqrt(dp**2 + h**2)

        # Compute quad's points
        # a direction point
        self._basisQuad[1][0]= self.dirMags[0] * cos(np.radians(self.dirAngles[0]))
        self._basisQuad[1][1]= self.dirMags[0] * sin(np.radians(self.dirAngles[0]))
        # b direction point
        self._basisQuad[2][0]= self.dirMags[1] * cos(np.radians(self.dirAngles[1]))
        self._basisQuad[2][1]= self.dirMags[1] * sin(np.radians(self.dirAngles[1]))
        # d direction point
        self._basisQuad[3][0]= self.dirMags[3] * cos(np.radians(self.dirAngles[3]))
        self._basisQuad[3][1]= self.dirMags[3] * sin(np.radians(self.dirAngles[3]))

        #update coordinates of basis
        self._updateBasis()
        
    def setApertureSize(self, apertureSize):
        self.aperture.setSize(apertureSize)

    def setAperturePos(self, aperturePos):
        self.aperture.setPos(aperturePos)

    def setGamma(self, gamma):
        self.gamma = gamma
        self._updateBasisQuad()
        self._updateXYs()
        if self._symbolicDotAngles==True:
            self._updateDotAngles()
        #print "gamma = %f" % gamma

    def setTheta(self, theta):
        self.theta = theta
        self._updateBasisQuad()
        self._updateXYs()
        if self._symbolicDotAngles==True:
            self._updateDotAngles()
        #print "theta = %f" % theta

    def setDotLum(self, dl, dotVariant=None):
        # set proper dot variant's data
        if dotVariant==None: self.dotLum = dl
        else: self.dotVariants[dotVariant-1]['dotLum']=dl
        self._updateDotLums()

    def _updateDotLums(self):                
        # prepare angles
        a = []
        a.append(self.dotLum)
        for v in self.dotVariants:
            if 'dotLum' in v:
                a.append(v['dotLum'])
        #print a        
        # set dot lum(s)
        g_rows = len(self.grid)
        if g_rows==0 or len(a)==1:
            visual.ElementArrayStim.setRgbs(self,a[0])
        else:
            lums = range(self.numDots)
            rows = self.dotsPerSide['b']
            cols = self.dotsPerSide['a']
            g_row = 0
            for r in range(rows):
                g_item = self.grid[g_row] #grid row description
                if not type(g_item) in [list, tuple]:
                    g_item = [g_item]
                g_cols = len(g_item)
                g_col = 0
                rowed = r*cols
                for c in range(cols):
                    lums[rowed+c] = a[g_item[g_col]]
                    g_col += 1
                    if g_col >= g_cols: g_col=0
                g_row += 1
                if g_row >= g_rows: g_row=0
            visual.ElementArrayStim.setRgbs(self,lums)
            #print lums

    def setDotAngle(self, da, dotVariant=None):
        # set proper dot variant's data
        if dotVariant==None: self.dotAngle = da
        else: self.dotVariants[dotVariant]['dotAngle']=da
        self._updateDotAngles()

    def _updateDotAngles(self):                
        # prepare angles
        a = []
        a.append(self._dereferenceDotAngle(self.dotAngle))
        for v in self.dotVariants:
            if 'dotAngle' in v:
                a.append(self._dereferenceDotAngle(v['dotAngle']))
        #print a        
        # set dot angle(s)
        g_rows = len(self.grid)
        if g_rows==0 or len(a)==1:
            visual.ElementArrayStim.setOris(self,a[0])
        else:
            #oris = np.ndarray(shape=(self.numDots,1), dtype=float)
            oris = range(self.numDots)
            rows = self.dotsPerSide['b']
            cols = self.dotsPerSide['a']
            g_row = 0
            for r in range(rows):
                g_item = self.grid[g_row] #grid row description
                if not type(g_item) in [list, tuple]:
                    g_item = [g_item]
                g_cols = len(g_item)
                g_col = 0
                rowed = r*cols
                for c in range(cols):
                    oris[rowed+c] = a[g_item[g_col]]
                    g_col += 1
                    if g_col >= g_cols: g_col=0
                g_row += 1
                if g_row >= g_rows: g_row=0
            visual.ElementArrayStim.setOris(self,oris)
            #print "calced oris: got there!"
            #print oris

    def _dereferenceDotAngle(self, da):
        if type(da)==str:
            if da=="a": a=self.dirAngles[0]
            elif da=="b": a=self.dirAngles[1]
            elif da=="c": a=self.dirAngles[2]
            elif da=="d": a=self.dirAngles[3]
            else: a = da
            self._symbolicDotAngles=True
        else: a=da    
        return a    

    def setDotSize(self, size):
        self.dotSize = size
        #self._updateDotSpacing()
        visual.ElementArrayStim.setSizes(self,size)
        #self._updateXYs()
        #print "dotSize = %f" % size

    def setADistance(self, aDist):
        self.aDistance = aDist
        self._updateBasisQuad()
        self._updateDotSpacing()
        self._updateXYs()
        #print "aDistance = %f" % aDist

    def setBAAspect(self, baAspect):
        self.baAspect = baAspect
        self._updateBasisQuad()
        self._updateDotSpacing()
        self._updateXYs()
        if self._symbolicDotAngles==True:
            self._updateDotAngles()
        #print "baAspect = %f" % baAspect

    def setContrast(self, contr):
        visual.ElementArrayStim.setContrs(self,contr)
        #print "contrast = %f" % contr

    def _updateXYs(self):    
        self.dotXYs = self._calcLattice()
        visual.ElementArrayStim.setXYs(self,self.dotXYs)

    def _updateDotSpacing(self):    
        self.dotSpacing['a'] = self.aDistance #*self.dotSize
        self.dotSpacing['b'] = self.aDistance*self.baAspect #*self.dotSize

    def _updateBasis(self):
        vertices = [[0,0],[0,0]]
        for i in range(2):
            vertices[1]=self._basisQuad[i+1]
            self.basis[i].setVertices(vertices)
        vertices.append(self._basisQuad[1])
        vertices[0]=vertices[1]
        vertices[1]=self._basisQuad[3]
        self.basis[2].setVertices(vertices)

    def _calcLattice(self):
        rows = self.dotsPerSide['b']
        cols = self.dotsPerSide['a']
        halfXspan = cols/2 - 0.5*(1-cols%2)
        #print halfXspan
        halfYspan = rows/2 - 0.5*(1-rows%2)
        #print halfYspan
        ys, xs = np.mgrid[-halfYspan:halfYspan:1j*rows, -halfXspan:halfXspan:1j*cols]
        xs = xs.reshape(self.numDots,1)
        ys = ys.reshape(self.numDots,1)
        #print xs
        #print ys
        sinG = sin(np.radians(self.gamma))
        cosG = cos(np.radians(self.gamma))
        sinT = sin(np.radians(self.theta))
        cosT = cos(np.radians(self.theta))
        aDist = self.dotSpacing['a']
        bDist = self.dotSpacing['b']
        dotXYs = np.ndarray(shape=(self.numDots,2))
        for i in range(self.numDots):
            x = xs[i]*aDist + ys[i]*bDist*cosG #x coord with gamma angle applied
            y = ys[i]*bDist*sinG # same for y
            # theta rotation
            dotXYs[i][0] = x * cosT - y * sinT 
            dotXYs[i][1] = x * sinT + y * cosT
        return dotXYs    

    def makeTanh(self, res, a, c, x, val_min, val_max):
        span = val_max - val_min
        x_dist, y_dist = np.mgrid[-x:x:1j*res, 0:res]
        tex = (1 + a * cos(pi*x_dist/2)) * (tanh(x_dist+c) - tanh(x_dist-c))
        tex /= (np.max(tex) / span) #normalize and scale
        tex += val_min #shift to start with val_min instead of 0
        return tex



    """

    LATTICE DATA ACCESS 

    """

    def getBgndColor(self):
        return self.bgnd

    def pointInsideAperture(self, pt):
        s = self.apertureSize/2
        return all((pt[0]>-s, pt[0]<s, pt[1]>-s, pt[1]<s))

    def getDotParameter(self, parm, dotVariant=None):
        if parm=='dotAlternation':
            if self.dotAlternation==None: da='no'
            else: da=self.dotAlternation
            return da
        else:
            if dotVariant==None: dv='None'
            else: dv = dotVariant
            s = "getDotParameter(): wrong parameter (%s) variant: %s" % (parm, str(dv))
            prim = (dotVariant==None or dotVariant==0)
            if prim:
                if parm=='dotSize': return self.dotSize
                elif parm=='dotLum': return self.dotLum
                elif parm=='dotAngle': return self.dotAngle
            else:
                try: return self.dotVariants[dotVariant-1][parm]
                except: print s
        
    def getLatticeParameter(self, parm):
        if parm=='theta': return self.theta
        elif parm=='baAspect': return self.baAspect
        elif parm=='gamma': return self.gamma
        elif parm=='lenA': return self.dirMags[0]
        elif parm=='angleA': return self.dirAngles[0]
        elif parm=='lenB': return self.dirMags[1]
        elif parm=='angleB': return self.dirAngles[1]
        elif parm=='lenC': return self.dirMags[2]
        elif parm=='angleC': return self.dirAngles[2]
        elif parm=='lenD': return self.dirMags[3]
        elif parm=='angleD': return self.dirAngles[3]
        else: print "getlatticeParameter(): unknown parameter name (%s)" % parm
        


    """

    LATTICE DRAWING 

    """

    def draw(self):
        if self.apertureVisible==True:
            self.aperture.enable()
            self._abgnd.draw()
        else: self.aperture.disable()
        
        visual.ElementArrayStim.draw(self)
        if self.drawBasis == True:
            for b in self.basis:
                b.draw()

        if self.apertureVisible==True:
            self.aperture.disable()

    def drawMask(self, repDur):
        clock = core.Clock()
        count=0
        #while count<reps:
        while repDur>=clock.getTime():
            if not count%self._pnum:
                shuffle(self._plutx)
                i=0
            # permute dot positions
            ndx = self._plutx[i]
            perm = self._perm[ndx]
            XYs = self.dotXYs + perm
            visual.ElementArrayStim.setXYs(self, XYs)
            self.draw()
            self.win.flip()
            #core.wait(repDur)
            i+=1
            count+=1    
        #restore unpermuted positions    
        visual.ElementArrayStim.setXYs(self, self.dotXYs)
            

    def showAperture(self, show):
        self.apertureVisible = show




class LatticeController():        
    """
    Class to generate lattice parameters in a controlled fashion
    Three types of controllers can be created:
    A) Control. This one is for an individual lattice parameter
    B) Matrix. This is a cross-combination of two or more Controls.
               Each Control defines one dimension of a matrix.
    C) Block. This represents an explicit sequence of multi-parameter blocks
              to control each subsequent experiment trial
    """
    def __init__(self,
                 lattice
                 ):
        """
        :Parameters:
        
            lattice :
                instance of dot lattice to control
        """
        self.lat=lattice
        self.ctrl={}
        self.mat={}
        self.blockseqs={}



    """

    CONFIGURATION 

    """

    def resetControl(self):
        self.ctrl={}

    def addBlocksToSequence(self, seqName, blocks):
        if seqName in self.blockseqs:
            seq = self.blockseqs[seqName]
            lblk = len(seq['parnames'])
            lins = len(blocks)
            #check whether parmValues contains more than one block
            count=0
            if type(blocks[0]) in [tuple, list]:
                if lblk == len(blocks[0]): #it seems like there is more than one block
                     count=lins
            if count==0 and lins == lblk:
                blocks=[blocks]
                count=1
            if count==0:
                print "addBlockToSequence(): check parameters!"
            else:    
                vals = seq['blocks']
                for i in range(count):
                    vals.append(blocks[i])
                #print seq
        else:
            print "addBlockToSequence(): unknown sequence name: %s!" % seqName
            

    def createBlockSequence(self, seqName, parmNames):
        """
        Blocks holds combinations of given parameters' values.
        Unlike in matrix the combinations are not exhaustive,
        they are not derived from the controlled parameters dictionary,
        but are explicitely provided by user.
        'parmNames' should be a tuple or list of parameter names, e.g. ('dotSize','dotLum2',gamma')
        """

        if type(parmNames) in [list, tuple]:
            count=len(parmNames)
        elif type(parmNames)==str:
            parmNames=[parmNames]
            count=1
        else:
            print "createBlockSequence(): parmNames should be specified as tuple, list or a string"
            return
        if seqName in self.blockseqs:
            print "createBlockSequence(): sequence %s already exists!" % seqName
            return
        else:
            seq = {
                'parnames': parmNames,
                'blocks': []
                }
            self.blockseqs[seqName]=seq
            #print self.blockseqs
            
        
    def finalizeSequence(self, seqName, numCopies=1, shuffled=True):
        if seqName in self.blockseqs:
            seq = self.blockseqs[seqName]
            blkcnt = len(seq['blocks'])
            seqlen = blkcnt*numCopies
            #set lookup order                
            lut=range(seqlen)
            if numCopies>1:
                for i in range(seqlen):
                    lut[i] %= blkcnt
            if shuffled==True: shuffle(lut)
            seq['order']=lut
            seq['current']=0
        else:
            print "addBlockToSequence(): unknown sequence name: %s!" % seqName
        


    def setControl(self,
                   parName,
                   parms
                   ):
        """
        Sets controlling context for specific lattice parameter
        """
        #create controller item
        good=True
        #set defaults
        parMin=None
        parMax=None
        parStep=None
        stepCount=None
        parTable=None
        shuffled=True
        method='regular'
        #update from dictionary
        if 'parMin' in parms: parMin = parms['parMin']
        if 'parMax' in parms: parMax = parms['parMax']
        if 'parStep' in parms: parStep = parms['parStep']
        if 'stepCount' in parms: stepCount = parms['stepCount']
        if 'parTable' in parms: parTable = parms['parTable']
        if 'shuffled' in parms: shuffled = parms['shuffled']
        if 'method' in parms: method = parms['method']
        #set ranges based on known limitations
        if parTable==None:
            if parName=="gamma":
                if parMin==None: parMin=60.0
                if parMax==None: parMax=90.0
            elif parName=="theta":
                if parMin==None: parMin=0.0
                if parMax==None: parMax=180.0
            elif parName=="dotAngle":
                if parMin==None: parMin=0.0
                if parMax==None: parMax=180.0
            elif parName=="baAspect":
                if parMin==None: parMin=1.0
            elif parName=="contrast":
                if parMax==None: parMax=1.0

            if method=='random':
                if all((stepCount!=None, parMin!=None, parMax!=None)):
                    r = abs(parMax-parMin)
                    parTable=np.random.rand(stepCount)
                    parTable *= r
                    parTable += parMin
                else:
                    print "setControl(): not enough data (count, min, max) for random generation!"
                    good=False
            elif parStep!=None or stepCount!=None:
                if parMin!=None and parMax!=None:
                    if parStep!=None: parTable=np.arange(parMin,parMax+parStep,parStep)
                    else: parTable = np.linspace(parMin,parMax,stepCount,True)
                else:
                    print "setControl(): not enough data (min, max) for random generation!"
                    good=False
            else:
                print "setControl(): not enough data (step, count) for random generation!"
                good=False

        if good==True:
            count = len(parTable)
            if count==0: print "setControl(): zero length of parameter table!"
            else:
                lut=range(count)
                if shuffled==True:
                    shuffle(lut)
                controlUnit = {
                    'order': lut,
                    'table': parTable,
                    'current': 0,
                    'shuffled': shuffled
                    }
                self.ctrl[parName]=controlUnit
                #print self.ctrl

    def setControlMatrix(self,
                         matName,
                         parms):
        """
        Creates exhaustive combination of given parameters' values
        """
        if type(parms) in [list, tuple]:
            #number of matricized parameters is expected to be 2 or 3
            count = len(parms)
            if not all((count>=2, count<=3)):
                print "setControlMatrix: number of parms should be 2 or 3"
                return
            permcount=1
            names=[]
            shape=[]
            for p in parms:
                #make sure all parameters are known
                try:
                    c = self.ctrl[p]
                    l = len(c['table'])
                    names.append(p)
                    shape.append(l)
                    permcount *= l
                except:
                    print"setControlMatrix(): unknown parameter name %s" % p
                    return
            ndxar = np.ndarray(shape=(permcount,count), dtype=int)
            if count==2:
                istep = shape[1]
                jstep=1
            else:    
                istep = shape[1]*shape[2]
                jstep = shape[2]
            for i in range(shape[0]):
                i_ndx = i*istep
                for j in range(shape[1]):
                    j_ndx = j*jstep
                    if count == 2:
                        ndx=i_ndx+j_ndx
                        ndxar[ndx][0]=i
                        ndxar[ndx][1]=j
                    else:
                        for k in range(shape[2]):
                            ndx=i_ndx+j_ndx+k
                            ndxar[ndx][0]=i
                            ndxar[ndx][1]=j
                            ndxar[ndx][2]=k
            #set lookup order                
            lut=range(permcount)
            shuffle(lut)
            #add to the matrix dictionary
            m = {
                'parnames': names,
                'indices':  ndxar,
                'order':    lut,
                'current':  0
                }
            self.mat[matName]=m
        else:
            print "setControlMatrix(): parms should be specified as tuple or list"
                         

    """

    DATA ACCESS 

    """

    def getSequenceBlockCount(self, seqName):
        if seqName in self.blockseqs:
            seq = self.blockseqs[seqName]
            return len(seq['blocks'])
        else: 
            print "getSequenceBlockCount(): unknown sequence name: %s!" % seqName
            
        
    def getPermutationsCount(self):
        count=1
        for k, v in self.ctrl.iteritems():
            count *= len(v['table'])
        return count

    def getMatrixPermutationsCount(self, matName):
        """
        Returns number of permutations this matrix contains
        """
        try:
            m = self.mat[matName]
            return len(m['order'])
        except:
            print "getMatrixPermutationsCount(): unknown matrix name '%s'" % matName
        

    def getSequenceLength(self, seqName):
        if seqName in self.blockseqs:
            seq = self.blockseqs[seqName]
            return len(seq['order'])
        else: 
            print "getSequenceLength(): unknown sequence name: %s!" % seqName

    def getLongestSequenceLength(self):
        ln = 0
        for k in self.blockseqs.iterkeys():
            seql = self.getSequenceLength(k)
            ln = max(ln,seql)
        return ln    
        

    """

    CONTROL 

    """

    def rewindAll(self):
        #rewind blocks
        for b in self.blockseqs.itervalues():
            b['current']=0
        #rewind matrices
        for m in self.mat.itervalues():
            m['current']=0
        #rewind individual controls
        self.rewindParms()
        
        
    def rewindParms(self, parms=None):
        if parms==None:
            for v in self.ctrl.itervalues():
                v['current']=0
        else:
            t = type(parms)
            if t==str:
                try:
                    c = self.ctrl[parms]
                    c['current']=0
                except:
                    print "rewind(): non-existent parameter %s!" % parms
            elif t in [list, tuple]:   
                for n in parms:
                    try:
                        c = self.ctrl[n]
                        c['current']=0
                    except:
                        print "rewind(): non-existent parameter %s!" % n


    def nextInControl(self, parNames, ndx='auto'):
        """
        Sets lattice parameter(s) to the currently expected value(s) and advances the control pointer
        """
        if type(parNames) == str:
            self._nextInControl(parNames,ndx)
        elif type(parNames) in [list, tuple]:
            for n in parNames:
                self._nextInControl(n,ndx)

    def _nextInControl(self, parName, index):
        """
        Sets lattice parameter(s) to the currently expected value(s) and advances the control pointer
        """
        try:
            c = self.ctrl[parName]

            if index=='auto': ndx = c['current']
            else: ndx=index
            
            isLast = ((len(c['table'])-ndx) == 1)
            i = c['order'][ndx]
            v = c['table'][i]
            self._setLatticeParm(parName,v)
            #update control pointer    
            if isLast==True:
                if c['shuffled']==True: shuffle(c['order'])
                c['current']=0
            else: c['current']=ndx+1
        except:
            print "_nextInControl(): unknown parameter '%s'" % parName

    def _setLatticeParm(self, parName, parVal):
        if parName=="gamma":
            self.lat.setGamma(parVal)
        elif parName=="theta":
            self.lat.setTheta(parVal)
        elif parName=="dotAngle":
            self.lat.setDotAngle(parVal)
        elif parName=="dotSize":
            self.lat.setDotSize(parVal)
        elif parName=="aDistance":
            self.lat.setADistance(parVal)
        elif parName=="baAspect":
            self.lat.setBAAspect(parVal)
        elif parName=="contrast":
            self.lat.setContrast(parVal)
        elif parName.find("dotLum") != -1:
            sn = parName[-1]
            n = None
            if sn.isdigit(): n = int(sn)-1
            #print "lum=%d" % n
            self.lat.setDotLum(parVal,n)
        elif parName=="dotAlternation":
            self.lat.setDotAlternation(parVal)
        
    def nextInMatrix(self, matName):
        """
        Sets lattice parameter(s) to the combination of controled parameters stored in matrix
        """
        if type(matName) == str:
            try:
                m = self.mat[matName]
                n = m['parnames']
                #print n
                ind = m['current']
                isLast = ((len(m['order'])-ind) == 1)
                ndx = m['order'][ind]
                for i in range(len(n)):
                    self.nextInControl(n[i],m['indices'][ndx][i])
                #update control pointer    
                if isLast==True: m['current']=0
                else: m['current'] = ind+1
            except:
                print "nextInMatrix(): unknown matrix name '%s'" % matName
        else: 
            print "nextInMatrix(): given matrix name parameter is not a string!"


    def nextInSequence(self, seqName):
        if seqName in self.blockseqs:
            seq = self.blockseqs[seqName]
            parms = seq['parnames']
            vals = seq['blocks'][seq['order'][seq['current']]]
            for i in range(len(parms)):
                self._setLatticeParm(parms[i],vals[i])
            #advance pointer
            seq['current'] += 1
            if seq['current']==len(seq['order']): #wrap around
                seq['current']=0
        else: 
            print "nextInSequence(): unknown sequence name: %s!" % seqName
        



class Aperture:
    """Used to create a shape (usually circular) to restrict a stimuli
    visibility area
    """
    def __init__(self, win, size_px, pos_px=(0,0)):
        self.win=win
        self._quad=GLU.gluNewQuadric() 
        #self.winSize = win.monitor.getSizePix()
        #print self.winSize
        if win.winType == 'pyglet':
            self.winSize = win.winHandle.get_size()
        self.winAspect = self.winSize[1] / float(self.winSize[0])
        self.setSize(size_px, False)
        self.setPos(pos_px)

    def reset(self):
        self.enable()
        GL.glClearStencil(0)
        GL.glClear(GL.GL_STENCIL_BUFFER_BIT)
        
        GL.glPushMatrix()
        GL.glTranslate(self.pos[0], self.pos[1], 0)

        GL.glDisable(GL.GL_LIGHTING)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glDepthMask(GL.GL_FALSE)
        
        GL.glStencilFunc(GL.GL_NEVER, 0, 0)
        GL.glStencilOp(GL.GL_INCR, GL.GL_INCR, GL.GL_INCR)
        GL.glColor3f(1.0,1.0,1.0)
        GL.glScalef(self.winAspect, 1.0, 1.0)
        GLU.gluDisk(self._quad, 0, self.size, 120, 2)
        GL.glStencilFunc(GL.GL_EQUAL, 1, 1)
        GL.glStencilOp(GL.GL_KEEP, GL.GL_KEEP, GL.GL_KEEP)

        GL.glPopMatrix()

    def setSize(self, size_px, needReset=True):
        self.size = size_px / float(self.winSize[1])
        if needReset: self.reset()

    def setPos(self, pos_px, needReset=True):
        self.pos = (pos_px[0]/float(self.winSize[0]), pos_px[1]/float(self.winSize[1]))
        #print pos_px
        #print self.pos
        if needReset: self.reset()

    def enable(self):
        GL.glEnable(GL.GL_STENCIL_TEST)

    def disable(self):
        GL.glDisable(GL.GL_STENCIL_TEST)




class FixationScreen():
    """Class for presenting pre-stimulus fixation
    """
    def __init__(self,
                 win,
                 session,
                 fixType="cross",
                 pos=(0.0, 0.0),
                 size=20,
                 lineWidth=1.0,
                 lineColor='white'
                 ):
        """
        :Parameters:
        
            win :
                a :class:`~psychopy.visual.Window` object (required)
                                 
            fixType :
                type of fixation shape
                                 
            size : 
                sizes of fixation shape

        """
        self.ss = session
        session.setFixation(self)
        self.dur = session.getExpdata(session.fixDur)
        
        if type(size) in [tuple,list]:
            self.size = np.array(size,float)
        else:
            self.size = np.array((size,size),float)

        self.shapes=[]

        if fixType=="cross":
            
            xy = [[-self.size[0]/2,0],[self.size[0]/2,0]]
            h = visual.ShapeStim(win=win, units="pix", vertices=xy, lineWidth=lineWidth, lineColor=lineColor,
                                  closeShape=False, pos=pos)
            self.shapes.append(h)
            
            xy = [[0,-self.size[1]/2],[0,self.size[1]/2]]
            v = visual.ShapeStim(win=win, units="pix", vertices=xy, lineWidth=lineWidth, lineColor=lineColor,
                                  closeShape=False, pos=pos)
            self.shapes.append(v)

    def draw(self):
        for shape in self.shapes:
            shape.draw()



        
class _ResponseVariant(visual.ShapeStim):
    """Class for presenting specific response variant
    """
    def __init__(self,
                 win,
                 vertices,
                 name=None,
                 pos=(0.0, 0.0),
                 size=128,
                 ori=0.0,
                 lineWidth=3.0,
                 lineColor='white',
                 fillColor='grey'
                 ):
        """
        :Parameters:
        
            win :
                a :class:`~psychopy.visual.Window` object (required)
                                 
            name : 
                name of an option for easier identification

            pos:
                position of the center
                
            size:
                diameter
                
            ori:
                angular orientation
                
        """
        #init base object
        visual.ShapeStim.__init__(self,win=win, units='pix', vertices=vertices,
                                  lineWidth=lineWidth, lineColor=lineColor, fillColor=fillColor,
                                  closeShape=True, ori=ori, pos=pos)
        self.name = name
        self.varPos = pos
        self.size = size
        self._setBoundBox()

    def setVariantPos(self, pos):
        visual.ShapeStim.setPos(self, pos)
        self.varPos = pos
        self._setBoundBox()

    def _setBoundBox(self):
        radius=self.size/2
        x0, y0 = self.pos
        self.boundBox = ((x0-radius,y0-radius),(x0+radius,y0+radius))

    def setName(self, name):
        self.name=name

    def isPointWithin(self, point):
        bmin = self.boundBox[0]
        bmax = self.boundBox[1]
        if all((point[0]>=bmin[0], point[1]>=bmin[1], point[0]<=bmax[0], point[1]<=bmax[1])):
            return True
        else: return False

    

class ResponseScreen():
    """Class for presenting response variants and collecting responses from user
    """
    def __init__(self,
                 win,
                 session,
                 lineColor = 'black',
                 fillColor = -0.5,
                 lineWidth = 3
                 ):
        """
        :Parameters:
        
            win :
                a :class:`~psychopy.visual.Window` object (required)
                                 
            lat :
                a :class:`DotLatticeStim` object (required)
                                 
            variantGap : 
                horizontal and vertical gaps between response items centers

            variantSize:
                size of the response item 

        """
        self.win=win
        self.ss=session
        session.setResponseObject(self)
        self.lattice=session.getLattice()
        self.variantSize=session.getExpdata(session.respSize)
        self.varNames = ('a', 'b', 'c', 'd')
        self._lutx = range(4)
        self.selected = -1
        self.prevsel = self.selected
        self.winBgnd = session.getExpdata(session.respBgnd)

        variantGap =session.getExpdata(session.respDist)
        if type(variantGap) in [int, float]:
           self.variantGap= (variantGap, variantGap) 
        else: self.variantGap=variantGap
        
        #position 4 response variants
        x, y = np.mgrid[-self.variantGap[0]/2:self.variantGap[0]/2:1j*2, -self.variantGap[1]/2:self.variantGap[1]/2:1j*2]
        x = x.reshape(4,1)
        y = y.reshape(4,1)
        self.XY = np.ndarray(shape=(4,2))
        for i in range(4):
            self.XY[i][0]=x[i]
            self.XY[i][1]=y[i]
            
        edges=64
        radius=self.variantSize/2
        vert = np.ndarray(shape=(edges+2,2))
        #x0, y0 = pos
        #circle vertices coords
        for i in range(edges):
            angle = -2.0*pi*i/edges
            vert[i][0]= radius * cos(angle)
            vert[i][1]= radius * sin(angle)
        #ori line coords
        vert[edges][0]= radius
        vert[edges][1]= 0
        vert[edges+1][0]= -radius
        vert[edges+1][1]= 0
        #calc oris and shuffle the order   
        self.variants = [_ResponseVariant(win=win, vertices=vert, name=self.varNames[i], size=self.variantSize, 
                         lineColor=lineColor, fillColor=fillColor, lineWidth=lineWidth) 
                         for i in range(4)]
        vert2 = np.resize(vert,(edges,2))
        vert2 *= 0.85
        self.selvariant = visual.ShapeStim(win=win, units='pix', vertices=vert2,
                                  lineWidth=lineWidth, lineColor=fillColor+0.5)
        #print self.variants
        self.updateVariants()


    def getCurVariantName(self):
        return self.varNames[self.selected]

    def getVariantsOrder(self):
        return self.varNames[self._lutx[0]]+self.varNames[self._lutx[1]]+self.varNames[self._lutx[2]]+self.varNames[self._lutx[3]]

    def updateVariants(self, updateObjects=True):
        
        shuffle(self._lutx)
        #print self.order

        if updateObjects == True:
            for i in range(4):
                ndx= self._lutx[i]
                self.variants[ndx].setName(self.varNames[ndx])
                self.variants[ndx].setVariantPos(self.XY[i])
                self.variants[ndx].setOri(-self.lattice.dirAngles[ndx])
                #self.variants[ndx].setLineColor(self.colors[ndx])

    def blinkSelected(self):
        for i in range(6):
            self.draw(i%2)
            self.win.flip()
            core.wait(0.05)
        self.prevsel=-2 #trick to avoid option being selected automatically (TO_DO: remove this when mouse cursor becomes positionable)    

    def draw(self, drawSel=True):
        for v in self.variants:
            v.draw()
        if drawSel and self.selected != -1:
            self.selvariant.draw()

    def variantForPoint(self, point):
        sel = -1
        for i in range(4):
            ndx= self._lutx[i]
            if self.variants[ndx].isPointWithin(point) == True:
                sel = ndx
        if sel != -1 and self.prevsel==-2:
            sel = -1
        elif sel != self.prevsel:
            self.prevsel = sel
            if sel != -1:
                self.selvariant.setPos(self.variants[sel].pos)
        self.selected=sel        
        return sel    



    
        

