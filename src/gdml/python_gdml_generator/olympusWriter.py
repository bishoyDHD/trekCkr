# Defines.cookerWriter class, which extends the writer class that comes packaged with ROOT.  These classes are designed to be used in python scripts for
# generating dense GDML documents in a user-readable way.
#
# To import the module "writer", you must add $ROOTSYS/geom/gdml to your PYTHONPATH environment variable.  Verify that the file writer.py is in that folder.
#
#
#
# VERSION 1.0
# March 13, 2012
# C. O'Connor
#
# Written on a Mac running Python 2.6.1 and GCC 4.2.1 on Darwin
###########################################################################################################################################################

import writer

class.cookerWriter(writer.writer):

    # Initialization function

    def __init__(self, fname):
        """Override the writer constructor with an identical version so that we have an.cookerWriter constructor."""

        self.gdmlfile = fname
        self.define = ['define',{},[]]
        self.materials = ['materials',{},[]]
        self.solids = ['solids',{},[]]
        self.structure = ['structure',{},[]]
        self.document = ['gdml',{'xmlns:gdml':"http://cern.ch/2001/Schemas/GDML",
                                 'xmlns:xsi':"http://www.w3.org/2001/XMLSchema-instance",
                                 'xsi:noNamespaceSchemaLocation':"schema/gdml.cooker.xsd"},
                         [self.define, self.materials, self.solids, self.structure]]
        self.ignoreSolidList = []
        self.ignoreVolumeList = []
        self.worldname = ''

    # Ways to copy one.cookerWriter to make another that is identical but doesn't point to the same memory.
    # cloneList() should be called internally by copy(), not directly by the user.

    def cloneList(self, originalList, cloneList):
       """Turn the list copyList into a clone of originalList and return nothing.  Recursively call the function on any elements that are lists themselves.  This function ensures that copyList is not originalList.  If you wanted that, you could've just written copyList = originalList instead of calling this function."""
        
       for i in range(len(originalList)):                    # Loop over all members of originalList.
           if type(originalList[i])==list:                   # 1. For list elements:
               if i >= len(cloneList):                       #
                   cloneList.append([])                      #    Append if we're past the end of copyList...
               else:                                         #
                   cloneList[i] = []                         #    ...or assign if we're not.
               self.cloneList(originalList[i], cloneList[i]) #    Since this element was a list, call copyList recursively.
           else:                                             # 2. For other elements:
               if i >= len(cloneList):                       #
                   cloneList.append(originalList[i])         #    Append if we're past the end of copyList...
               else:                                         #
                   cloneList[i] = originalList[i]            #    ...or assign if we're not.

    def copy(self):
        """Produce a copy of the.cookerWriter who called the function and return the copy."""

        copy =.cookerWriter(self.gdmlfile)              # Give copy the same name as what it's copying,
        copy.worldname = self.worldname                  # and all of the same attributes.
                                                         #
        self.cloneList(self.define, copy.define)         # 
        self.cloneList(self.materials, copy.materials)   # Note that copy.document does not need to be set,
        self.cloneList(self.solids, copy.solids)         # since it is already defined appropriately just
        self.cloneList(self.structure, copy.structure)   # from the constructor call.
                                                         #
        copy.ignoreSolidList = self.ignoreSolidList[:]   # Directly copy the names listed in the ignore lists.
        copy.ignoreVolumeList = self.ignoreVolumeList[:] #
                                                         #
        return copy                                      # Return copy.

    # Why doesn't writer.py have any set functions at all?

    def setFile(self, fname):
        self.gdmlfile = fname

    def setWorldName(self, name):
        self.worldname = name

    # Functions for adding things that writer objects don't know how to add, or that are formatted differently

    def addElcone(self, name, dx, dy, zcut, zmax):
        self.solids[2].append(['elcone',{'name':name, 'dx':dx, 'dy':dy, 'zcut':zcut, 'zmax':zmax, 'lunit':'cm'}, []])

    def addIsotope(self, name, n, z, a):
        self.materials[2].append(['isotope', {'name':name, 'N':n, 'Z':z}, [['atom', {'unit':'g/mole', 'value':a}, []]]])

    def addElementFromIsotopes(self, name, isotopes):
        subel = []
        for iso in isotopes.keys():
            subel.append(['fraction', {'n':isotopes[iso], 'ref':iso}, []])

        self.materials[2].append(['element', {'name':name}, subel])

    def addMaterialWithAttributes(self, name, state, attribs, elems):
        subel = []
        for att in attribs.keys():
            subel.append([att, {'unit':attribs[att][0], 'value':attribs[att][1]}, []])
            
        for el in elems.keys():
            subel.append(['fraction', {'n':elems[el], 'ref':el}, []])

        self.materials[2].append(['material', {'name':name, 'state':state}, subel])

    def addMaterialG4Al(self, name='G4_Al', state='solid', a=26.9815, z=13, rho=2.699, mee=166):
        self.materials[2].append(['material', {'name':name, 'Z':z, 'state':state},
                                  [['D', {'unit':'g/cm3', 'value':rho}, []],
                                   ['atom', {'unit':'g/mol', 'value':a}, []],
                                   ['MEE', {'unit':'eV', 'value':mee}, []]]])

    def addVolumeWithNamedDaughters(self, volume, solid, material, color, sensdet, copyno, daughters):
        subels = [['materialref',{'ref':material},[]],
                  ['solidref',{'ref':solid},[]]]
        subels.append(['auxiliary',{'auxtype':'Color', 'auxvalue':color},[]])
        if sensdet!='':
            subels.append(['auxiliary',{'auxtype':'SensDet', 'auxvalue':sensdet},[]])
            subels.append(['auxiliary',{'auxtype':'CopyNo',  'auxvalue':copyno}, []])
        for child in daughters:
            subsubels = [['volumeref',{'ref':child[1]},[]],
                         ['positionref',{'ref':child[2]},[]]]
            if child[2]!='':
                subsubels.append( ['rotationref',{'ref':child[3]},[]] )

            subels.append( ['physvol',{'name':child[0]}, subsubels])

        used = 0
        self.structure[2].append(['volume',{'name':volume}, subels, used])

    # Shouldn't be necessary, but it's nice to have a removal handle for solids and volumes

    def removeSolid(self, name):
        while self.hasSolid(name):
            index = 0
            for solid in self.solids[2]:
                if name==solid[1]['name']:
                    index = self.solids[2].index(solid)
            self.solids[2].remove(self.solids[2][index])

    def removeVolume(self, name):
        while self.hasVolume(name):
            index = 0
            for volume in self.structure[2]:
                if name==volume[1]['name']:
                    index = self.structure[2].index(volume)
            self.structure[2].remove(self.structure[2][index])


    # Check functions to see if something has already been added to the.cookerWriter object or not

    def hasPosition(self, name):
        for pos in self.define[2]:
            if pos[0]=='position':
                if name==pos[1]['name']: return True
        return False

    def hasRotation(self, name):
        for rot in self.define[2]:
            if rot[0]=='rotation':
                if name==rot[1]['name']: return True
        return False
    
    def hasSolid(self, name):
        for solid in self.solids[2]:
            if name==solid[1]['name']: return True
        return False

    def hasVolume(self, name):
        for volume in self.structure[2]:
            if name==volume[1]['name']: return True
        return False

    def hasVolumeWithDaughterVolume(self, name, daughtervolume):
        if not self.hasVolume(name): return False
        index = 0
        for volume in self.structure[2]:
            if name==volume[1]['name']:
                index = self.structure[2].index(volume)
        for daughter in self.structure[2][index][2]:
            if daughter[0]=='physvol':
                if daughter[2][0][1]['ref']==daughtervolume: return True
        return False

    # "Ignore" functions for letting some.cookerWriter objects skip over parts of the geometry that they don't want

    def ignoreSolid(self, name):
        if name not in self.ignoreSolidList: self.ignoreSolidList.append(name)

    def ignoreSolids(self, solids):
        for name in solids:
            self.ignoreSolid(name)

    def stopIgnoringSolid(self, name):
        if name in self.ignoreSolidList: self.ignoreSolidList.remove(name)

    def ignoreNoSolid(self):
        self.ignoreSolidList = []

    def isIgnoringSolid(self, name):
        if name in self.ignoreSolidList: return True
        else: return False

    def ignoreVolume(self, name):
        if name not in self.ignoreVolumeList: self.ignoreVolumeList.append(name)

    def ignoreVolumes(self, volumes):
        for name in volumes:
            self.ignoreVolume(name)

    def stopIgnoringVolume(self, name):
        if name in self.ignoreVolumeList: self.ignoreVolumeList.remove(name)

    def ignoreNoVolume(self):
        self.ignoreVolumeList = []

    def isIgnoringVolume(self, name):
        if name in self.ignoreVolumeList: return True
        else: return False

    #  Modified functions for using millimeters, like a sensible person, rather than writer.py's forced centimeters, because it's part of ROOT 

    def addPositionInMM(self, name, x, y, z):
        self.define[2].append(['position',{'name':name, 'x':x, 'y':y, 'z':z, 'unit':'mm'},[]])

    def addBoxInMM(self, name, dx, dy, dz):
        self.solids[2].append(['box',{'name':name, 'x':dx, 'y':dy, 'z':dz, 'lunit':'mm'},[]])
	
    def addParaboloidInMM(self, name, rlo, rhi, dz):
        self.solids[2].append(['paraboloid',{'name':name, 'rlo':rlo, 'rhi':rhi, 'dz':dz, 'lunit':'mm'},[]])
	
    def addArb8InMM(self, name, v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, v5x, v5y, v6x, v6y, v7x, v7y, v8x, v8y, dz):
        self.solids[2].append(['arb8',{'name':name, 'v1x':v1x, 'v1y':v1y, 'v2x':v2x, 'v2y':v2y, 'v3x':v3x, 'v3y':v3y, 'v4x':v4x, 'v4y':v4y, 'v5x':v5x, 'v5y':v5y, 'v6x':v6x, 'v6y':v6y, 'v7x':v7x, 'v7y':v7y, 'v8x':v8x, 'v8y':v8y, 'dz':dz, 'lunit':'mm'},[]])

    def addSphereInMM(self, name, rmin, rmax, startphi, deltaphi, starttheta, deltatheta):
        self.solids[2].append(['sphere',{'name':name, 'rmin':rmin, 'rmax':rmax,
                                         'startphi':startphi, 'deltaphi':deltaphi,
                                         'starttheta':starttheta, 'deltatheta':deltatheta,
                                         'aunit':'deg', 'lunit':'mm'},[]])

    def addConeInMM(self, name, z, rmin1, rmin2, rmax1, rmax2, sphi, dphi):
        self.solids[2].append(['cone',{'name':name, 'z':z, 'rmin1':rmin1, 'rmin2':rmin2,
                                       'rmax1':rmax1, 'rmax2':rmax2,
                                       'startphi':sphi, 'deltaphi':dphi, 'lunit':'mm', 'aunit':'deg'}, []] )

    def addParaInMM(self, name, x, y, z, alpha, theta, phi):
        self.solids[2].append(['para',{'name':name, 'x':x, 'y':y, 'z':z,
                                       'alpha':alpha, 'theta':theta, 'phi':phi, 'lunit':'mm', 'aunit':'deg'}, []] )

    def addTrapInMM(self, name, z, theta, phi, y1, x1, x2, alpha1, y2, x3, x4, alpha2):
        self.solids[2].append(['trap', {'name':name, 'z':z, 'theta':theta, 'phi':phi,
                                        'x1':x1, 'x2':x2, 'x3':x3, 'x4':x4,
                                        'y1':y1, 'y2':y2, 'alpha1':alpha1, 'alpha2':alpha2, 'lunit':'mm', 'aunit':'deg'}, []])
					
    def addTwistedTrapInMM(self, name, z, theta, phi, y1, x1, x2, alpha1, y2, x3, x4, alpha2, twist):
        self.solids[2].append(['twistTrap', {'name':name, 'z':z, 'theta':theta, 'phi':phi,
                                             'x1':x1, 'x2':x2, 'x3':x3, 'x4':x4,
                                             'y1':y1, 'y2':y2, 'alpha1':alpha1, 'alpha2':alpha2, 'twist':twist, 'aunit':'deg', 'lunit':'mm'}, []])

    def addTrdInMM(self, name, x1, x2, y1, y2, z):
        self.solids[2].append(['trd',{'name':name, 'x1':x1, 'x2':x2,
                                      'y1':y1, 'y2':y2, 'z':z, 'lunit':'mm'}, []])

    def addTubeInMM(self, name, rmin, rmax, z, startphi, deltaphi):
        self.solids[2].append(['tube',{'name':name, 'rmin':rmin, 'rmax':rmax,
                                       'z':z, 'startphi':startphi, 'deltaphi':deltaphi, 'lunit':'mm', 'aunit':'deg'},[]])
				       
    def addCutTubeInMM(self, name, rmin, rmax, z, startphi, deltaphi, lowX, lowY, lowZ, highX, highY, highZ):
        self.solids[2].append(['cutTube',{'name':name, 'rmin':rmin, 'rmax':rmax,
                                          'z':z, 'startphi':startphi, 'deltaphi':deltaphi,
					  'lowX':lowX, 'lowY':lowY, 'lowZ':lowZ, 'highX':highX, 'highY':highY, 'highZ':highZ, 'lunit':'mm', 'aunit':'deg'},[]])

    def addPolyconeInMM(self, name, startphi, deltaphi, zplanes):
        zpls = []
        for zplane in zplanes:
            zpls.append( ['zplane',{'z':zplane[0], 'rmin':zplane[1], 'rmax':zplane[2]},[]] )
        self.solids[2].append(['polycone',{'name':name,
                                           'startphi':startphi, 'deltaphi':deltaphi, 'lunit':'mm', 'aunit':'deg'}, zpls])

    def addTorusInMM(self, name, r, rmin, rmax, startphi, deltaphi):
        self.solids[2].append( ['torus',{'name':name, 'rtor':r, 'rmin':rmin, 'rmax':rmax,
                                         'startphi':startphi, 'deltaphi':deltaphi, 'lunit':'mm', 'aunit':'deg'},[]] )

    def addPolyhedraInMM(self, name, startphi, deltaphi, numsides, zplanes):
        zpls = []
        for zplane in zplanes:
            zpls.append( ['zplane',{'z':zplane[0], 'rmin':zplane[1], 'rmax':zplane[2]},[]] )
        self.solids[2].append(['polyhedra',{'name':name,
                                            'startphi':startphi, 'deltaphi':deltaphi,
                                            'numsides':numsides, 'lunit':'mm', 'aunit':'deg'}, zpls])
					    
    def addXtrusionInMM(self, name, vertices, sections):
        elems = []
	for vertex in vertices:
	    elems.append( ['twoDimVertex',{'x':vertex[0], 'y':vertex[1]},[]] )
	for section in sections:
	    elems.append( ['section',{'zOrder':section[0], 'zPosition':section[1], 'xOffset':section[2], 'yOffset':section[3], 'scalingFactor':section[4]},[]] )
	self.solids[2].append(['xtru',{'name':name, 'lunit':'mm'}, elems])

    def addEltubeInMM(self, name, x, y, z):
        self.solids[2].append( ['eltube', {'name':name, 'dx':x, 'dy':y, 'dz':z, 'lunit':'mm'},[]] )

    def addHypeInMM(self, name, rmin, rmax, inst, outst, z):
        self.solids[2].append( ['hype', {'name':name, 'rmin':rmin, 'rmax':rmax,
                                         'inst':inst, 'outst':outst, 'z':z, 'lunit':'mm', 'aunit':'deg'},[]] )

    def addElconeInMM(self, name, dx, dy, zcut, zmax):
        self.solids[2].append(['elcone',{'name':name, 'dx':dx, 'dy':dy, 'zcut':zcut, 'zmax':zmax, 'lunit':'mm'}, []])

    def addPosInMM(self, subels, type, name, v):
        if v[0]!=0.0 or v[1]!=0.0 or v[2]!=0.0:
            subels.append( [type,{'name':name, 'x':v[0], 'y':v[1], 'z':v[2], 'unit':'mm'},[]] )

    def addPosByName(self, subels, type, name):
        if name!='zeropos':
            subels.append( [type,{'ref':name},[]] )

    def addRotByName(self, subels, type, name):
        if name!='zerorot':
            subels.append( [type,{'ref':name},[]] )

    def addUnionInMM(self, name, lname, lposname, ltr, lrotname, lrot, rname, rposname, rtr, rrotname, rrot):
        subels = [['first',{'ref':lname},[]],
                ['second',{'ref':rname},[]]]
        if self.hasPosition(rposname): self.addPosByName(subels, 'position', rposname)
        else: self.addPosInMM(subels, 'position', rposname, rtr)
        if self.hasRotation(rrotname): self.addRotByName(subels, 'rotation', rrotname)
        else: self.addRot(subels, 'rotation', rrotname, rrot)
        if self.hasPosition(lposname): self.addPosByName(subels, 'firstposition', lposname)
        else: self.addPosInMM(subels, 'firstposition', lposname, ltr)
        if self.hasRotation(lrotname): self.addRotByName(subels, 'firstrotation', lrotname)
        else: self.addRot(subels, 'firstrotation', lrotname, lrot)
        self.solids[2].append( ['union',{'name':name}, subels])

    def addSubtractionInMM(self, name, lname, lposname, ltr, lrotname, lrot, rname, rposname, rtr, rrotname, rrot):
        subels = [['first',{'ref':lname},[]],
                  ['second',{'ref':rname},[]]]
        if self.hasPosition(rposname): self.addPosByName(subels, 'position', rposname)
        else: self.addPosInMM(subels, 'position', rposname, rtr)
        if self.hasRotation(rrotname): self.addRotByName(subels, 'rotation', rrotname)
        else: self.addRot(subels, 'rotation', rrotname, rrot)
        if self.hasPosition(lposname): self.addPosByName(subels, 'firstposition', lposname)
        else: self.addPosInMM(subels, 'firstposition', lposname, ltr)
        if self.hasRotation(lrotname): self.addRotByName(subels, 'firstrotation', lrotname)
        else: self.addRot(subels, 'firstrotation', lrotname, lrot)
        self.solids[2].append( ['subtraction',{'name':name}, subels])

    def addIntersectionInMM(self, name, lname, lposname, ltr, lrotname, lrot, rname, rposname, rtr, rrotname, rrot):
        subels = [['first',{'ref':lname},[]],
                  ['second',{'ref':rname},[]]]
        if self.hasPosition(rposname): self.addPosByName(subels, 'position', rposname)
        else: self.addPosInMM(subels, 'position', rposname, rtr)
        if self.hasRotation(rrotname): self.addRotByName(subels, 'rotation', rrotname)
        else: self.addRot(subels, 'rotation', rrotname, rrot)
        if self.hasPosition(lposname): self.addPosByName(subels, 'firstposition', lposname)
        else: self.addPosInMM(subels, 'firstposition', lposname, ltr)
        if self.hasRotation(lrotname): self.addRotByName(subels, 'firstrotation', lrotname)
        else: self.addRot(subels, 'firstrotation', lrotname, lrot)
        self.solids[2].append( ['intersection',{'name':name}, subels])
