# Function library for writing multiple GDML files:
#     1..cooker_RootNoTor.gdml  --- for the EventDisplay with the toroid removed
#     2..cooker_Fitting.gdml    --- for use with the track fitting plugins
#     2..cooker_Root.gdml       --- for the EventDisplay
#     3..cooker_Geant4.gdml     --- for the Monte Carlo
#
# The file.cookerWriter.py should be stored in the same folder as this script.
# To import the module.cookerWriter, you must add $ROOTSYS/geom/gdml to your PYTHONPATH environment variable.  Verify that the file writer.py is in that folder.
# This expects a set of text files describing component geometries to be stored in the same directory.
#
#
#
# VERSION 1.1
# April 16, 2012
# C. O'Connor
#
# Written on a Mac running Python 2.6.1 and GCC 4.2.1 on Darwin
###########################################################################################################################################################

import.cookerWriter

###########################################################################################################################################################
# set up "define" using:
#     addPositionInMM(name, x, y, z)
#     addRotation(name, x, y, z)
###########################################################################################################################################################

def buildPositionsRotations(owriter, posfile, rotfile):
    
    # Open the data file for positions.  Its structure must match what is expected.
    # All lengths are assumed to be input in mm.
    
    infile = open(posfile, 'r')
    f = infile.readlines()
    i = 1

    # Scan the whole contents of the file

    while i < len(f):
        
        l = f[i].split()
        if l==[]:
            i += 1
            continue
        elif l[0]=='Position':
            i += 1
            continue
        else:
            # interpret a named (x,y,z) position
            assert len(l)==5
            name = l[0]
            x = float(l[1])
            y = float(l[2])
            z = float(l[3])
            if l[4]!='mm':
                print 'Position ' + name + ' not in millimeters --- please put it into millimeters in the positions source file! The GDML produced will be wrong.'
            owriter.addPositionInMM(name, x, y, z)
                
        i += 1

    # Close the data file

    infile.close()

    # Open the data file for rotations.  Its structure must match what is expected.
    # All angles are assumed to be input in degrees.
    
    infile = open(rotfile, 'r')
    f = infile.readlines()
    i = 1

    # Scan the whole contents of the file

    while i < len(f):
        
        l = f[i].split()
        if l==[]:
            i += 1
            continue
        elif l[0]=='Rotation':
            i += 1
            continue
        else:
            # interpret a named (x,y,z) rotation
            assert len(l)==5
            name = l[0]
            x = float(l[1])
            y = float(l[2])
            z = float(l[3])
            if l[4]!='deg':
                print 'Rotation ' + name + ' not in degrees --- please put it into degrees in the positions source file! The GDML produced will be wrong.'
            owriter.addRotation(name, x, y, z)
                
        i += 1

    # Close the data file

    infile.close()

###########################################################################################################################################################
# set up "materials" using:
#     addMaterial(name, a, z, rho) unit of density (rho) is fixed by a GDML default; it's PROBABLY g/cm^3 but I'm not 100% sure
#     addMixture(name, rho, elems) where elems is a Dictionary; e.g., 90:10 Ar:CO2 gas would have elems = {'Argon':0.9, 'CarbonDioxide':0.1}
#     addElement(symbol, name, z, a)
#
#
#
#     Additional methods defined for.cookerWriter:
#
#     addIsotope(name, n, z, a)
#     addElementFromIsotopes(name, isotopes) where isotopes is a Dictionary of the same format as elems in addMixture
#     addMaterialWithAttributes(name, state, attribs, elems) where attribs is a Dictionary; e.g., attribs = {attrib1:[unit1, value1], attrib2:[unit2, value2], ...} and elems is a Dictionary as in addMixture
#         Note: attributes should be chosen from the list 'MEE', 'D', 'P', 'T', 'atom' with corresponding units (in.cooker) of 'eV', 'g/cm3', 'pascal', 'K', 'g/mol' respectively
#     addMaterialG4Al(name='G4_Al', state='solid', a=26.9815, z=13, rho=2.699, mee=166)
#         Note: this material was formatted differently from any other, so I've made a tailored constructor just for it
###########################################################################################################################################################

def buildMaterials(owriter, matfile, nongaseous):

    # Open the data file for materials.  Its structure must match what is expected.
    
    infile = open(matfile, 'r')
    f = infile.readlines()
    i = 1

    # Scan the whole contents of the file

    while i < len(f):
        
        l = f[i].split()
        if l==[]:
            i += 1
            continue
        materialtype = l[0]
        materialtype = materialtype[:len(materialtype)-1]
        if materialtype: name = l[1]
        i += 1
        l = f[i].split()

        # Parse and add isotopes
        
        if materialtype=='G4Isotope':
            
            n = int(l[2])
            i += 1
            l = f[i].split()
            z = int(l[2])
            i += 1
            l = f[i].split()
            a = float(l[2])
            owriter.addIsotope(name, n, z, a)
            i += 1

        # Parse and add elements

        elif materialtype=='G4Element':

            assert len(l)==3
            i += 1
            l = f[i].split()

            # add each isotope component
            isotopes = []
            while l!=[]:
                assert len(l)==2
                isoname = l[0]
                isoname = isoname[:len(isoname)-1]
                fraction = float(l[1])
                isotopes.append((isoname, fraction))
                i += 1
                if i < len(f):
                    l = f[i].split()
                else:
                    break

            # write the complete element info
            isodict = dict(isotopes)
            owriter.addElementFromIsotopes(name, isodict)

        # Parse and add materials

        elif materialtype=='G4Material':

            # Parse and add the uniquely defined material G4_Al

            assert len(l)==2
            state = l[1]
            if state!='gas':
                nongaseous.append(name)
            i += 1
            l = f[i].split()

            # add each attribute
            attributes = []
            while l[0]!='Elements':
                assert len(l)==4
                att = l[0]
                value = float(l[2])
                unit = l[3]
                unit = unit[1:len(unit)-1]
                attributes.append((att, [unit, value]))
                i += 1
                l = f[i].split()

            i += 1
            l = f[i].split()

            # add each element component
            elements = []
            while l!=[]:
                assert len(l)==2
                elem = l[0]
                elem = elem[:len(elem)-1]
                fraction = float(l[1])
                elements.append((elem, fraction))
                i += 1
                if i < len(f):
                    l = f[i].split()
                else:
                    break

            # write the complete material info
            attdict = dict(attributes)
            elemdict = dict(elements)
            owriter.addMaterialWithAttributes(name, state, attdict, elemdict)

        i += 1

    # Close the data file

    infile.close()

###########################################################################################################################################################
# set up "solids" using: EVERYTHING IS IN CENTIMETERS AND DEGREES as fixed by writer.py
#
#     Methods defined in writer that aren't needed for.cooker:
#
#     addReflSolid(name, solid, dx, dy, dz, sx, sy, sz, rx, ry, rz)
#     addParaboloid(name, rlo, rhi, dz)
#     addArb8(name, v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, v5x, v5y, v6x, v6y, v7x, v7y, v8x, v8y, dz)
#     addSphere(name, rmin, rmax, startphi, deltaphi, starttheta, deltatheta)
#     addPara(name, x, y, z, alpha, theta, phi)
#     addTwistedTrap(name, z, theta, phi, y1, x1, x2, alpha1, y2, x3, x4, alpha2, twist)
#     addCutTube(name, rmin, rmax, z, startphi, deltaphi, lowX, lowY, lowZ, highX, highY, highZ)
#     addPolycone(name, startphi, deltaphi, zplanes) where zplanes is a List of Lists; e.g., zplanes = [[zplane1z, zplane1rmin, zplane1rmax], [zplane2z, zplane2rmin, zplane2rmax], ...]
#     addTorus(name, r, rmin, rmax, startphi, deltaphi)
#     addPolyhedra(name, startphi, deltaphi, numsides, zplanes) where zplanes is a List of Lists; ""
#     addXtrusion(name, vertices, sections): where vertices is a List of Lists; e.g., vertices = [[vertex1x, vertex1y], [vertex2x, vertex2y], ...] and sections is a List of Lists; e.g., sections = [[section1zOrder, section1zPosition, section1xOffset, section1yOffset, section1scalingFactor], ...]
#     addHype(name, rmin, rmax, inst, outst, z)
#
#
#
#     Methods we will use:
#
#     addBox(name, dx, dy, dz)
#     addCone(name, z, rmin1, rmin2, rmax1, rmax2, sphi, dphi)
#     addTrap(name, z, theta, phi, y1, x1, x2, alpha1, y2, x3, x4, alpha2)
#     addTrd(name, x1, x2, y1, y2, z)
#     addTube(name, rmin, rmax, z, startphi, deltaphi)
#     addEltube(name, x, y, z)
#
#
#
#     Combined solids:
#
#     addUnion(name, lname, ltr, lrot, rname, rtr, rrot) where ltr, lrot, rtr, and rrot are each Lists containing three doubles, "l" referring to the first solid and "r" to the second
#     addSubtraction(name, lname, ltr, lrot, rname, rtr, rrot) where ""
#     addIntersection(name, lname, ltr, lrot, rname, rtr, rrot) where ""
#
#
#
#     Additional methods defined for.cookerWriter:
#
#     addElcone(name, dx, dy, zcut, zmax)
#     hasSolid(name)
#
#
#
#     Add "InMM" to the end of a function name to make it use better units; the default is centimeters.
###########################################################################################################################################################

def buildSolids(owriter, solfile, windowscale=1):

    # Open the data file for solids.  Its structure must match what is expected.
    
    infile = open(solfile, 'r')
    f = infile.readlines()
    i = 1

    # Scan the whole contents of the file

    while i < len(f):
        
        l = f[i].split()
        if l==[]:
            i += 1
            continue
        solidtype = l[0]
        solidtype = solidtype[:len(solidtype)-1]
        name = ''
        if solidtype and len(l)>1: name = l[1]
        i += 1
        if i >= len(f): break
        l = f[i].split()

        if not name: continue

        # Don't add any solids we already have
        if owriter.hasSolid(name): continue

        # Don't add any solids we're ignoring for this file
        if owriter.isIgnoringSolid(name):
            print 'Solid ' + name + " skipped because it doesn't belong in " + owriter.gdmlfile
            continue

        # Parse and add simple solids

        if solidtype in ['G4Box', 'G4Cone', 'G4Trap', 'G4Trd', 'G4Tube', 'G4Eltube', 'G4Elcone', 'G4GenericTrap']:

            # add each parameter needed for this solid
            args = []
            while l!=[]:
                assert len(l)==4
                argval = l[2]
                args.append(float(argval))
                argunit = l[3]
                argunit = argunit[1:len(argunit)-1]
                if argunit not in ['mm', 'deg']:
                    print 'Solid ' + name + ' skipped because its units are wrong --- please put it into millimeters and degrees in the solids source file! The GDML produced will be wrong.'
                    owriter.ignoreSolid(name)
                i += 1
                if i >= len(f): break
                l = f[i].split()

            # scale the target chamber windows
            if name=='TC_Win_Box0':
                args[1] = args[1]*windowscale
                args[2] = args[2]*windowscale
            elif name=='Win_1_solid0':
                args[0] = args[0]*windowscale
                args[1] = args[1]*windowscale
            elif name in ['TC_Win_Cyl0','Win_2_solid0']:
                args[1] = args[1]*windowscale

            # write the appropriate solid
            if solidtype=='G4Box':
                assert len(args)==3
                owriter.addBoxInMM(name, args[0], args[1], args[2])                                                                           # in mm or degrees
            if solidtype=='G4GenericTrap':
                assert len(args)==17
                owriter.addArb8InMM(name, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13], args[14], args[15], args[16]) # in mm or degrees
            elif solidtype=='G4Cone':
                assert len(args)==7
                owriter.addConeInMM(name, args[0], args[1], args[2], args[3], args[4], args[5], args[6])                                      # in mm or degrees
            elif solidtype=='G4Trap':
                assert len(args)==11
                owriter.addTrapInMM(name, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10]) # in mm or degrees
            elif solidtype=='G4Trd':
                assert len(args)==5
                owriter.addTrdInMM(name, args[0], args[1], args[2], args[3], args[4])                                                          # in mm or degrees
            elif solidtype=='G4Tube':
                assert len(args)==5
                owriter.addTubeInMM(name, args[0], args[1], args[2], args[3], args[4])                                                         # in mm or degrees
            elif solidtype=='G4Eltube':
                assert len(args)==3
                owriter.addEltubeInMM(name, args[0], args[1], args[2])                                                                         # in mm or degrees
            elif solidtype=='G4Elcone':
                assert len(args)==4
                owriter.addElconeInMM(name, args[0], args[1], args[2], args[3])                                                                # in mm or degrees

        # Parse and add compound solids

        elif solidtype in ['G4UnionSolid', 'G4SubtractionSolid', 'G4IntersectionSolid']:

            assert len(l)==3
            aname = l[2]
            i += 1
            l = f[i].split()
            assert len(l)==3
            bname = l[2]
            aposname = 'zeropos'
            arotname = 'zerorot'
            bposname = 'zeropos'
            brotname = 'zerorot'
            apos = [0.0, 0.0, 0.0] # relative position of solid a
            arot = [0.0, 0.0, 0.0] # relative rotation of solid a
            bpos = [0.0, 0.0, 0.0] # relative position of solid b
            brot = [0.0, 0.0, 0.0] # relative rotation of solid b
            i += 1
            if i >= len(f): break
            l = f[i].split()

            if owriter.hasSolid(aname) and owriter.hasSolid(bname):
                
                while l!=[]:

                    # add relative position
                    if l[0]=='Position':
                        bposname = l[6]
                        if not owriter.hasPosition(bposname):
                            i += 1
                            l = f[i].split()
                            assert len(l)==8
                            posx = l[4]
                            posx = posx[1:len(posx)-1]
                            posy = l[5]
                            posy = posy[:len(posy)-1]
                            posz = l[6]
                            posz = posz[:len(posz)-1]
                            posunit = l[7]
                            posunit = posunit[1:len(posunit)-1]
                            if posunit!='mm':
                                print 'Solid ' + name + ' skipped because its units are wrong --- please put it into millimeters and degrees in the solids source file! The GDML produced will be wrong.'
                                owriter.ignoreSolid(name)
                            bpos[0] = float(posx) # in mm
                            bpos[1] = float(posy) # in mm
                            bpos[2] = float(posz) # in mm
                        i += 1
                        if i >= len(f): break
                        l = f[i].split()

                    # add relative rotation
                    elif l[0]=='Rotation':
                        brotname = l[6]
                        if not owriter.hasRotation(brotname):
                            i += 1
                            l = f[i].split()
                            assert len(l)==8
                            rotx = l[4]
                            rotx = rotx[1:len(rotx)-1]
                            roty = l[5]
                            roty = roty[:len(roty)-1]
                            rotz = l[6]
                            rotz = rotz[:len(rotz)-1]
                            rotunit = l[7]
                            rotunit = rotunit[1:len(rotunit)-1]
                            if rotunit!='deg':
                                print 'Solid ' + name + ' skipped because its units are wrong --- please put it into millimeters and degrees in the solids source file! The GDML produced will be wrong.'
                                owriter.ignoreSolid(name)
                            brot[0] = float(rotx) # in mm
                            brot[1] = float(roty) # in mm
                            brot[2] = float(rotz) # in mm
                        i += 1
                        if i >= len(f): break
                        l = f[i].split()

                # scale the target chamber windows
                if name in ['Win_3_solid0','Win_solid0']:
                    bpos[0] = bpos[0]*windowscale
                elif name in ['TC_F4_solid0','TC_F5_solid0','TC_Frame_solid0']:
                    bpos[2] = bpos[2]*windowscale

                # write the appropriate solid
                if solidtype=='G4UnionSolid': owriter.addUnionInMM(name, aname, aposname, apos, arotname, arot, bname, bposname, bpos, brotname, brot)                 # in mm or degrees
                elif solidtype=='G4SubtractionSolid': owriter.addSubtractionInMM(name, aname, aposname, apos, arotname, arot, bname, bposname, bpos, brotname, brot)   # in mm or degrees
                elif solidtype=='G4IntersectionSolid': owriter.addIntersectionInMM(name, aname, aposname, apos, arotname, arot, bname, bposname, bpos, brotname, brot) # in mm or degrees

            # if a component solid is missing because we're intentionally ignoring it, we should probably ignore this compound solid too
            elif owriter.isIgnoringSolid(aname) or owriter.isIgnoringSolid(bname):        
                owriter.ignoreSolid(name)
                print 'Compound solid ' + name + ' skipped because its component ' + aname + ' or ' + bname + " doesn't belong in " + owriter.gdmlfile
            else:
                print "Can't find solids " + aname + ' or ' + bname + ' to make ' + name

        i += 1

    # Close the data file

    infile.close()

###########################################################################################################################################################
# set up "structures" using:
#     addVolume(name, solid, material, daughters) where daughters is a List of Lists; e.g., daughters = [['wire0vol', 'wire0pos', 'wire0rot'], ['wire1vol', 'wire1pos', 'wire1rot'], ...]
#     addAssembly(name, daughters) where daughters is a List of Lists; ""
###########################################################################################################################################################

def buildVolumes(owriter, volfile, buildworld=True, gasify=[]):

    # Open the data file for volumes.  Its structure must match what is expected.
    
    infile = open(volfile, 'r')
    f = infile.readlines()
    i = 1

    # Scan the whole contents of the file

    while i < len(f):
        
        l = f[i].split()
        if l==[]:
            i += 1
            continue
        volumetype = l[0]
        name = ''
        sensdet = ''
        color = 'f999999999'
        copyno = '0'
        if volumetype=='Logical': name = l[2]
        i += 1
        if i >= len(f): break
        l = f[i].split()

        # Parse and add volumes

        if not name:
            i += 1
            continue

        if owriter.hasVolume(name):
            i += 1
            continue
            
        assert l[0]=='Solid:'
        solid = l[1]
        i += 1
        l = f[i].split()
        assert l[0]=='Material:'
        material = l[1]
        daughters = []
        i += 1
        l = f[i].split()
        assert l[0]=='Color:'
        color = l[1]
        i += 1
        l = f[i].split()

        # Check for Sensitive Detector
        if l!=[]:
            if l[0]=='Sensitive':
                sensdet = l[2]
                i += 1
                l = f[i].split()

        # Get CopyNo

        if name.find('log') >= 0:

            if name.find('WC_wire') >= 0:
                wireNumberInCell = name[name.find('wire')+4]
                cellNumberInSector = name[name.find('log')+3:]
                if name.find('_2_') >= 0: sectorNumber = '1'
                else: sectorNumber = '0'
                copyno = str(477*int(sectorNumber) + 3*int(cellNumberInSector) + int(wireNumberInCell)) 

            else:
                copyno = name[name.find('log')+3:]

        # Don't build the World volume unless buildworld==True
        if name==owriter.worldname:
            if buildworld:
                pass
            else:
                print 'Not building the World volume yet'
                continue

        # Don't add any volumes we're ignoring for this file
        if owriter.isIgnoringSolid(solid):
            owriter.ignoreVolume(name)
            print 'Volume ' + name + " skipped because its solid doesn't belong in " + owriter.gdmlfile
            continue
        elif owriter.isIgnoringVolume(name):
            print 'Volume ' + name + " skipped because it doesn't belong in " + owriter.gdmlfile
            continue
            
        # Gasify anything that needs to be gasified
        elif material in gasify and name!='Win_log0':
            material = 'G4_AIR'
            print 'Volume ' + name + " G4_AIRified because it isn't made out of gas"
                
        while l!=[]:

            # add daughter volumes
            if l[0]=='Daughter' and l[1]=='(Physical':
                daughtername = l[3]
                i += 1
                l = f[i].split()
                daughtervolume = l[3]
                daughterposition = 'zeropos'
                daughterrotation = 'zerorot'
                i += 1
                l = f[i].split()

                # add positions and rotations for daughter volumes (relative to mother volume)
                while l!=[] and l[1]!='(Physical':

                    if l[1]=='Position:':
                        daughterposition = l[2]

                    elif l[1]=='Rotation:':
                        daughterrotation = l[2]

                    i += 1
                    if i >= len(f): break
                    l = f[i].split()

                if owriter.hasVolume(daughtervolume): daughters.append([daughtername, daughtervolume, daughterposition, daughterrotation])
                elif owriter.isIgnoringVolume(daughtervolume): print 'Daughter Volume ' + daughtername + ' in ' + name + " skipped because it doesn't belong in " + owriter.gdmlfile
                else: print "Can't find volume " + daughtervolume + ' to make ' + name
                    
                if i >= len(f): break

        # Don't add any volumes we already have, unless they have new daughters
        if owriter.hasVolume(name):
            alldaughters = True
            for daughter in daughters:
                if not owriter.hasVolumeWithDaughterVolume(name, daughter[1]):
                    alldaughters = False
            if alldaughters:
                continue
            else:
                owriter.removeVolume(name)

        #print 'Adding: ' + name + ' ' + str(daughters)
        owriter.addVolumeWithNamedDaughters(name, solid, material, color, sensdet, copyno, daughters)
                
        i += 1

    # Close the data file

    infile.close()
