# Driver script to write multiple GDML files:
#     1..cooker_RootNoTor.gdml  --- for the EventDisplay with the toroid removed
#     2..cooker_Fitting.gdml    --- for use with the track fitting plugins
#     2..cooker_Root.gdml       --- for the EventDisplay
#     3..cooker_Geant4.gdml     --- for the Monte Carlo
#
# The files generateGDML.py and.cookerWriter.py should be stored in the same folder as this script.
# To import the module.cookerWriter, you must add $ROOTSYS/geom/gdml to your PYTHONPATH environment variable.  Verify that the file writer.py is in that folder.
# This expects a set of text files describing component geometries to be stored in the same directory.
#
#
#
# VERSION 1.2
# May 24, 2012
# C. O'Connor
#
# Written on a Mac running Python 2.6.1 and GCC 4.2.1 on Darwin
###########################################################################################################################################################

import sys
from gdmlModule import *

###########################################################################################################################################################
#    Run as a script
#
#
#
#    Optional arguments:
#
#        -tag tag
#        Append tag to the end of each generated filename.
#
#        -remove volume1 volume2 ...
#        Remove the listed logical volumes from the geometry and write GDML files that omit them.  Any number of logical volumes may be listed.
#
#        -world worldname
#        Use the volume worldname as your geometry's world volume.  This shouldn't be used for.cooker as the default is correct.
#
#        -posfile posfile
#        Use positions as defined in the file posfile (include path).
#
#        -rotfile rotfile
#        Use rotations as defined in the file rotfile (include path).
#
#        -solfile solfile
#        Use solids as defined in the file solfile (include path).
#
#        -volfile volfile
#        Use volumes as defined in the file volfile (include path).  This applies to the Root and Geant4 GDML files but does not override the volfile used for the FittingGeant4 file.
#
#        -nominal
#        Use the posfile and rotfile that describe the nominal.cooker geometry.
#
#        -survey1
#        Use the posfile and rotfile that describe the.cooker geometry based on the first round of surveys.
#
#        -survey2
#        Use the solfile, volfile, posfile and rotfile that describe the.cooker geometry after the second SYMB survey.  This is the default.
###########################################################################################################################################################

print "\n\nWelcome!  Let's generate some beautiful GDML files together.\n\n"

# set defaults

tag = ''
worldname='World_log0'
omit = []

txtpath = '../human_readable/' # relative path
matfile = txtpath + 'materials.txt'
solfile = txtpath + 'solids2.txt'
volfile = txtpath + 'volumes2.txt'
posfile = txtpath + 'positions_survey2.txt'
rotfile = txtpath + 'rotations_survey2.txt'

###############################################################################
# Handle command line options
#
#

if  '-nominal' in sys.argv:
    solfile = txtpath + 'solids1.txt'
    volfile = txtpath + 'volumes1.txt'
    posfile = txtpath + 'positions_nominal.txt'
    rotfile = txtpath + 'rotations_nominal.txt'

if  '-survey1' in sys.argv:
    solfile = txtpath + 'solids1.txt'
    volfile = txtpath + 'volumes1.txt'
    posfile = txtpath + 'positions_survey1.txt'
    rotfile = txtpath + 'rotations_survey1.txt'

if  '-survey2' in sys.argv:
    solfile = txtpath + 'solids2.txt'
    volfile = txtpath + 'volumes2.txt'
    posfile = txtpath + 'positions_survey2.txt'
    rotfile = txtpath + 'rotations_survey2.txt'

if '-tag' in sys.argv:
    tag = sys.argv[sys.argv.index('-tag')+1]

if '-remove' in sys.argv:
    firstvol = sys.argv.index('-remove') + 1
    lastvol = len(sys.argv) - 1
    for i in range(firstvol, len(sys.argv)):
        if sys.argv[i][0]=='-':
            lastvol = i - 1
            break
    if firstvol<=lastvol:
        for i in range(firstvol, lastvol+1):
            omit.append(sys.argv[i])
        print 'Omitting ' + str(omit) + '\n'

if '-world' in sys.argv:
    worldname = sys.argv[sys.argv.index('-world')+1]

if  '-posfile' in sys.argv:
    posfile = sys.argv[sys.argv.index('-posfile')+1]

if  '-rotfile' in sys.argv:
    rotfile = sys.argv[sys.argv.index('-rotfile')+1]

if  '-solfile' in sys.argv:
    solfile = sys.argv[sys.argv.index('-solfile')+1]

if  '-volfile' in sys.argv:
    volfile = sys.argv[sys.argv.index('-volfile')+1]

###############################################################################
# Build GDML files
#
#

# gdml files to be built are represented by these.cookerWriter objects:
#
#    gdmlRoot
#    gdmlGeant4
#    gdmlFittingGeant4
#    gdmlBigFittingGeant4
#
# modes built by bashscript:
#
#    nominal
#    survey1
#    survey2
#    nominal_WC_cells
#    survey1_WC_cells
#    survey2_WC_cells

# build the gdmlRoot.cookerWriter

gdmlfileRoot = '../gdml_files/Root' + tag + '.gdml'
gdmlRoot =.cookerWriter.cookerWriter(gdmlfileRoot)
gdmlRoot.worldname = worldname

print '\nGenerating ' + gdmlRoot.gdmlfile
print 'World volume: ' + gdmlRoot.worldname

print '\nBuilding positions and rotations...'
print 'Using positions as defined in ' + posfile
print 'Using rotations as defined in ' + rotfile
buildPositionsRotations(gdmlRoot, posfile, rotfile)
    
print '\nBuilding materials...'
print 'Using materials as defined in ' + matfile
nongaseous = []
buildMaterials(gdmlRoot, matfile, nongaseous)

print '\nBuilding solids...'
print 'Using solids as defined in ' + solfile
gdmlRoot.ignoreNoSolid()
gdmlRoot.ignoreSolid('Rear_Flare0')
gdmlRoot.ignoreSolid('Plane_box_solid0')
gdmlRoot.ignoreSolids(['Big_TC_solid0','Big_TC_Win_Box0','Big_Win_1_solid0','Big_Cell_Outer0','Big_Cell_Inner0',
                       'Big_WC_temp_solid0','Big_WC_temp_solid1','Big_WC_gas_temp_solid0','Big_WC_gas_temp_solid1','Big_WC_frame_temp_solid0','Big_WC_frame_temp_solid1',
                       'Big_Chamber_box_solid0','Big_Plane_box_solid0','Big_Superlayer_box_solid0','Big_Layer_box_solid0','Big_Cell_box_solid0','Big_Wire_box_solid0',
                       'Big_WC_Win_solid0','Big_WC_Win_solid1'])
for i in range(0,16):
    gdmlRoot.ignoreSolids(['WS_innersection' + str(i) + '_solid0','WS_outersection' + str(i) + '_solid0'])
buildSolids(gdmlRoot, solfile)

# copy the gdmlRoot.cookerWriter and build the gdmlFittingGeant4.cookerWriter

gdmlfileFittingGeant4 = '../gdml_files/FittingGeant4' + tag + '.gdml'
print '\nCopying ' + gdmlRoot.gdmlfile + ' content to ' + gdmlfileFittingGeant4
gdmlFittingGeant4 = gdmlRoot.copy()
gdmlFittingGeant4.worldname = worldname
gdmlFittingGeant4.setFile(gdmlfileFittingGeant4)

# copy the gdmlRoot.cookerWriter and build the gdmlBigFittingGeant4.cookerWriter

gdmlfileBigFittingGeant4 = '../gdml_files/BigFittingGeant4' + tag + '.gdml'
print '\nCopying ' + gdmlRoot.gdmlfile + ' content to ' + gdmlfileBigFittingGeant4
gdmlBigFittingGeant4 = gdmlRoot.copy()
gdmlBigFittingGeant4.worldname = worldname
gdmlBigFittingGeant4.setFile(gdmlfileBigFittingGeant4)

# finish building the gdmlRoot.cookerWriter

print '\nNow back to working on ' + gdmlRoot.gdmlfile

print '\nBuilding volumes...'
print 'Using volumes as defined in ' + volfile
gdmlRoot.ignoreNoVolume()
gdmlRoot.ignoreVolumes(omit)
buildVolumes(gdmlRoot, volfile, False)

# copy the gdmlRoot.cookerWriter and build the gdmlGeant4.cookerWriter

gdmlfileGeant4 = '../gdml_files/Geant4' + tag + '.gdml'
print '\nCopying ' + gdmlRoot.gdmlfile + ' content to ' + gdmlfileGeant4
gdmlGeant4 = gdmlRoot.copy()
gdmlGeant4.worldname = worldname
gdmlGeant4.setFile(gdmlfileGeant4)

print '\nGenerating ' + gdmlGeant4.gdmlfile + ' ...'
print 'World volume: ' + gdmlGeant4.worldname

# add the elcone solid and related volumes for the collimator

print '\nBuilding the elcone and related compound solids...'
print 'Using solids as defined in ' + solfile
gdmlGeant4.ignoreNoSolid()
gdmlGeant4.ignoreSolid('Plane_box_solid0')
gdmlGeant4.ignoreSolids(['Big_TC_solid0','Big_TC_Win_Box0','Big_Win_1_solid0','Big_Cell_Outer0','Big_Cell_Inner0',
                       'Big_WC_temp_solid0','Big_WC_temp_solid1','Big_WC_gas_temp_solid0','Big_WC_gas_temp_solid1','Big_WC_frame_temp_solid0','Big_WC_frame_temp_solid1',
                       'Big_Chamber_box_solid0','Big_Plane_box_solid0','Big_Superlayer_box_solid0','Big_Layer_box_solid0','Big_Cell_box_solid0','Big_Wire_box_solid0',
                       'Big_WC_Win_solid0','Big_WC_Win_solid1'])
buildSolids(gdmlGeant4, solfile)

print '\nBuilding the collimator volumes...'
print 'Using volumes as defined in ' + volfile
gdmlGeant4.removeVolume('TC_log0')
gdmlGeant4.removeVolume('WS_log0')
gdmlGeant4.ignoreNoVolume()
gdmlGeant4.ignoreVolumes(omit)
buildVolumes(gdmlGeant4, volfile, False)

# go back to gdmlRoot and add the invisible WC planes

print '\nReturning to ' + gdmlRoot.gdmlfile + ' for a moment...'
print '\nAdding invisible WC planes'
gdmlRoot.ignoreNoSolid()
gdmlRoot.ignoreSolid('Rear_Flare0')
gdmlRoot.ignoreSolids(['Big_TC_solid0','Big_TC_Win_Box0','Big_Win_1_solid0','Big_Cell_Outer0','Big_Cell_Inner0',
                       'Big_WC_temp_solid0','Big_WC_temp_solid1','Big_WC_gas_temp_solid0','Big_WC_gas_temp_solid1','Big_WC_frame_temp_solid0','Big_WC_frame_temp_solid1',
                       'Big_Chamber_box_solid0','Big_Plane_box_solid0','Big_Superlayer_box_solid0','Big_Layer_box_solid0','Big_Cell_box_solid0','Big_Wire_box_solid0',
                       'Big_WC_Win_solid0','Big_WC_Win_solid1'])
for i in range(0,16):
    gdmlRoot.ignoreSolids(['WS_innersection' + str(i) + '_solid0','WS_outersection' + str(i) + '_solid0'])
buildSolids(gdmlRoot, solfile)

for i in range(0,3):
    gdmlRoot.removeVolume('WC_chamber_log' + str(i))
    gdmlRoot.removeVolume('WC_chamber_2_log' + str(i))
gdmlRoot.removeVolume('WC_gas_log0')
gdmlRoot.removeVolume('WC_gas_2_log0')
gdmlRoot.removeVolume('WC_log0')
gdmlRoot.removeVolume('WC_2_log0')
gdmlRoot.ignoreNoVolume()
gdmlRoot.ignoreVolumes(omit)
buildVolumes(gdmlRoot, volfile, False)

# build the gdmlFittingGeant4.cookerWriter

print '\nGenerating ' + gdmlFittingGeant4.gdmlfile
print 'World volume: ' + gdmlFittingGeant4.worldname

fitvolfile = txtpath + 'volumes_fitting.txt'

print '\nBuilding volumes...'
print 'Using volumes as defined in ' + fitvolfile + ' (automatic for this GDML file)'
gdmlFittingGeant4.ignoreNoVolume()
gdmlFittingGeant4.ignoreVolumes(omit)
buildVolumes(gdmlFittingGeant4, fitvolfile, False)

# build the gdmlBigFittingGeant4.cookerWriter

print '\nGenerating ' + gdmlBigFittingGeant4.gdmlfile
print 'World volume: ' + gdmlBigFittingGeant4.worldname

bigfitvolfile = txtpath + 'volumes_big_fitting.txt'

print '\nBuilding the Big solids...'
print 'Using solids as defined in ' + solfile
gdmlBigFittingGeant4.ignoreNoSolid()
buildSolids(gdmlBigFittingGeant4, solfile)

print '\nBuilding volumes...'
print 'Using volumes as defined in ' + bigfitvolfile + ' (automatic for this GDML file)'
gdmlBigFittingGeant4.ignoreNoVolume()
gdmlBigFittingGeant4.ignoreVolumes(omit)
buildVolumes(gdmlBigFittingGeant4, bigfitvolfile, False)
    
# for each GDML file:
#     add the World volume, incorporating whichever volumes are appropriate
#     call addSetup(name, version, world)
#     write the file

olympusWriters = {gdmlRoot: volfile,
                  gdmlGeant4: volfile,
                  gdmlFittingGeant4: fitvolfile,
                  gdmlBigFittingGeant4: bigfitvolfile}
sub = 0

for owriter in.cookerWriters:
    print '\nBuilding the World volume for ' + owriter.gdmlfile
    buildVolumes(owriter,.cookerWriters[owriter])
    owriter.addSetup('Default', '3.' + str(sub), owriter.worldname)
    print 'Writing ' + owriter.gdmlfile
    owriter.writeFile()
    sub += 1

print "\nLooks like we made it... after all!"
