# import libraries
import arcpy
from arcpy.sa import *
from arcpy import env

__author__ = 'senait senay'

''' The functions rescale50to100, rescalecustom and rescale100to0 permanently add a column on the respective processed
files 'gdd', 'lcc' and 'lev'. If fixing subsequent steps is necessary simply continue from the lookup step where
it copies out the permanently processed "NORMALIZED" column from each of these component landscapes. The temporary
cleaning part of the script deletes all the rasters except the original 'gdd', 'lcc' and 'lev' files and the final
result 'landscape1' and 'landscape2'. By using the cleaning section of the script and manually deleting the two final
landscape files one can re-run the script from the look up step forward. This script requires a full license
ArcInfo license (ArcGIS advanced license) as it makes use of spatial analyst functions. However, the script can
easily be adapted to work with rddal, ogr functions.Either GRASS or QGIS are easy to use tools for this purpose.'''

# set the environments
arcpy.env.overwriteOutput = True
# env.workspace = arcpy.GetParameter(0)
env.workspace = "path_to_folder_RawLandscapeComponents"
home = env.workspace


# Define a function to list data in the home workspace

def catws():
    desc = arcpy.Describe(env.workspace)
    if hasattr(desc, "name"):
        print "Name:        " + desc.name
    if hasattr(desc, "dataType"):
        print "DataType:    " + desc.dataType
    if hasattr(desc, "catalogPath"):
        print "CatalogPath: " + desc.catalogPath

    # Examine children and print their name and dataType
    print "Children:"
    for child in desc.children:
        print "\t%s = %s" % (child.name, child.dataType)

# List input landscape Components and the column name that has landscape values (Make sure these three files exist)
gdd = "csi_gdd50"
lcc = "csi_lcc50"
lev = "csi_lev50"
value = "VALUE"  # one column name is listed for all three components as they have the same column name


# function to rescale rater values of between 0.5-1 for GDD

def rescale50to100(rast, val, home):
    arcpy.AddField_management(rast, "NORMALIZED", "DOUBLE")
    scratchTable = home + "\\temptable"
    arcpy.Statistics_analysis(rast, scratchTable, [[str(val), "MIN"], [str(val), "MAX"]])
    cursor = arcpy.SearchCursor(scratchTable, ("MAX_" + str(val), "MIN_" + str(val)))
    for row in cursor:
        maxNum = row.getValue('{0}'.format("MAX_" + str(val)))
        minNum = row.getValue('{0}'.format("MIN_" + str(val)))
    del cursor, row
    arcpy.Delete_management(scratchTable)
    cursor = arcpy.UpdateCursor(rast, ("VALUE", "NORMALIZED"))
    for row in cursor:
        number = row.VALUE
        if number == 0:
            row.NORMALIZED == 0
        else:
            row.NORMALIZED = ((number - (minNum + 1)) / (maxNum - (minNum + 1))) * (1 - 0.5) + 0.5
        cursor.updateRow(row)
    del cursor, row


# function to rescale rater values of between 0 - 1 for LCC

def rescalecustom(rast):
    arcpy.AddField_management(rast, "NORMALIZED", "DOUBLE")
    cursor = arcpy.UpdateCursor(rast, ("VALUE", "NORMALIZED"))
    for row in cursor:
        number = row.VALUE
        if number == 1:
            row.NORMALIZED = float(0.9)
        elif number == 2:
            row.NORMALIZED = float(0.8)
        elif number == 3:
            row.NORMALIZED = float(0.5)
        elif number == 4:
            row.NORMALIZED = float(0.3)
        elif number == 5:
            row.NORMALIZED = float(0.1)
        else:
            row.NORMALIZED = float(0.0)
        cursor.updateRow(row)
    del cursor, row


#  function to rescale rater values of between 1-0 for LEV (inverse)

def rescale100to0(rast, value, home):
    arcpy.AddField_management(rast, "NORMALIZED", "DOUBLE")
    scratchTable = home + "\\temptable"
    arcpy.Statistics_analysis(rast, scratchTable, [[str(value), "MIN"], [str(value), "MAX"]])
    cursor = arcpy.SearchCursor(scratchTable, ("MAX_" + str(value), "MIN_" + str(value)))
    for row in cursor:
        maxNum = row.getValue('{0}'.format("MAX_" + str(value)))
        minNum = row.getValue('{0}'.format("MIN_" + str(value)))
    del cursor, row
    arcpy.Delete_management(scratchTable)
    cursor = arcpy.UpdateCursor(rast, ("VALUE", "NORMALIZED"))
    for row in cursor:
        number = row.VALUE
        row.NORMALIZED = ((number - minNum - maxNum + minNum) / (maxNum - minNum))
        cursor.updateRow(row)
    del cursor, row


arcpy.CheckOutExtension("spatial")  # check out a spatial analyst licence

# Local variables to save landscape component rescaled values
gdd_surv = "gdd_surv"
lcc_surv = "lcc_surv"
lev_surv = "lev_surv"

# Process GDD
rescale50to100(gdd, value, home)  # the rescaled values are permanently added on the layer gdd in column "NORMALIZATION"
gddl = Lookup(gdd, "NORMALIZED")
gddl.save(gdd_surv)

# Process LCC
rescalecustom(lcc)
lccl = Lookup(lcc, "NORMALIZED")  # the rescaled values are permanently added on the layer lcc in column "NORMALIZATION"
lccl.save(lcc_surv)

# Process LEV
lev_temp1 = "lev_temp1"

rescale100to0(lev, value, home)
levl = Lookup(lev, "NORMALIZED")  # the rescaled values are permanently added on the layer lev in column "NORMALIZATION"
levl.save(lev_temp1)
print arcpy.GetRasterProperties_management(lev_temp1, "MINIMUM")
print arcpy.GetRasterProperties_management(lev_temp1, "MAXIMUM")
outAbs = Abs(lev_temp1)
outAbs.save(lev_surv)
print arcpy.GetRasterProperties_management(lev_surv, "MINIMUM")
print arcpy.GetRasterProperties_management(lev_surv, "MAXIMUM")


# Define a function to rescale the final landscapes between 0 and 1 after the components are merged


def rescale0to1(rast, value, home):
    arcpy.AddField_management(rast, "NORMALIZED", "DOUBLE")
    scratchTable = home + "\\temptable"
    arcpy.Statistics_analysis(rast, scratchTable, [[str(value), "MIN"], [str(value), "MAX"]])
    cursor = arcpy.SearchCursor(scratchTable, ("MAX_" + str(value), "MIN_" + str(value)))
    for row in cursor:
        maxNum = row.getValue('{0}'.format("MAX_" + str(value)))
        minNum = row.getValue('{0}'.format("MIN_" + str(value)))
    del cursor, row
    arcpy.Delete_management(scratchTable)
    cursor = arcpy.UpdateCursor(rast, ("VALUE", "NORMALIZED"))
    for row in cursor:
        number = row.VALUE
        row.NORMALIZED = ((number - minNum) / (maxNum - minNum))
        cursor.updateRow(row)
    del cursor, row


#  -------------   PREPARE THE TWO LANDSCAPES --------------------------------#

# -----------------------Landscape 1------------------------------------------#
ls1 = "ls_1raw"

outTimes1 = arcpy.sa.Times(gdd_surv, lcc_surv)  # Multiply gdd by lcc
outTimes1.save(ls1)

# Check min and max
print arcpy.GetRasterProperties_management(ls1, "MINIMUM")
print arcpy.GetRasterProperties_management(ls1, "MAXIMUM")

# multiply the layer by 100 to get survival percentage
LS_1 = "LS_1"
LS_1r = "LS_1r"
LS1_final = "ls1_final"
landscape1 = "landscape1"

outTimes3 = Times(ls1, 100)
outTimes3.save(LS_1)
OutInt = Int(LS_1)
OutInt.save(LS_1r)
rescale0to1(LS_1r, value, home)
lslook = Lookup(LS_1r, "NORMALIZED")
lslook.save(LS1_final)
outTimes5 = Times(LS1_final, 100)
outTimes5.save(landscape1)

print arcpy.GetRasterProperties_management(landscape1, "MINIMUM")
print arcpy.GetRasterProperties_management(landscape1, "MAXIMUM")

# -----------------------Landscape 2------------------------------------------#
ls2 = "ls_2raw"

outTimes2 = arcpy.sa.Times(ls1, lev_surv)  # Multiply ls_1raw (gdd * lcc) by lev
outTimes2.save(ls2)
# check min and max
print arcpy.GetRasterProperties_management(ls2, "MINIMUM")
print arcpy.GetRasterProperties_management(ls2, "MAXIMUM")

# multiply the layer  by 100 to get survival percentage
LS_2 = "LS_2"
LS_2r = "LS_2r"
LS2_final = "ls2_final"
landscape2 = "landscape2"

outTimes4 = arcpy.sa.Times(ls2, 100)
outTimes4.save(LS_2)
OutInt = Int(LS_2)
OutInt.save(LS_2r)
rescale0to1(LS_2r, value, home)
lslook2 = Lookup(LS_2r, "NORMALIZED")
lslook2.save(LS2_final)
outTimes6 = Times(LS2_final, 100)
outTimes6.save(landscape2)

print arcpy.GetRasterProperties_management(landscape2, "MINIMUM")
print arcpy.GetRasterProperties_management(landscape2, "MAXIMUM")

arcpy.CheckInExtension("spatial")  # check in the spatial analyst licence you checked out

# Export the two final landscapes from  ESRI GRID file to TIF
FinalFolder = home + "\FinalFolder"  # make sure this folder exists in the home folder
arcpy.RasterToOtherFormat_conversion("landscape1;landscape2","FinalFolder", "TIFF")

# Clean up
del_list = ("gdd_surv", "lcc_surv","lev_surv", "lev_temp1","ls1","LS_1","LS_1r",
            "LS1_final", "ls2", "LS_2", "LS_2r", "LS2_final")

for raster in del_list:
    arcpy.Delete_management(raster)

