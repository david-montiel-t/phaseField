#!/usr/bin/env python3

import sys
from visit import *

#DeleteAllPlots()

# Step 1: Open a database (the whole .vtu time series)
dbname="solution-*.vtu database"
OpenDatabase(dbname)

no_ops=21

#Step 2: Add contour plot to obtain grid information and dimensionality
#Add Contour plot
for j in range(0,no_ops):
  op_name="n"+str(j)
  AddPlot("Pseudocolor", op_name)
  AddOperator("Threshold")
  ThresholdAtts = ThresholdAttributes()
  ThresholdAtts.lowerBounds = (0.5)
  #ThresholdAtts.upperBounds = (j+0.5)
  SetOperatorOptions(ThresholdAtts)
  DrawPlots()

  print("Saving Order parameters %d" % (j))
  SaveWindowAtts = SaveWindowAttributes()
  SaveWindowAtts.family = 0
  SaveWindowAtts.fileName = "n"+str(j)
  SetSaveWindowAttributes(SaveWindowAtts)
  SaveWindow()
  DeleteAllPlots()
CloseDatabase(dbname)

sys.exit()
