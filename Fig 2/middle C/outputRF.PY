from abaqus import *
from abaqusConstants import *
from viewerModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
odb = session.odbs['C:/Users/Wildcat/Documents/2023 sCBM/UEL/UEL-main/Lagrange Multiplier/40 dof/TPB_notch_M7_good/PNAS/medium/sample.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    'RF', NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=("PART-1-1.DL", ))
xyp = session.XYPlot('XYPlot-2')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
xy1 = session.xyDataObjects['_RF:RF2 PI: PART-1-1 N: 20']
xy2 = session.xyDataObjects['_RF:RF2 PI: PART-1-1 N: 21']
xy3 = session.xyDataObjects['_RF:RF2 PI: PART-1-1 N: 22']
xy4 = session.xyDataObjects['_RF:RF2 PI: PART-1-1 N: 293']
xy5 = session.xyDataObjects['_RF:RF2 PI: PART-1-1 N: 294']
xy6 = -sum((xy1, xy2, xy3, xy4, xy5))
xyp = session.xyPlots['XYPlot-2']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
c1 = session.Curve(xyData=xy6)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)

# Find the maximum value
max_value = np.max(xy6)

# Print or use the maximum value
print("Maximum value in xy6:", max_value)