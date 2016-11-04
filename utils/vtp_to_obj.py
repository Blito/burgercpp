#!/usr/bin/env python

import vtk
import os

reader = vtk.vtkXMLPolyDataReader()
path = "C:/Simeco/simeco/test/data/ircad11_large/body/body_0_triangle.vtp"
reader.SetFileName(path)
reader.Update()

polydata = reader.GetOutput()
#print(polydata)

cell_array = polydata.GetPolys()
idList = vtk.vtkIdList()

while (cell_array.GetNextCell(idList)):
	for cell_id in range(idList.GetNumberOfIds()):
		point_id = idList.GetId(cell_id)
		point = polydata.GetPoint(point_id)
		print point[0], point[1], point[2],
	print


#interactor.Start()