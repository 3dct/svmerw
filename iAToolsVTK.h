/************************************  SVMERW  ****************************************
* **  Command line tool for segmenting multi-modal 3D data with a Support Vector   ** *
* **  Machine (SVM) and the Extended Random Walker (ERW)                           ** *
***************************************************************************************
* Copyright (C) 2017 Bernhard Fr�hler                                                 *
***************************************************************************************
* This program is free software: you can redistribute it and/or modify it under the   *
* terms of the GNU General Public License as published by the Free Software           *
* Foundation, either version 3 of the License, or (at your option) any later version. *
*                                                                                     *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY     *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     *
* PARTICULAR PURPOSE.  See the GNU General Public License for more details.           *
*                                                                                     *
* You should have received a copy of the GNU General Public License along with this   *
* program.  If not, see http://www.gnu.org/licenses/                                  *
***************************************************************************************
* Contact: FH O� Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*       Stelzhamerstra�e 23, 4600 Wels / Austria, Email: bernhard.froehler@fh-wels.at *
**************************************************************************************/
#pragma once

#include <vtkSmartPointer.h>

class vtkImageData;

class QString;

// image creation:
void DeepCopy(vtkSmartPointer<vtkImageData> input, vtkSmartPointer<vtkImageData> output);
vtkSmartPointer<vtkImageData> AllocateImage(vtkSmartPointer<vtkImageData> img);
vtkSmartPointer<vtkImageData> AllocateImage(int vtkType, int dimensions[3], double spacing[3]);

// image I/O (using ITK methods of iAITKIO)
void StoreImage(vtkSmartPointer<vtkImageData> image, QString const & filename, bool useCompression);
vtkSmartPointer<vtkImageData> ReadImage(QString const & filename);
