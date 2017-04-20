/************************************  SVMERW  ****************************************
* **  Command line tool for segmenting multi-modal 3D data with a Support Vector   ** *
* **  Machine (SVM) and the Extended Random Walker (ERW)                           ** *
***************************************************************************************
* Copyright (C) 2017 Bernhard Fröhler                                                 *
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
* Contact: FH OÖ Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*       Stelzhamerstraße 23, 4600 Wels / Austria, Email: bernhard.froehler@fh-wels.at *
**************************************************************************************/
#include "iAToolsVTK.h"

#include "iAITKIO.h"

#include <vtkImageData.h>

#include "iAConnector.h"

void DeepCopy(vtkSmartPointer<vtkImageData> input, vtkSmartPointer<vtkImageData> output)
{
	output->DeepCopy(input);
}

vtkSmartPointer<vtkImageData> AllocateImage(int vtkType, int dimensions[3], double spacing[3])
{
	vtkSmartPointer<vtkImageData> result = vtkSmartPointer<vtkImageData>::New();
	result->SetDimensions(dimensions);
	result->AllocateScalars(vtkType, 1);
	result->SetSpacing(spacing);
	return result;
}

vtkSmartPointer<vtkImageData> AllocateImage(vtkSmartPointer<vtkImageData> img)
{
	return AllocateImage(img->GetScalarType(), img->GetDimensions(), img->GetSpacing());
}


void StoreImage(vtkSmartPointer<vtkImageData> image, QString const & filename, bool useCompression)
{
	iAConnector con;
	con.SetImage(image);
	iAITKIO::ScalarPixelType pixelType = con.GetITKScalarPixelType();
	iAITKIO::writeFile(filename, con.GetITKImage(), pixelType, useCompression);
}

vtkSmartPointer<vtkImageData> ReadImage(QString const & filename)
{
	iAConnector con;
	iAITKIO::ScalarPixelType pixelType;
	iAITKIO::ImagePointer img = iAITKIO::readFile(filename, pixelType);
	con.SetImage(img);
	return con.GetVTKImage();
}
