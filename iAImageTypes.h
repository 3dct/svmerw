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
#pragma once

#include <itkImage.h>

#include <vtkSmartPointer.h>

class vtkImageData;

// TODO: data type -> short/char
typedef int LabelPixelType;
typedef itk::Image<LabelPixelType, 3> LabelImageType;
typedef LabelImageType::Pointer LabelImagePointer;

typedef itk::Image<double, 3> ProbabilityImageType;
typedef ProbabilityImageType::Pointer ProbabilityImagePointer;

typedef vtkSmartPointer<vtkImageData> PriorModelImagePointer;

typedef int PreviewPixelType;
typedef itk::Image<PreviewPixelType, 3> PreviewImageType;
typedef PreviewImageType::Pointer PreviewImagePointer;

static const int m_DIM = 3;
typedef itk::ImageBase< m_DIM > ImageBaseType;
typedef ImageBaseType::Pointer ImagePointer;
typedef ImagePointer ClusterImageType;
