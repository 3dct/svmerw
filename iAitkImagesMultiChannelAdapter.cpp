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
#include "iAitkImagesMultiChannelAdapter.h"

iAvtkImagesMultiChannelAdapter::iAvtkImagesMultiChannelAdapter(size_t width, size_t height, size_t depth):
	m_coordConv(width, height, depth)
{
}

	
void iAvtkImagesMultiChannelAdapter::AddImage(vtkSmartPointer<vtkImageData> img)
{
	int extent[6];
	img->GetExtent(extent);
	assert((extent[1]-extent[0]+1) == m_coordConv.GetWidth() &&
		(extent[3]-extent[2]+1) == m_coordConv.GetHeight() &&
		(extent[5]-extent[4]+1) == m_coordConv.GetDepth());
	m_images.push_back(img);
}

size_t iAvtkImagesMultiChannelAdapter::size() const
{
	return m_coordConv.GetVertexCount();
}

size_t iAvtkImagesMultiChannelAdapter::channelCount() const
{
	return m_images.size();
}

QSharedPointer<iASpectrumType const> iAvtkImagesMultiChannelAdapter::get(size_t voxelIdx) const
{
	return QSharedPointer<iASpectrumType const>(new iADirectAccessSpectrumType(*this, voxelIdx));
}

iASpectrumDataType iAvtkImagesMultiChannelAdapter::get(size_t voxelIdx, size_t channelIdx) const
{
	iAImageCoordinate coords = m_coordConv.GetCoordinatesFromIndex(voxelIdx);
	iASpectrumDataType value = m_images[channelIdx]->GetScalarComponentAsDouble(coords.x, coords.y, coords.z, 0);
	return value;
}
