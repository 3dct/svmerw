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
#include "iASpectrumType.h"

iASpectralVoxelData::iASpectralVoxelData(): m_maxSum(-1)
{}

iASpectralVoxelData::~iASpectralVoxelData()
{}



iASpectrumDataType iASpectrumType::operator[](size_t channelIdx) const
{
	return get(channelIdx);
}

QSharedPointer<iASpectrumType const> iASpectrumType::normalized() const
{
	QSharedPointer<iAStandaloneSpectrumType> result(new iAStandaloneSpectrumType(size()));
	iASpectrumDataType sum = 0;
	for(iASpectrumType::IndexType i = 0; i<size(); ++i)
	{
		sum += get(i);
	}
	for(iASpectrumType::IndexType i = 0; i<size(); ++i)
	{
		result->set(i, get(i) / sum);
	}
	return result;
}



iADirectAccessSpectrumType::iADirectAccessSpectrumType(iASpectralVoxelData const & data, size_t voxelIdx):
	m_data(data),
	m_voxelIdx(voxelIdx)
{}

iASpectrumDataType iADirectAccessSpectrumType::get(size_t channelIdx) const
{
	iASpectrumDataType value = m_data.get(m_voxelIdx, channelIdx);
	return value;
}

iASpectrumType::IndexType iADirectAccessSpectrumType::size() const
{
	return m_data.channelCount();
}



iAStandaloneSpectrumType::iAStandaloneSpectrumType(IndexType size):
	m_data(size)
{}

iASpectrumDataType iAStandaloneSpectrumType::get(size_t idx) const
{
	return m_data[idx];
}

iASpectrumType::IndexType iAStandaloneSpectrumType::size() const
{
	return m_data.size();
}

void iAStandaloneSpectrumType::set(iASpectrumType::IndexType idx, iASpectrumDataType value)
{
	m_data[idx] = value;
}
