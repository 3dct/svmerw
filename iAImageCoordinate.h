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

typedef int iAVoxelIndexType;

/**
 * Helper struct for containing 3D image coordinates
 */
struct iAImageCoordinate
{
public:
	enum IndexOrdering
	{
		RowColDepMajor, // x, y(, z)
		ColRowDepMajor, // y, x(, z)
	};
	iAImageCoordinate();
	iAImageCoordinate(iAVoxelIndexType x, iAVoxelIndexType y, iAVoxelIndexType z);
	iAVoxelIndexType x, y, z;
};

bool operator==(iAImageCoordinate const & a, iAImageCoordinate const & b);

/** conversion utilty class
 * to convert from a up to 3-dimensional image index (x, y, z) to a "flat" index
 */

class iAImageCoordConverter
{
public:
	iAImageCoordConverter(iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth=1,
		iAImageCoordinate::IndexOrdering ordering=iAImageCoordinate::RowColDepMajor
	);
	iAImageCoordinate GetCoordinatesFromIndex(iAVoxelIndexType index) const;
	iAVoxelIndexType GetIndexFromCoordinates(iAImageCoordinate coords) const;
	iAVoxelIndexType GetVertexCount() const;
	static iAImageCoordinate GetCoordinatesFromIndex(
		iAVoxelIndexType index,
		iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth,
		iAImageCoordinate::IndexOrdering ordering);
	static iAVoxelIndexType GetIndexFromCoordinates(
		iAImageCoordinate coords,
		iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth,
		iAImageCoordinate::IndexOrdering ordering);
	iAVoxelIndexType GetWidth() const;
	iAVoxelIndexType GetHeight() const;
	iAVoxelIndexType GetDepth() const;
private:
	iAVoxelIndexType m_width;
	iAVoxelIndexType m_height;
	iAVoxelIndexType m_depth;
	iAImageCoordinate::IndexOrdering m_ordering;
};
