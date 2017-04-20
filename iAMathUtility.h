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

#include <cassert>

/**
 * make sure the given value is inside the given interval
 * @param min the minimum value which should be returned
 * @param max the maximum value which should be returned
 * @param val the value to check
 * @return min if the given value smaller than min, max if
 *         given value bigger than max, value if it is in
 *         between min and max
 */
template <typename T>
T clamp(T min, T max, T val)
{
	return (val < min) ? min : ((val > max) ? max : val);
}

/**
 * map value from given interval to "norm" interval [0..1]
 * @param minSrcVal minimum value of source interval
 * @param maxSrcVal maximum value of source interval
 * @param value a value in source interval
 * @return if norm was in interval [minSrcVal..maxSrcVal], the
 *     corresponding mapped value in interval [0..1]
 */
template <typename SrcType>
double mapToNorm(SrcType minSrcVal, SrcType maxSrcVal, SrcType value)
{
	assert (value >= minSrcVal && value <= maxSrcVal);
	if (maxSrcVal == minSrcVal)
	{	// to prevent division by 0
		return 0;
	}
	double returnVal = static_cast<double>(value - minSrcVal) / (maxSrcVal - minSrcVal);
	assert(returnVal >= 0 && returnVal <= 1);
	return returnVal;
}

/**
 * map value from "norm" interval [0..1] to the given interval
 * @param minDstVal minimum value of destination interval
 * @param maxDstVal maximum value of destination interval
 * @param norm a value in interval [0..1]
 * @return if norm was in [0..1], the corresponding mapped value
 *     in interval [minDstVal..maxDstVal]
 */
template <typename DstType>
DstType mapNormTo(DstType minDstVal, DstType maxDstVal, double norm)
{
	double returnVal = norm * (maxDstVal - minDstVal) + minDstVal;
	assert (returnVal >= minDstVal && returnVal <= maxDstVal);
	return returnVal;
}

/**
  * map value from one interval to another
  * @param minSrcVal minimum value of source interval
  * @param maxSrcVal maximum value of source interval
  * @param minDstVal minimum value of destination interval
  * @param maxDstVal maximum value of destination interval
  * @param value a value in source interval
  * @return if value was in interval [minSrcVal..maxSrcVal], the
  *     corresponding mapped value in interval [minDstVal..maxDstVal]
  */
template <typename SrcType, typename DstType>
DstType mapValue(SrcType minSrcVal, SrcType maxSrcVal, DstType minDstVal, DstType maxDstVal, SrcType value)
{
	return mapNormTo(minDstVal, maxDstVal, mapToNorm(minSrcVal, maxSrcVal, value));
}

/**
  * round a number to the nearest integer representation (by "round half away from zero" method)
  * @param number the number to round
  * @return the rounded number
  */
template <typename T>
T round(T const & number)
{
	return number < 0.0 ? std::ceil(number - 0.5) : std::floor(number + 0.5);
}
