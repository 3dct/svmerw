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

#include <cassert>

/**
 * Class representing a generic (single-parameter) function as can be passed into
 * the functional boxplot calculation
 */
template <typename ArgType, typename ValType>
class iAFunction
{
public:
	iAFunction();
	iAFunction(ArgType size);
	iAFunction(iAFunction<ArgType, ValType> const & other);
	~iAFunction();
	typedef ArgType ArgumentType;
	typedef ValType ValueType;
	
	ValType get(ArgType a) const;
	void set(ArgType a, ValType);
	ValType operator[](ArgType a) const;
	//ValType & operator[](ArgType a);
	ArgType size() const;
	void init(ArgType size);

private:
	ArgumentType m_size;
	ValType * m_data;
};


template <typename ArgType, typename ValType>
iAFunction<ArgType, ValType>::iAFunction():
	m_size(0),
	m_data(0)
{}

template <typename ArgType, typename ValType>
iAFunction<ArgType, ValType>::iAFunction(ArgType size)
{
	init(size);
}

template <typename ArgType, typename ValType>
iAFunction<ArgType, ValType>::iAFunction(iAFunction<ArgType, ValType> const & other)
{
	init(other.m_size);
	for (int i=0; i<other.m_size; ++i)
	{
		m_data[i] = other.m_data[i];
	}
}

template <typename ArgType, typename ValType>
void iAFunction<ArgType, ValType>::init(ArgType size)
{
	m_size = size;
	m_data = new ValType[m_size];
}

template <typename ArgType, typename ValType>
iAFunction<ArgType, ValType>::~iAFunction()
{
	delete [] m_data;
}

template <typename ArgType, typename ValType>
ValType iAFunction<ArgType, ValType>::get(ArgType idx) const
{
	assert (m_data);
	assert (idx >= 0 && idx < m_size);
	return m_data[idx];
}

template <typename ArgType, typename ValType>
void iAFunction<ArgType, ValType>::set(ArgType idx, ValType v)
{
	assert (m_data);
	assert (idx >= 0 && idx < m_size);
	m_data[idx] = v;
}

template <typename ArgType, typename ValType>
ValType iAFunction<ArgType, ValType>::operator[](ArgType idx) const
{
	assert (m_data);
	assert (idx >= 0 && idx < m_size);
	return m_data[idx];
}

/*
template <typename ArgType, typename ValType>
ValType & iAFunction<ArgType, ValType>::operator[](ArgType idx)
{
	assert (m_data);
	assert (idx >= 0 && idx < m_size);
	return m_data[idx];
}
*/

template <typename ArgType, typename ValType>
ArgType iAFunction<ArgType, ValType>::size() const
{
	return m_size;
}
