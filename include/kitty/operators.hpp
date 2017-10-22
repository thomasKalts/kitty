/* kitty: C++ truth table library
 * Copyright (C) 2017  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file operators.hpp
  \brief Implements operator shortcuts to operations

  \author Mathias Soeken
*/

#pragma once

#include "operations.hpp"

namespace kitty
{

/*! \brief Operator for unary_not */
template<typename TT>
inline TT operator~( const TT& tt )
{
  return unary_not( tt );
}

/*! \brief Operator for binary_and */
template<typename TT>
inline TT operator&( const TT& first, const TT& second )
{
  return binary_and( first, second );
}

/*! \brief Operator for binary_and and assign */
template<typename TT>
inline void operator&=( TT& first, const TT& second )
{
  first = binary_and( first, second );
}

/*! \brief Operator for binary_or */
template<typename TT>
inline TT operator|( const TT& first, const TT& second )
{
  return binary_or( first, second );
}

/*! \brief Operator for binary_or and assign */
template<typename TT>
inline void operator|=( TT& first, const TT& second )
{
  first = binary_or( first, second );
}

/*! \brief Operator for binary_xor */
template<typename TT>
inline TT operator^( const TT& first, const TT& second )
{
  return binary_xor( first, second );
}

/*! \brief Operator for binary_xor and assign */
template<typename TT>
inline void operator^=( TT& first, const TT& second )
{
  first = binary_xor( first, second );
}

/*! \brief Operator for equal */
template<typename TT>
inline bool operator==( const TT& first, const TT& second )
{
  return equal( first, second );
}

/*! \brief Operator for less_than */
template<typename TT>
inline bool operator<( const TT& first, const TT& second )
{
  return less_than( first, second );
}

} // namespace kitty