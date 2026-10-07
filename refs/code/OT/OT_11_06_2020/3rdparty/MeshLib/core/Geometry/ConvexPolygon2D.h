#ifndef _MESHLIB_CONVEX_POLYGON_2D_H_
#define _MESHLIB_CONVEX_POLYGON_2D_H_

/*!
*      \file ConvexPolygon2D.h
*      \brief convex planar polygon
*	   \author David Gu
*      \date 10/15/2020
*
*/
#include <vector>
#include "Segment2D.h"
#include "Ray2D.h"

namespace MeshLib{

	class ConvexPolygon2D
	{
	public:
		std::vector<CPoint2> Corners;

	public:
		ConvexPolygon2D(std::vector<CPoint2> corners)
		{
			Corners = corners;
		}
	};

}

#endif
