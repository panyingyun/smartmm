#ifndef _MESHLIB_PLANE_3D_H_
#define _MESHLIB_PLANE_3D_H_

/*!
*      \file Plane3D.h
*      \brief a plane is R3
*	   \author David Gu
*      \date 10/14/2020
*
*/
#include <vector>
#include "Segment3D.h"
#include "Ray3D.h"
#include "Point.h"
#include "Polygon2D.h"
#include "Polygon3D.h"

namespace MeshLib {

	/*!
	*	\brief 3D plane
	*
	*/
	class CPlane3D
	{
	public:

		/*!
		*	CPlane Constructor
		*/
		CPlane3D() {};
		~CPlane3D() {};

		//construct local orthonomal frame
		CPlane3D(CPoint a, CPoint b, CPoint c)
		{
			m_o = a;
			CPoint n = (b - a) ^ (c - a);
			n /= n.norm();
			m_height = m_o * n;
			m_e[0] = (b - a) / (b - a).norm();
			m_e[2] = n;
			m_e[1] = m_e[2] ^ m_e[0];
		};

		//transform a 3D point to the local coordiantes of the plane
		CPoint transform(CPoint& p)
		{
			CPoint d = p - m_o;
			return CPoint(d * m_e[0], d * m_e[1], d * m_e[2]);
		};

		/*
		*	if a point is above the plane, 
		*   +1 above
		*	-1 below
		*	0  on
		*/
		
		int  onTop(CPoint& p)
		{
			double d = (p - m_o)*m_e[2];
			if (d > 0) return +1;
			if (d < 0) return -1;
			return 0;
		}

		/*!
		*	Intersection between a plane and line segment
		*   if it intersects, then return true; q is the intersection point
		*/
		bool intersection(CSegment3D& seg, CPoint & q)
		{
			CPoint s = seg.start();
			CPoint e = seg.end();

			int ds = onTop(s);
			int de = onTop(e);
			if (ds * de != -1) return false;

			double t = (m_height - (m_e[2]*e))/((s - e) * m_e[2]);
			q = s * t + e * (1 - t);
			return true;
		}

		/*!
		*	given (x,y), find the z-value on the plane
		*/
		double z(CPoint2 p)
		{
			CPoint n = m_e[2];
			return (m_height - n[0] * p[0] - n[1] * p[1]) / n[2];
		}

		/*!
		*	inverse_transformation
		*/
		CPoint inverse_transform(CPoint p)
		{
			return m_o + (m_e[0] * p[0]) + (m_e[1] * p[1]) + (m_e[2] * p[2]);
		}
	protected:
		CPoint  m_e[3];
		CPoint  m_o;
		double  m_height;
	};

};
#endif


