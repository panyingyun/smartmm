#ifndef _MESHLIB_SPHERICAL_POLYGON_3D_H_
#define _MESHLIB_SPHERICAL_POLYGON_3D_H_

/*!
*      \file SphericalPolygon3D.h
*      \brief Spherical polygon
*	   \author David Gu
*      \date 05/03/2020
*
*/
#include <vector>
#include "Segment3D.h"
#include "Ray3D.h"
#include "Point.h"
#include "Polygon2D.h"
#include "Polygon3D.h"
#include "SphericalTriangle.h"

namespace MeshLib {


	/*!
	*	\brief spacial polygon
	*
	*/
	class CSphericalPolygon
	{
	public:

		/*!
		*	CPolygon constructor
		*/
		CSphericalPolygon() {};
		~CSphericalPolygon() {};

		/*! add one segment in sequencial order */
		void add(CSegment3D& seg)
		{
			m_edges.push_back(seg);
		};

		/*! spherical polygon area */
		double area();

		/*! clip by the unit disk*/
		int finite_clip(CSphericalPolygon& inner_poly) { return 0; };

		/*! the edges*/
		std::vector<CSegment3D>& edges() { return m_edges; };

		/*! centeter of mass */
		CPoint mass_center();


	protected:
		std::vector<CSegment3D> m_edges;
	};

	/*! the total area of a planar polygon */

inline	double CSphericalPolygon::area()
	{
		double s = 0;
		std::vector<CPoint> pts;
		for (size_t i = 0; i < m_edges.size(); i++)
		{
			pts.push_back(m_edges[i].start());
		}

		for (size_t i = 1; i < m_edges.size() - 1; i++)
		{
			CSphericalTriangle tri(pts[0], pts[i], pts[i + 1]);
			s += tri.area();
		}
		return s;
	};

inline	CPoint CSphericalPolygon::mass_center()
	{
		std::vector<CPoint> pts;
		for (size_t i = 0; i < m_edges.size(); i++)
		{
			CSegment3D& s = m_edges[i];
			pts.push_back(s.start());
		}

		CPoint  center(0, 0, 0);
		double  area = 0;

		for (int j = 1; j < pts.size() - 1; j++)
		{
			CPoint c = (pts[0] + pts[j] + pts[j + 1]) / 3;
			CPoint e = pts[j + 0] - pts[0];
			CPoint d = pts[j + 1] - pts[0];
			CPoint n = e ^ d;
			double s = n.norm() / 2;
			area += s;
			center += c * s;
		}

		center /= area;

		return center;
	}

};
#endif


