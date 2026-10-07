#ifndef _MESHLIB_SPHERICAL_TRIANGLE_H_
#define _MESHLIB_SPHERICAL_TRIANGLE_H_

/*!
*      \file SphericalTriangle.h
*      \brief Spherical Triangle
*	   \author David Gu
*      \date 05/04/2020
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
	*	\brief spacial polygon
	*
	*/
	class CSphericalTriangle
	{
	public:

		/*!
		*	CPolygon constructor
		*/
		CSphericalTriangle(CPoint& A, CPoint& B, CPoint& C);
		CSphericalTriangle(CPoint pt[]);
		~CSphericalTriangle() {};

		/*! spherical polygon area */
		double area();

		/*! centeter of mass */
		CPoint mass_center();


	protected:
		CPoint m_p[3];
		double m_length[3];
		double m_angle[3];
	};

	/*! the total area of a planar polygon */

inline	double CSphericalTriangle::area()
	{
		//return m_angle[0] + m_angle[1] + m_angle[2] - PI;
		return m_angle[0] + m_angle[1] + m_angle[2] - 3.141592653589793;
	};

inline	CPoint CSphericalTriangle::mass_center()
	{
		CPoint  center(0, 0, 0);
		return center;
	}

inline	CSphericalTriangle::CSphericalTriangle(CPoint pt[])
	{
		m_p[0] = pt[0]; m_p[1] = pt[1]; m_p[2] = pt[2];
		double cs[3], sn[3];
		//calculate the lengths and angles
		for (int i = 0; i < 3; i++)
		{
			cs[i] = m_p[(i + 1) % 3] * m_p[(i + 2) % 3];
			sn[i] = sqrt(1 - cs[i] * cs[i]);
			m_length[i] = acos(cs[i]);
		}
		for (int i = 0; i < 3; i++)
		{
			double c = m_length[(i + 0) % 3];
			double a = m_length[(i + 1) % 3];
			double b = m_length[(i + 2) % 3];
			double C = (cos(a) * cos(b) - cos(c)) / (sin(a) * sin(b));
			m_angle[i] = acos(C);
		}
	}

inline	CSphericalTriangle::CSphericalTriangle(CPoint &A, CPoint &B, CPoint &C)
	{
		m_p[0] = A; m_p[1] = B; m_p[2] = C;

		double cs[3], sn[3];
		//calculate the lengths and angles
		for (int i = 0; i < 3; i++)
		{
			cs[i] = m_p[(i + 1) % 3] * m_p[(i + 2) % 3];
			sn[i] = sqrt(1 - cs[i] * cs[i]);
			m_length[i] = acos(cs[i]);
		}
		for (int i = 0; i < 3; i++)
		{
			double c = m_length[(i + 0) % 3];
			double a = m_length[(i + 1) % 3];
			double b = m_length[(i + 2) % 3];
			double C = (cos(c)-cos(a) * cos(b)) / (sin(a) * sin(b));
			m_angle[i] = acos(C);
		}
	}

};
#endif


