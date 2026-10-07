#ifndef _MESHLIB_SEGMENT_2D_H_
#define _MESHLIB_SEGMENT_2D_H_

/*!
 *      \file Segment 2D.h
 *      \brief Planar segment
 *	   \author David Gu
 *      \date 05/01/2020
 *
 */
#include "Circle.h"
#include "Point2.h"
#include "predicates.h"
#include <vector>

namespace MeshLib
{

/*!
 *	\brief CSegment2D class, a finite Segment on the plane
 *
 */
class CSegment2D
{
  public:
    /*!
     *	CSegment2D Constructor,
     *   \param s starting point
     *   \param e ending   point
     */

    CSegment2D(CPoint2 st, CPoint2 ed)
    {
        m_start = st;
        m_end = ed;
    };
    CSegment2D(const CSegment2D& segment)
    {
        m_start = segment.m_start;
        m_end = segment.m_end;
    };
    CSegment2D(){};

    ~CSegment2D(){};

    bool intersect(CCircle& circle, std::vector<CPoint2>& intersection_points);
    bool intersect(CSegment2D& sg, CPoint2& q);

    /* start point */
    CPoint2& start() { return m_start; };
    CPoint2& end() { return m_end; };

    /*! verify which side the point is on the segment
     *
     * +1 left
     * 0  on
     * -1 right
     */
    int side(CPoint2 pt)
    {
        double p[2], s[2], e[2];

        p[0] = pt[0];
        p[1] = pt[1];
        s[0] = m_start[0];
        s[1] = m_start[1];
        e[0] = m_end[0];
        e[1] = m_end[1];

        double d = orient2d(p, s, e);
        if (d > 0)
            return +1;
        if (d < 0)
            return -1;

        return 0;
    };

  protected:
    CPoint2 _intersection(CPoint2& cp1, CPoint2& cp2, CPoint2& s, CPoint2& e);

  protected:
    CPoint2 m_start;
    CPoint2 m_end;
};

inline bool CSegment2D::intersect(CCircle& circ, std::vector<CPoint2>& intersection_point)
{
    double A = (m_start - m_end).norm2();
    double B = 2 * ((m_start - m_end) * (m_end - circ.c()));
    double C = (m_end - circ.c()).norm2() - circ.r() * circ.r();

    double D = B * B - 4 * A * C;
    if (D < 0)
        return false;

    double t = (-B - sqrt(D)) / (2 * A);

    if (t >= 0 && t <= 1)
    {
        CPoint2 p = m_start * t;
        CPoint2 q = m_end * (1 - t);

        intersection_point.push_back(p + q);
    }
    t = (-B + sqrt(D)) / (2 * A);
    if (t >= 0 && t <= 1)
    {
        intersection_point.push_back((m_start * t) + (m_end * (1 - t)));
    }

    return intersection_point.size() > 0;
};

inline bool CSegment2D::intersect(CSegment2D& sg, CPoint2& intersection_point)
{
    CPoint2 s1 = sg.start();
    CPoint2 e1 = sg.end();

    CPoint2 s2 = m_start;
    CPoint2 e2 = m_end;

    // verify if the poings s1 and e1 are at different sides of this segment
    if (side(s1) * side(e1) != -1)
        return false;
    // verify if the poings s2 and e2 are at different sides of sg
    if (sg.side(s2) * sg.side(e2) != -1)
        return false;

    intersection_point = _intersection(s1, e1, s2, e2);

    return true;
};

// calculate intersection point
CPoint2 CSegment2D::_intersection(CPoint2& cp1, CPoint2& cp2, CPoint2& s, CPoint2& e)
{
    CPoint2 dc = {cp1[0] - cp2[0], cp1[1] - cp2[1]};
    CPoint2 dp = {s[0] - e[0], s[1] - e[1]};

    double n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0];
    double n2 = s[0] * e[1] - s[1] * e[0];
    double n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0]);

    return {(n1 * dp[0] - n2 * dc[0]) * n3, (n1 * dp[1] - n2 * dc[1]) * n3};
}

}; // namespace MeshLib
#endif
