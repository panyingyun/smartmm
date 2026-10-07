#ifndef _MESHLIB_POLYGON_2D_H_
#define _MESHLIB_POLYGON_2D_H_

/*!
 *      \file Polygon2D.h
 *      \brief planar polygon
 *	   \author David Gu
 *      \date 05/01/2020
 *
 */
#include "Ray2D.h"
#include "Segment2D.h"
#include "predicates.h"
#include <deque>
#include <vector>

namespace MeshLib
{

/*!
 *	\brief planar polygon
 *
 */
class CPolygon2D
{
  public:
    /*!
     *	CPolygon constructor
     */
    CPolygon2D(){};
    ~CPolygon2D(){};

    /*! add one segment in sequencial order */
    void add(CSegment2D& seg) { m_edges.push_back(seg); };

    /*! polygon area */
    double area();

    /*! clip by the unit disk*/
    int finite_clip(CPolygon2D& inner_poly);

    /*! the edges*/
    std::vector<CSegment2D>& edges() { return m_edges; };

    /*! the polygon is an infinite cell, clipped by the unit circle */
    int infinite_clip(CPolygon2D& inner_poly);

    /*! centeter of mass */
    CPoint2 mass_center();

    /*! verify if a point is inside the polygon
     *   assume the polygon is convex
     *
     */
    bool insideConvexPolygon(CPoint2& pt)
    {

        for (size_t i = 0; i < m_edges.size(); i++)
        {
            int side = m_edges[i].side(pt);
            // +1 left, -1 right, 0 on
            if (side != 1)
                return false;
        }

        return true;
    }

    /*!
     *	intersection between a ray and the convex polygon
     */
    bool intersectConvexPolygon(CRay2D& ray, CPoint2& p);

    /*!
     *	intersection between the convex polygon and a segment
     *  Input segment  : seg
     *  Output segment : intersection
     */
    bool intersectConvexPolygon(CSegment2D& seg, CSegment2D& intersection);

  protected:
    std::vector<CSegment2D> m_edges;
};

/*! the total area of a planar polygon */

inline double CPolygon2D::area()
{
    double s = 0;
    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment2D& seg = m_edges[i];
        s += seg.start() ^ seg.end();
    }
    return s / 2.0;
};

/*! the finite polygon is clipped by the unit disk */
inline int CPolygon2D::finite_clip(CPolygon2D& inner_poly)
{
    double epsilon = 1e-30;

    CCircle circ(CPoint2(0, 0), 1);
    std::vector<int> on_circle_idx;

    std::vector<CPoint2> inner_pts;
    std::vector<CPoint2> outer_pts;

    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment2D& seg = m_edges[i];
        CPoint2& pt = seg.start();

        double R = pt.norm2();

        if (R < 1 - epsilon)
        {
            inner_pts.push_back(pt);
        }
        if (R > 1 + epsilon)
        {
            outer_pts.push_back(pt);
        }
        if (R <= 1 + epsilon && R >= 1 - epsilon)
        {
            inner_pts.push_back(pt);
            outer_pts.push_back(pt);
        }
        std::vector<CPoint2> intersection_points;
        bool cross = seg.intersect(circ, intersection_points);
        if (!cross)
            continue;

        for (size_t j = 0; j < intersection_points.size(); j++)
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points[j]);
            outer_pts.push_back(intersection_points[j]);
        }
    }

    // completely outside
    if (inner_pts.empty())
        return -1;

    // completely inside
    if (outer_pts.empty())
        return +1;

    int n = (int) inner_pts.size();
    int on_circle_indx = 0;

    for (int i = 0; i < (int) on_circle_idx.size(); i++)
    {
        if ((on_circle_idx[i] + 1) % n == on_circle_idx[(i + 1) % on_circle_idx.size()])
        {
            on_circle_indx = on_circle_idx[i];
        }
    }

    for (int i = 0; i < n; i++)
    {
        CSegment2D seg(inner_pts[(i + on_circle_indx) % n], inner_pts[(i + 1 + on_circle_indx) % n]);
        inner_poly.add(seg);
    }

    return 0;
}

/*! the finite polygon is clipped by the unit disk */
inline int CPolygon2D::infinite_clip(CPolygon2D& inner_poly)
{
    double epsilon = 1e-30;

    CCircle circ(CPoint2(0, 0), 1);
    std::vector<int> on_circle_idx;

    std::vector<CPoint2> inner_pts;
    std::vector<CPoint2> outer_pts;

    {
        CSegment2D& sg = m_edges.front();
        CPoint2 st = sg.start();
        CPoint2 dr = sg.end() - sg.start();
        CRay2D ray(st, dr);
        std::vector<CPoint2> intersection_points;
        bool cross = ray.intersect(circ, intersection_points);

        if (!intersection_points.empty())
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points.front());
            outer_pts.push_back(intersection_points.front());
        }
    }

    for (size_t i = 1; i < m_edges.size(); i++)
    {
        CSegment2D& seg = m_edges[i];
        CPoint2& pt = seg.start();

        double R = pt.norm2();

        if (R < 1 - epsilon)
        {
            inner_pts.push_back(pt);
        }
        if (R > 1 + epsilon)
        {
            outer_pts.push_back(pt);
        }
        if (R <= 1 + epsilon && R >= 1 - epsilon)
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(pt);
            outer_pts.push_back(pt);
        }

        if (i == m_edges.size() - 1)
            continue;

        std::vector<CPoint2> intersection_points;
        bool cross = seg.intersect(circ, intersection_points);
        if (!cross)
            continue;

        for (size_t j = 0; j < intersection_points.size(); j++)
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points[j]);
            outer_pts.push_back(intersection_points[j]);
        }
    }

    {
        CSegment2D& sg = m_edges.back();
        CPoint2 st = sg.start();
        CPoint2 dr = sg.end() - sg.start();
        CRay2D ray(st, dr);
        std::vector<CPoint2> intersection_points;
        bool cross = ray.intersect(circ, intersection_points);

        if (!intersection_points.empty())
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points.front());
            outer_pts.push_back(intersection_points.front());
        }
    }

    // completely outside
    if (inner_pts.empty())
        return -1;

    // completely inside
    if (outer_pts.empty())
        return +1;

    int n = (int) inner_pts.size();
    int on_circle_indx = 0;

    for (int i = 0; i < (int) on_circle_idx.size(); i++)
    {
        if ((on_circle_idx[i] + 1) % n == on_circle_idx[(i + 1) % on_circle_idx.size()])
        {
            on_circle_indx = on_circle_idx[i];
        }
    }

    for (int i = 0; i < n; i++)
    {
        CSegment2D seg(inner_pts[(i + on_circle_indx) % n], inner_pts[(i + 1 + on_circle_indx) % n]);
        inner_poly.add(seg);
    }

    /*
                    CSegment2D & seg = inner_poly.edges().front();
                    CPoint2 sp = seg.start();
                    CPoint2 ep = seg.end();
                    std::cout << sp.norm() << " " << ep.norm() << std::endl;
    */
    return 0;
}

inline CPoint2 CPolygon2D::mass_center()
{
    if (m_edges.empty())
        return CPoint2(0, 0);

    std::vector<CPoint2> pts;
    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment2D& s = m_edges[i];
        pts.push_back(s.start());
    }

    CPoint2 center(0, 0);
    double area = 0;

    for (int j = 1; j < pts.size() - 1; j++)
    {
        CPoint2 c = (pts[0] + pts[j] + pts[j + 1]) / 3;
        double s = (pts[j] - pts[0]) ^ (pts[j + 1] - pts[0]) / 2;
        area += s;
        center += c * s;
    }

    center /= area;

    return center;
}

bool inline CPolygon2D::intersectConvexPolygon(CRay2D& ray, CPoint2& p)
{
    CPoint2 d = ray.direction();
    d /= d.norm();

    CSegment2D seg(ray.start(), ray.start() + d * 100);

    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment2D& sg = m_edges[i];
        if (sg.intersect(seg, p))
            return true;
    }
    return false;
}

bool inline CPolygon2D::intersectConvexPolygon(CSegment2D& seg, CSegment2D& intersection)
{
    std::vector<CPoint2> inter_pts;
    CPoint2 s = seg.start();
    CPoint2 e = seg.end();

    bool s_inside = insideConvexPolygon(s);
    bool e_inside = insideConvexPolygon(e);

    if (s_inside && e_inside)
    {
        intersection = seg;
        return true;
    }

    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CPoint2 inter_pt;
        CSegment2D& sg = m_edges[i];
        if (sg.intersect(seg, inter_pt))
        {
            inter_pts.push_back(inter_pt);
        }
    }

    if (s_inside && !e_inside)
    {
        intersection.start() = s;
        intersection.end() = inter_pts.front();
        return true;
    }

    if (!s_inside && e_inside)
    {
        intersection.start() = e;
        intersection.end() = inter_pts.front();
        return true;
    }

    if (inter_pts.empty())
        return false;

    intersection.start() = inter_pts.front();
    intersection.end() = inter_pts.back();

    return true;
}

// check if a point is on the LEFT side of an edge
bool inside(CPoint2& p0, CPoint2& p1, CPoint2& p2)
{
    double a[2] = {p0[0], p0[1]};
    double b[2] = {p1[0], p1[1]};
    double c[2] = {p2[0], p2[1]};

    return orient2d(a, b, c) > 0;

    // return (p2[1] - p1[1]) * p[0] + (p1[0] - p2[0]) * p[1] + (p2[0] * p1[1] - p1[0] * p2[1]) < 0;
};

// calculate intersection point
CPoint2 intersection(CPoint2& cp1, CPoint2& cp2, CPoint2& s, CPoint2& e)
{
    CPoint2 dc = {cp1[0] - cp2[0], cp1[1] - cp2[1]};
    CPoint2 dp = {s[0] - e[0], s[1] - e[1]};

    double n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0];
    double n2 = s[0] * e[1] - s[1] * e[0];
    double n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0]);

    return {(n1 * dp[0] - n2 * dc[0]) * n3, (n1 * dp[1] - n2 * dc[1]) * n3};
}

// Sutherland-Hodgman clipping
void SutherlandHodgman(const std::vector<CPoint2>& subjectPolygon, const std::vector<CPoint2> clipPolygon, std::vector<CPoint2>& outputPolygon)
{

    CPoint2 cp1, cp2, s, e;
    std::vector<CPoint2> inputPolygon;

    // copy subject polygon to new polygon and set its size
    outputPolygon = subjectPolygon;

    // newPolygonSize = subjectPolygonSize;

    for (size_t j = 0; j < clipPolygon.size(); j++)
    {
        // copy new polygon to input polygon & set counter to 0
        inputPolygon = outputPolygon;

        outputPolygon.clear();

        // get clipping polygon edge
        cp1 = clipPolygon[j];
        cp2 = clipPolygon[(j + 1) % clipPolygon.size()];

        for (size_t i = 0; i < inputPolygon.size(); i++)
        {
            // get subject polygon edge
            s = inputPolygon[i];
            e = inputPolygon[(i + 1) % inputPolygon.size()];

            // Case 1: Both vertices are inside:
            // Only the second vertex is added to the output list
            if (inside(s, cp1, cp2) && inside(e, cp1, cp2))
                outputPolygon.push_back(e);

            // Case 2: First vertex is outside while second one is inside:
            // Both the point of intersection of the edge with the clip boundary
            // and the second vertex are added to the output list
            else if (!inside(s, cp1, cp2) && inside(e, cp1, cp2))
            {
                outputPolygon.push_back(intersection(cp1, cp2, s, e));
                outputPolygon.push_back(e);
            }

            // Case 3: First vertex is inside while second one is outside:
            // Only the point of intersection of the edge with the clip boundary
            // is added to the output list
            else if (inside(s, cp1, cp2) && !inside(e, cp1, cp2))
                outputPolygon.push_back(intersection(cp1, cp2, s, e));

            // Case 4: Both vertices are outside
            else if (!inside(s, cp1, cp2) && !inside(e, cp1, cp2))
            {
                // No vertices are added to the output list
            }
        }
    }
}

/*!
 *	subject polygon maybe not convex
 *   clip polygon is the background polygon, which must be convex
 *
 */
void SutherlandHodgman(CPolygon2D& subjectPolygon, CPolygon2D& clipPolygon, CPolygon2D& newPolygon)
{

    std::vector<CPoint2> subjectPoints;
    std::vector<CPoint2> clipPoints;
    std::vector<CPoint2> newPoints;

    for (size_t i = 0; i < subjectPolygon.edges().size(); i++)
    {
        subjectPoints.push_back(subjectPolygon.edges()[i].start());
    }
    for (size_t i = 0; i < clipPolygon.edges().size(); i++)
    {
        clipPoints.push_back(clipPolygon.edges()[i].start());
    }

    if (subjectPoints.size() > clipPoints.size())
        SutherlandHodgman(subjectPoints, clipPoints, newPoints);
    else
        SutherlandHodgman(clipPoints, subjectPoints, newPoints);

    for (size_t i = 0; i < newPoints.size(); i++)
    {
        CPoint2 sp = newPoints[i];
        CPoint2 ep = newPoints[(i + 1) % newPoints.size()];
        CSegment2D seg(sp, ep);
        newPolygon.add(seg);
    }
};

/*!
 *	Intersection between two convex polygon
 */
void convex_polygon_intersection(CPolygon2D& p1, CPolygon2D& p2, CPolygon2D& I)
{

    std::vector<std::pair<CPoint2, bool>> pts1;
    for (size_t i = 0; i < p1.edges().size(); i++)
    {
        CSegment2D& s = p1.edges()[i];
        CPoint2 p = s.start();
        bool inside = p2.insideConvexPolygon(p);
        pts1.push_back(std::pair<CPoint2, bool>(p, inside));
    }
    size_t head1, tail1;
    for (size_t i = 0; i < pts1.size(); i++)
    {
        bool is0 = pts1[(i + pts1.size() - 1) % pts1.size()].second;
        bool is1 = pts1[i].second;
        bool is2 = pts1[(i + 1) % pts1.size()].second;

        if (!is0 && is1)
            head1 = i;
        if (!is2 && is1)
            tail1 = i;
    }

    std::vector<std::pair<CPoint2, bool>> pts2;
    for (size_t i = 0; i < p2.edges().size(); i++)
    {
        CSegment2D& s = p2.edges()[i];
        CPoint2 p = s.start();
        bool inside = p1.insideConvexPolygon(p);
        pts2.push_back(std::pair<CPoint2, bool>(p, inside));
    }
    size_t head2, tail2;
    for (size_t i = 0; i < pts2.size(); i++)
    {
        bool is0 = pts2[(i + pts2.size() - 1) % pts2.size()].second;
        bool is1 = pts2[i].second;
        bool is2 = pts2[(i + 1) % pts2.size()].second;

        if (!is0 && is1)
            head2 = i;
        if (!is2 && is1)
            tail2 = i;
    }

    CSegment2D s2 = p2.edges()[(head2 + p2.edges().size() - 1) % p2.edges().size()];
    CSegment2D s1 = p1.edges()[tail1];
    CPoint2 q1;
    s1.intersect(s2, q1);

    CSegment2D t2 = p2.edges()[tail2];
    CSegment2D t1 = p1.edges()[(head1 + p1.edges().size() - 1) % p1.edges().size()];
    CPoint2 q2;
    t1.intersect(t2, q2);

    std::vector<CPoint2> pts;

    pts.push_back(q2);
    size_t i = head1;
    while (i != tail1)
    {
        pts.push_back(pts1[i].first);
        i = (i + 1) % pts1.size();
    }
    pts.push_back(q1);
    i = head2;
    while (i != tail2)
    {
        pts.push_back(pts2[i].first);
        i = (i + 1) % pts2.size();
    }

    for (size_t i = 0; i < pts.size(); i++)
    {
        I.edges().push_back(CSegment2D(pts[i], pts[(i + 1) % pts.size()]));
    }
}

}; // namespace MeshLib
#endif
