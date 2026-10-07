/*! \file PowerDynamicMesh.h
*   \brief Power Diagram, Power Delaunay Mesh
*   \author David Gu
*   \date   documented on 11/03/2020
*
*   Dynamic Mesh for power diagram and power Delaunay triangulation
*/
#ifndef  _POWER_DYNAMIC_MESH_H_
#define  _POWER_DYNAMIC_MESH_H_

#include <map>
#include <vector>
#include <stack>
#include <Eigen/Eigen>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/DynamicMesh.h"
#include "Mesh/Boundary.h"
#include "Mesh/Iterators.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"
#include "Geometry/Polygon2D.h"
#include "Geometry/Polygon3D.h"
#include "Geometry/Plane3D.h"
#include "Geometry/Transformation.h"
#include "Geometry/predicates.h"

namespace MeshLib
{
	/*! \brief CPowerFace class
	*
	*	Face class for power diagram
	*   Trait: face normal
	*/
	class CPowerFace : public CFace
	{
	protected:
		/*! face normal */
		CPoint m_normal;
		/*! face area  */
		double m_area;
		/*! face dual point */
		CPoint m_dual_point;
		/*! inside the target domain */
		bool m_inside;
		/*! normal map*/
		CPoint m_normal_map;
		/*! 2D transformation */
		CTransformation2D m_2d_transformation;
		/*! 3D transformation */
		CTransformation   m_transformation;

	public:
		/*! face normal */
		CPoint& normal() { return m_normal; };
		/*! face area */
		double& area() { return m_area; };
		/*! face dual point*/
		CPoint& dual_point() { return m_dual_point; };
		/*! inside the target domain*/
		bool& inside() { return m_inside; };

		/*! 2D transformation */
		CTransformation2D& planar_transformation() { return m_2d_transformation; };
		/*! 3D transformation */
		CTransformation& transformation() { return m_transformation; };

	};

	/*! \brief CPowerVertex class
	*
	*   Vertex class for power diagram
	*   Trait : vertex rgb color
	*/
	class CPowerVertex : public  CVertex
	{
	protected:
		/*! vertex rgb color */
		CPoint   m_rgb;
		/*! vertex convex hull point, Kantorovich potential phi*/
		CPoint   m_phi;

		/*! duall cell */
		CPolygon2D m_dual_cell2D;
		/*! dual cell in 3D*/
		CPolygon3D m_dual_cell3D;
		/*! normal map*/
		CPoint m_normal_map;

	public:
		/*! CViewerVertex Constructor */
		CPowerVertex()
		{
			m_rgb = CPoint(1, 1, 1); //default color is white
			m_normal = CPoint(1, 0, 0);
		}
		/*! vertex rgb color */
		CPoint& rgb() { return m_rgb; };
		/*! vertex Kantorovich potential phi*/
		CPoint& phi() { return m_phi; };

		/*! read vertex rgb, uv from vertex string */
		void _from_string();
		/*! write vertex uv to vertex string */
		void _to_string();


		/*! dual cell*/
		CPolygon2D& dual_cell2D() { return m_dual_cell2D; };

		/*! vertex conjudate cell */
		CPolygon3D& dual_cell3D() { return m_dual_cell3D; };

		/*! normal map*/
		CPoint& normal_map() { return m_normal_map; };

	};

	// read vertex rgb, uv from vertex string 
	inline void CPowerVertex::_from_string()
	{
		CParser parser(m_string);

		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken* token = *iter;
			if (token->m_key == "uv")
			{
				token->m_value >> m_uv;
			}
			else if (token->m_key == "rgb")
			{
				token->m_value >> m_rgb;
			}
			else if (token->m_key == "normal")
			{
				token->m_value >> m_normal_map;
			}
			else if (token->m_key == "u")
			{
				token->m_value >> m_phi;
			}

		}
	};

	// read vertex rgb, uv from vertex string 
	inline void CPowerVertex::_to_string()
	{
		CParser parser(m_string);

		parser._removeToken("uv");
		parser._toString(m_string);

		std::stringstream iss;
		iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ") ";

		if (m_string.length() > 0)
		{
			m_string += " ";
		}
		m_string += iss.str();

	};

	/*! \brief CPowerEdge class
	*
	*   Edge class for power diagram
	*   Trait : Edge sharp
	*/
	class CPowerEdge : public  CEdge
	{
	protected:
		/*! edge sharp */
		bool   m_sharp;
		/*! processed */
		bool m_processed;
	public:
		/*! COTVertex Constructor */
		CPowerEdge()
		{
			m_sharp = false;
		}
		/*! vertex rgb color */
		bool& sharp() { return m_sharp; };
		/*! processed used in Lawson Edge Flip*/
		bool& processed() { return m_processed; };

		/*! read vertex rgb, uv from vertex string */
		void _from_string();

	};


	// read vertex rgb, uv from vertex string 
	inline void CPowerEdge::_from_string()
	{
		CParser parser(m_string);

		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken* token = *iter;
			if (token->m_key == "sharp")
			{
				m_sharp = true;
			}
			break;
		}
	};

	/*! \brief CPowerHalfEdge class
	 *
	 *   HalfEdge class for power diagram
	 *   Trait : corner uv, corner id
	 */
	class CPowerHalfEdge : public  CHalfEdge
	{
	public:
		/*! corner uv */
		CPoint2& uv() { return m_uv; };
		/*! corner id*/
		int& cid() { return m_cid; };

	protected:
		/*! halfedge uv */
		CPoint2   m_uv;
		/*! halfedge corner id*/
		int m_cid;
	};

	/*-------------------------------------------------------------------------------------------------------------------------------------

		Power Diagram Mesh

	--------------------------------------------------------------------------------------------------------------------------------------*/
	/*! \brief CPowerDiagramMesh class
	*
	*	mesh class for optimal transport
	*
	*/
	template<typename V, typename E, typename F, typename H>
	class CPowerDiagramDynamicMesh : public CDynamicMesh<V, E, F, H>
	{
	public:

		typedef V CVertex;
		typedef E CEdge;
		typedef F CFace;
		typedef H CHalfEdge;

		typedef CBoundary<V, E, F, H> CBoundary;
		typedef CLoop<CVertex, CEdge, CFace, CHalfEdge> CLoop;
		typedef MeshVertexIterator<V, E, F, H> MeshVertexIterator;
		typedef MeshFaceIterator<V, E, F, H>   MeshFaceIterator;
		typedef MeshEdgeIterator<V, E, F, H>   MeshEdgeIterator;
		typedef VertexVertexIterator<V, E, F, H> VertexVertexIterator;
		typedef FaceVertexIterator<V, E, F, H> FaceVertexIterator;
		typedef VertexEdgeIterator<V, E, F, H> VertexEdgeIterator;
		typedef VertexFaceIterator<V, E, F, H> VertexFaceIterator;
		typedef VertexOutHalfedgeIterator<V, E, F, H> VertexOutHalfedgeIterator;
		typedef FaceHalfedgeIterator<V, E, F, H>  FaceHalfedgeIterator;

	protected:
		//compute the face dual point
		//each face is a plane z=ax+by-c, the dual point is (a,b,c)
			/*
			*	dual point of a face, Ax+By -C = Z, (A,B,C)
			*
			*   <p_k, (x,y,-1)> = z, k = 0, 1, 2
			*
			*	<p_1-p_0, (x,y,-1)> = 0,
			*   <p_2-p_0, (x,y,-1)> = 0,
			*
			*   (x,y,-1) parallel to the normal
			*/
		virtual void __face_dual_point(CFace* pf)
		{
			CPoint p[3];
			int i = 0;
			for (FaceVertexIterator fviter(pf); !fviter.end(); fviter++)
			{
				CVertex* pv = *fviter;
				p[i++] = pv->phi();
			}
			CPoint n = (p[1] - p[0]) ^ (p[2] - p[0]);

			double A = -n[0] / n[2];
			double B = -n[1] / n[2];
			double C = p[0] * n / n[2];

			pf->dual_point() = CPoint(A, B, -C);
			pf->inside() = (A * A + B * B < 1.0);

			n /= n.norm();
			pf->normal() = n;
		};

		//compute power cell clipped by polygon D
		virtual void __power_cell_clip(CVertex* pv, CPolygon2D& D)
		{
			std::vector<CPoint>  points;
			//whether the whole dual cell is inside the unit disk
			bool dual_cell_inside = true;

			for (VertexFaceIterator vfiter(pv); !vfiter.end(); vfiter++)
			{
				CFace* pf = *vfiter;
				points.push_back(pf->dual_point());
				if (!pf->inside()) dual_cell_inside = false;
			}

			CPolygon3D Poly;
			for (int i = 0; i < (int)points.size(); i++)
			{
				CSegment3D seg(points[i], points[(i + 1) % points.size()]);
				Poly.add(seg);
			}
			pv->dual_cell3D() = Poly;
			pv->dual_cell2D() = pv->dual_cell3D().project();

			//the whole dual cell is inside the target domain
			if (dual_cell_inside) return;

			//clip the dual cell with the unit circle
			//CPolygon2D inner_polygon;
			CPolygon2D inner_polygon;
			SutherlandHodgman(pv->dual_cell2D(), D, inner_polygon);
			pv->dual_cell2D() = inner_polygon;

			//compute the corresponding 3D conjugate cell
			CPolygon3D fpoly = pv->dual_cell3D();
			CPlane3D   plane = fpoly.plane();
			fpoly = transform(plane, fpoly);
			CPolygon3D dpoly;
			plane_cylinder_intersection(plane, D, dpoly);
			dpoly = transform(plane, dpoly);
			CPolygon2D inter_poly;

			auto fprj = fpoly.project();
			auto dprj = dpoly.project();
			SutherlandHodgman(fprj, dprj, inter_poly);
			pv->dual_cell3D() = inverse_transform(plane, inter_poly);
		};	
		//verify if an edge is local power Delaunay
		virtual bool __edge_local_power_Delaunay(CEdge* pe, bool nearest = true)
		{
			double p[4][3];
			CVertex* pV[4];

			CHalfEdge* ph = this->edgeHalfedge(pe, 0);
			pV[0] = this->halfedgeSource(ph);
			pV[1] = this->halfedgeTarget(ph);
			ph = this->halfedgeNext(ph);
			pV[2] = this->halfedgeTarget(ph);
			ph = this->edgeHalfedge(pe, 1);
			ph = this->halfedgeNext(ph);
			pV[3] = this->halfedgeTarget(ph);

			for (int k = 0; k < 4; k++)
			{
				CPoint& pt = pV[k]->phi();
				for (int i = 0; i < 3; i++)
				{
					p[k][i] = pt[i];
				}
			}

			// worst OT
			double sign = orient3d(p[0], p[1], p[2], p[3]);
			if ((nearest && sign > 0) || (!nearest && sign < 0))
				return false;
			return true;
		}
		//verify if an edge is flippable
		virtual bool __edge_flippable(CEdge* pe)
		{
			if (_boundary(pe))
			{
				return false;
			}

			double p[4][2];
			CVertex* pV[4];

			CHalfEdge* ph = this->edgeHalfedge(pe, 0);
			pV[0] = this->halfedgeSource(ph);
			pV[1] = this->halfedgeTarget(ph);
			ph = this->halfedgeNext(ph);
			pV[2] = this->halfedgeTarget(ph);
			ph = this->edgeHalfedge(pe, 1);
			ph = this->halfedgeNext(ph);
			pV[3] = this->halfedgeTarget(ph);

			for (int k = 0; k < 4; k++)
			{
				CPoint2& uv = pV[k]->uv();
				for (int i = 0; i < 2; i++)
				{
					p[k][i] = uv[i];
				}
			}

			if (orient2d(p[1], p[2], p[3]) < 0) //concave 
			{
				return false;
			}

			if (orient2d(p[0], p[3], p[2]) < 0) //concave 
			{
				return false;
			}

			return true;
		};


	public:
		/*!----------------------------------------------------------------------------------------------
		 * Handle infinite vertex
		-----------------------------------------------------------------------------------------------*/
		//add one infinity vertex, representing the dual support plane z=h
		void _add_infinity_vertex()
		{
			CBoundary bnd(this);
			CLoop* pL = bnd.loops().front();
			std::list<CHalfEdge*>& hes = pL->halfedges();

			//vertex at infinity
			CVertex* pV = this->createVertex(-1);
			pV->point() = CPoint(0, 0, 0);
			pV->uv() = CPoint2(0, 0);

			//the dual plane is a horizontal plane, z = 8.
			pV->phi() = CPoint(0, 0, -8);
			//generate the surrounding faces of the infinite vertex

			for (typename std::list<CHalfEdge*>::iterator hiter = hes.begin(); hiter != hes.end(); hiter++)
			{
				CHalfEdge* ph = *hiter;
				CVertex* pS = this->halfedgeSource(ph);
				CVertex* pT = this->halfedgeTarget(ph);
				std::vector<CVertex*> vs;
				vs.push_back(pV);
				vs.push_back(pT);
				vs.push_back(pS);
				int fid = this->numFaces();
				this->createFace(vs, -fid - 1);
			}

			//link one halfedge to the infinity vertex
			CHalfEdge* ph = hes.front();
			ph = this->halfedgeSym(ph);
			ph = this->halfedgeNext(ph);
			//vertexHalfedge(pV) = ph;
			pV->halfedge() = ph;

			for (MeshVertexIterator viter(this); !viter.end(); viter++)
			{
				CVertex* pV = *viter;
				pV->boundary() = false;
			}

			for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
			{
				CEdge* pE = *eiter;
				//pE->m_boundary = false;
			}

		}
		//whether a vertex is attached to the infinity vertex
		bool _infinite(CVertex* pv)
		{
			return pv->id() < 0;
		}
		//whether an edge is attached to the infinity vertex
		bool _infinite(CEdge* pe)
		{
			CVertex* pv = this->edgeVertex1(pe);
			if (_infinite(pv)) return true;
			pv = this->edgeVertex2(pe);
			return _infinite(pv);
		}
		//whether a halfedge is attached to the infinity vertex
		bool _infinite(CHalfEdge* ph)
		{
			CVertex* pv = this->halfedgeTarget(ph);
			if (_infinite(pv)) return true;
			pv = this->halfedgeSource(ph);
			return _infinite(pv);
		}
		//whether a face is attached to the infinity vertex
		bool _infinite(CFace* pf)
		{
			CHalfEdge* ph = this->faceHalfedge(pf);
			if (_infinite(ph)) return true;
			ph = this->halfedgeNext(ph);
			if (_infinite(ph)) return true;
			ph = this->halfedgeNext(ph);
			return _infinite(ph);
		}
		//verify if an edge is on the link of the inifinity vertex
		bool _boundary(CEdge* pe)
		{
			CHalfEdge* ph = this->edgeHalfedge(pe, 0);
			ph = this->halfedgeNext(ph);
			CVertex* pv = this->halfedgeTarget(ph);
			if (_infinite(pv)) return true;
			ph = this->edgeHalfedge(pe, 1);
			ph = this->halfedgeNext(ph);
			pv = this->halfedgeTarget(ph);
			return _infinite(pv);
		}
		//verify if a vertex is on the link of the infinity vertex
		bool _boundary(CVertex* pV)
		{
			for (VertexOutHalfedgeIterator vhiter(this, pV); !vhiter.end(); vhiter++)
			{
				CHalfEdge* ph = *vhiter;
				if (_infinite(ph)) return true;
			}
			return false;
		}


		/*!----------------------------------------------------------------------------------------------
		 * methods to compute the convex hull, Power Delaunay triangulation
		-----------------------------------------------------------------------------------------------*/

		//edge swap to compute convex hull
		bool _edge_swap(bool nearest = true)
		{
			int flipped_edges = 0;
			int non_flipplable = 0;

			while (true)
			{
				flipped_edges = 0;
				non_flipplable = 0;

				for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
				{
					CEdge* pe = *eiter;
					//if (this->isBoundary(pe))
					//skip all the edges connecting to the infinity vertex
					//skill the link of the infinity vertex
					if (_boundary(pe) || _infinite(pe)) continue;

					if (!__edge_local_power_Delaunay(pe))
					{
						if (__edge_flippable(pe))
						{
							swapEdge(pe);
							flipped_edges++;
						}
						else
							non_flipplable++;
					}
				}

				if (flipped_edges == 0) break;
			}

			if (non_flipplable > 0)
			{
				std::cout << "Nonconvex" << std::endl;
				return false;
			}
			std::cout << "Convex" << std::endl;
			return true;
		}
		//Charles Lawson algorithm, edge swap for nearest and farthest Delaunay Triangulation
		bool _Lawson_edge_swap(bool nearest = true)
		{
			int flipped_edges = 0;
			int non_flipplable = 0;

			std::queue<CEdge*> illegal_edges;

			//unmark all edges
			for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
			{
				CEdge* pe = *eiter;
				pe->processed() = false;
			}

			//enqueue illegal edges
			for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
			{
				CEdge* pe = *eiter;
				if (_boundary(pe) || _infinite(pe))	continue;
				if (!__edge_local_power_Delaunay(pe, nearest))
				{
					pe->processed() = true;
					illegal_edges.push(pe);
				}
			}
			//process all illegal edges
			while (!illegal_edges.empty())
			{
				CEdge* pe = illegal_edges.front();
				illegal_edges.pop();
				pe->processed() = false;

				if (!__edge_local_power_Delaunay(pe, nearest))
				{
					if (__edge_flippable(pe))
					{
						this->swapEdge(pe);

						std::vector<CEdge*> nbh_edges;

						CHalfEdge* ph = this->edgeHalfedge(pe, 0);
						nbh_edges.push_back(this->halfedgeEdge(this->halfedgeNext(ph)));
						nbh_edges.push_back(this->halfedgeEdge(this->halfedgePrev(ph)));
						CHalfEdge* ps = this->edgeHalfedge(pe, 1);
						nbh_edges.push_back(this->halfedgeEdge(this->halfedgeNext(ps)));
						nbh_edges.push_back(this->halfedgeEdge(this->halfedgePrev(ps)));

						for (size_t k = 0; k < nbh_edges.size(); k++)
						{
							CEdge* pE = nbh_edges[k];
							if (pE->processed()) continue;
							pE->processed() = true;
							illegal_edges.push(pE);
						}
					}
				}
			}

			//enqueue illegal edges
			for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
			{
				CEdge* pe = *eiter;
				if (_boundary(pe) || _infinite(pe))	continue;
				if (!__edge_local_power_Delaunay(pe, nearest))
				{
					std::cout << "Non Convex" << std::endl;
					return false;
				}
			}

			return true;
		}

		/*!----------------------------------------------------------------------------------------------
		 * methods to compute upper envelop, Power Voronoi Diagram
		-----------------------------------------------------------------------------------------------*/

		//Legendre transform
		virtual void _Legendre_transform(bool nearest)
		{

			for (MeshFaceIterator fiter(this); !fiter.end(); fiter++)
			{
				CFace* pf = *fiter;
				__face_dual_point(pf);
			}
		};
		//Legendre transform
		virtual void _Legendre_transform(bool nearest, CPolygon2D& D)
		{

			for (MeshFaceIterator fiter(this); !fiter.end(); fiter++)
			{
				CFace* pf = *fiter;
				__face_dual_point(pf);
				CPoint pt = pf->dual_point();
				CPoint2 pt2(pt[0], pt[1]);
				bool inside = D.insideConvexPolygon(pt2);
				pf->inside() = inside;
			}
		};

		/*!----------------------------------------------------------------------------------------------
		 * methods to clip power cells
		-----------------------------------------------------------------------------------------------*/

		//clip the power cells by D
		virtual void _power_cell_clip(CPolygon2D& D)
		{

			double total_area = 0;
			for (MeshVertexIterator viter(this); !viter.end(); viter++)
			{
				CVertex* pv = *viter;
				if (_infinite(pv)) continue;
				__power_cell_clip(pv, D);
				pv->dual_area() = pv->dual_cell2D().area();
				total_area += pv->dual_area();
			}
			std::cout << "Total dual area is " << total_area << std::endl;
		}
		//compute the mass center of the power cells
		virtual void _power_cell_center()
		{
			for (MeshVertexIterator viter(this); !viter.end(); viter++)
			{
				CVertex* pv = *viter;
				if (_infinite(pv)) continue;
				if (pv->dual_cell2D().edges().empty()) {
					std::cout << "Warning: dual cell empty" << pv->id() << std::endl;
					continue;
				}
				pv->uv() = pv->dual_cell2D().mass_center();
			}
		}

		//compute the Jacobian map from vertex->uv() to vertex->point()
		void _compute_Jacobian()
		{
			for (MeshFaceIterator fiter(this); !fiter.end(); fiter++)
			{
				CFace* pf = *fiter;
				std::vector<CPoint>  pts;
				std::vector<CPoint2> rts;
				std::vector<CPoint2> qts;
				for (FaceVertexIterator fviter(pf); !fviter.end(); fviter++)
				{
					CVertex* pV = *fviter;
					pts.push_back(pV->point());
					qts.push_back(pV->uv());
				}
				CTransformation T(pts[0], pts[1], pts[2]);
				CTransformation I = T.inverse();

				for (int i = 0; i < 3; i++)
				{
					CPoint q = T * pts[i];
					rts.push_back(CPoint2(q[0], q[1]));
				}

				CTransformation2D tr(qts, rts);
				pf->planar_transformation() = tr;
				pf->transformation() = I;

				/*
					for (int i = 0; i < 3; i++)
					{
						CPoint  R = pf->transformation() * (pf->planar_transformation() * qts[i]);
						std::cout << (pts[i] - R).norm() << std::endl;
					}
				*/
			}
		};
	};

	typedef CPowerDiagramDynamicMesh<CPowerVertex, CPowerEdge, CPowerFace, CPowerHalfEdge> CPDMesh;


}
#endif  //! _POWER_DYNAMIC_MESH_H_

