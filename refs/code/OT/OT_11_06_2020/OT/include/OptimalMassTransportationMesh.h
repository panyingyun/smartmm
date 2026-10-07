/*! \file  OptimalMassTransportationMesh.h
*   \brief Optimal Mass Transport DynamicMesh
*   \author David Gu
*   \date   documented on 11/03/2020
*
*   Dynamic Mesh for viewer
*/
#ifndef  _OPTIMAL_MASS_TRANSPORT_DYNAMIC_MESH_H_
#define  _OPTIMAL_MASS_TRANSPORT_DYNAMIC_MESH_H_

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
#include "PowerDiagramMesh.h"

#define T_PD_MESH CPowerDiagramDynamicMesh<V, E, F, H>

namespace MeshLib
{
	/*! \brief COTFace class
	*
	*	Face class for optimal mass transport
	*   Trait: face normal
	*/
	class COMTFace : public CPowerFace
	{
	protected:
	public:
	};

	/*! \brief COTVertex class
	*
	*   Vertex class for optimal mass transport
	*   Trait : vertex rgb color
	*/
	class COMTVertex : public  CPowerVertex
	{
	protected:
		/*! vertex target area */
		double   m_target_area;
		/*! power cell area */
		double   m_dual_area;
		/*! vertex index */
		int		 m_index;
		/*! update direction */
		double	 m_update_direction;
		/*! backup height */
		double	 m_backup_height;

	public:
		/*! CViewerVertex Constructor */
		COMTVertex() { m_index = 0; m_dual_area = 0; m_target_area = 0; }
		/*! vertex target area */
		double& target_area() { return m_target_area; };
		/*! vertex area */
		double& area() { return m_target_area; };
		/*! vertex dual area */
		double& dual_area() { return m_dual_area; };
		/*! vertex index*/
		int& index() { return m_index; };
		/*! update direction */
		double& update_direction() { return m_update_direction; };
		/*! backup height */
		double& backup_height() { return m_backup_height; };
	};


	/*! \brief COTEdge class
	*
	*   Edge class for optimal transport
	*   Trait : Edge sharp
	*/
	class COMTEdge : public  CPowerEdge
	{
	protected:
		/*! edge weight */
		double m_weight;
	public:
		/*! edge weight for Hessian matrix */
		double& weight() { return m_weight; };
	};

	/*! \brief COMTHalfEdge class
	 *
	 *   HalfEdge class for optimal mass transport
	 */
	class COMTHalfEdge : public  CPowerHalfEdge
	{
	public:
	protected:
	};

	/*-------------------------------------------------------------------------------------------------------------------------------------

		Optimal Mass Transport Mesh

	--------------------------------------------------------------------------------------------------------------------------------------*/
	/*! \brief COMTMesh class
	*
	*	mesh class for optimal mass transport
	*
	*/
	template<typename V, typename E, typename F, typename H>
	class COMTDynamicMesh : public CPowerDiagramDynamicMesh<V, E, F, H>
	{
	public:
		typedef typename T_PD_MESH::MeshVertexIterator MeshVertexIterator;
        typedef typename T_PD_MESH::MeshEdgeIterator MeshEdgeIterator;

		//newton's method for computing update direction
		virtual void _update_direction(CPolygon2D& D);
		//compute L2 error and relative error
		virtual double _error();

	protected:
		//compute one edge weight, intersecting with the polygon D
		virtual void __edge_weight(E* pe, CPolygon2D& D);
		//compute the edge weight
		virtual void _edge_weight(CPolygon2D& D);
		//compute the gradient
		virtual void _calculate_gradient();
		//compute the Hessian
		virtual void __calculate_hessian(Eigen::SparseMatrix<double>& m_hessian);
		//solve the Hessian system
		virtual int  __solve(Eigen::SparseMatrix<double>& m_hessian);

	};

	typedef MeshLib::COMTDynamicMesh<COMTVertex, COMTEdge, COMTFace, COMTHalfEdge> COMTMesh;

	/*! Compute L2 error and maximal relative error */
	template<typename V, typename E, typename F, typename H>
	double COMTDynamicMesh<V, E, F, H>::_error()
	{
		double max_error = -1e+10;
		double total_error = 0;

		for (MeshVertexIterator viter(this); !viter.end(); viter++)
		{
			V* pv = *viter;
			if (this->_infinite(pv)) continue;

			double   da = fabs(pv->target_area() - pv->dual_area());
			double error = da / pv->target_area();
			if (error > max_error)
			{
				max_error = error;
			}
			total_error += da * da;
		}
		std::cout << "Max relative error is " << max_error << " Total L2 error is " << total_error << std::endl;
		return max_error;
	}

	/*! compute the gradient of the energy */
	template<typename V, typename E, typename F, typename H>
	void COMTDynamicMesh<V, E, F, H>::_calculate_gradient()
	{
		for (MeshVertexIterator viter(this); !viter.end(); viter++)
		{
			V* pv = *viter;
			if (this->_infinite(pv)) continue;
			pv->dual_area() = pv->dual_cell2D().area();
			pv->update_direction() = pv->target_area() - pv->dual_area();
		}
	};

	//compute the edge weight for an interior edge
	template<typename V, typename E, typename F, typename H>
	void COMTDynamicMesh<V, E, F, H>::__edge_weight(E* pe, CPolygon2D& D)
	{
		F* pf1 = this->edgeFace1(pe);
		F* pf2 = this->edgeFace2(pe);

		V* pv1 = this->edgeVertex1(pe);
		V* pv2 = this->edgeVertex2(pe);

		CPoint2 duv1(pf1->dual_point()[0], pf1->dual_point()[1]);
		CPoint2 duv2(pf2->dual_point()[0], pf2->dual_point()[1]);

		CPoint2 uv1(pv1->phi()[0], pv1->phi()[1]);
		CPoint2 uv2(pv2->phi()[0], pv2->phi()[1]);

		CSegment2D seg(duv1, duv2);
		CSegment2D intersection;
		if (D.intersectConvexPolygon(seg, intersection))
		{
			pe->weight() = (intersection.start() - intersection.end()).norm() / (uv1 - uv2).norm();
			return;
		}
		pe->weight() = 0;

	};

	//compute the weights for all edges, boundary dual edges need to be clipped by D */
	template<typename V, typename E, typename F, typename H>
	void COMTDynamicMesh<V, E, F, H>::_edge_weight(CPolygon2D& D)
	{
		for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
		{
			E* pe = *eiter;
			__edge_weight(pe, D);
		}
	};

	/*! Compute Hessian Matrix */
	template<typename V, typename E, typename F, typename H>
	void COMTDynamicMesh<V, E, F, H>::__calculate_hessian(Eigen::SparseMatrix<double>& m_hessian)
	{
		m_hessian.resize(this->numVertices() - 1, this->numVertices() - 1);
		std::vector<Eigen::Triplet<double>> triplets;


		for (MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
		{
			E* pe = *eiter;
			V* pv0 = this->edgeVertex1(pe);
			V* pv1 = this->edgeVertex2(pe);

			if (this->_infinite(pv0) || this->_infinite(pv1)) continue;

			int idx0 = pv0->index();
			int idx1 = pv1->index();
			double h = pe->weight();

			triplets.push_back({ idx0, idx1,  h / 2 });
			triplets.push_back({ idx1, idx0,  h / 2 });
			triplets.push_back({ idx0, idx0, -h / 2 });
			triplets.push_back({ idx1, idx1, -h / 2 });
		}

		double epsilon = -1e-8;

		for (int i = 0; i < this->numVertices() - 1; i++)
		{
			triplets.push_back({ i,i, epsilon });
		}

		m_hessian.setFromTriplets(triplets.begin(), triplets.end());

		//std::cout << m_hessian;
	}

	/*! solve linear system */
	template<typename V, typename E, typename F, typename H>
	int COMTDynamicMesh<V, E, F, H>::__solve(Eigen::SparseMatrix<double>& m_hessian)
	{

		Eigen::VectorXd m_gradient;
		m_gradient.resize(this->numVertices() - 1);

		for (MeshVertexIterator viter(this); !viter.end(); viter++)
		{
			V* pv = *viter;
			if (this->_infinite(pv)) continue;
			m_gradient[pv->index()] = pv->update_direction();
		}

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		Eigen::VectorXd m_x;

		//hessian.coeffRef(0, 0) -= 1;
		solver.compute(m_hessian);
		m_x = solver.solve(m_gradient);
		// use conjugate gradient solver here
		//solver2.compute(hessian);
		//x = solver2.solve(gradient);
		m_x.array() -= m_x.mean();

		for (MeshVertexIterator viter(this); !viter.end(); viter++)
		{
			V* pv = *viter;
			if (this->_infinite(pv)) continue;
			pv->update_direction() = m_x[pv->index()];
		}

		//std::cout << m_hessian << std::endl;

		//std::cout << m_x << std::endl;
		return 0;
	}

	//use Newton's method to compute the update direction
	template<typename V, typename E, typename F, typename H>
	void COMTDynamicMesh<V, E, F, H>::_update_direction(CPolygon2D& D)
	{
		_calculate_gradient();
		_edge_weight(D);
		Eigen::SparseMatrix<double> hessian;
		__calculate_hessian(hessian);
		__solve(hessian);
	}

	//damping algorithm for computing OT Map
	template<class COMTMesh>
	void _OT_damping(COMTMesh& mesh, CPolygon2D& D, double step, const bool nearest)
	{
		//set the step length
		double step_length = (nearest) ? step : -step;

		//backup the current height
		for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
		{
			typename COMTMesh::CVertex* pv = *viter;
			if (mesh._infinite(pv)) continue;
			pv->backup_height() = pv->phi()[2];
		}

		//make sure all the samples with heights are on the convex position
		while (true)
		{
			for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
			{
				typename COMTMesh::CVertex* pv = *viter;
				if (mesh._infinite(pv)) continue;
				pv->phi()[2] += step_length * pv->update_direction();
			}
			//construct convex hull
			if (mesh._Lawson_edge_swap(nearest)) break;
			//roll back
			std::cout << "Roll Back" << std::endl;

			for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
			{
				typename COMTMesh::CVertex* pv = *viter;
				if (mesh._infinite(pv)) continue;
				pv->phi()[2] = pv->backup_height();
			}
			step_length /= 2.0;
		}
		//make sure all power cells are non-empty
		while (true)
		{
			mesh._Legendre_transform(nearest, D);
			mesh._power_cell_clip(D);

			bool legal = true;
			for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
			{
				typename COMTMesh::CVertex* pv = *viter;
				if (mesh._infinite(pv)) continue;
				if (fabs(pv->dual_area()) < pv->target_area() * 1e-10) {
					legal = false;
					break;
				}
			}
			if (legal) break;
			//roll back
			std::cout << "Empty power cell. Roll Back" << std::endl;
			for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
			{
				typename COMTMesh::CVertex* pv = *viter;
				if (mesh._infinite(pv)) continue;
				pv->phi()[2] = pv->backup_height();
			}
			step_length /= 2.0;
			for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
			{
				typename COMTMesh::CVertex* pv = *viter;
				if (mesh._infinite(pv)) continue;
				pv->phi()[2] += step_length * pv->update_direction();
			}
			mesh._Lawson_edge_swap(nearest);
		}
	};

	//Newton's method for computing OT map
	template<class COMTMesh>
	void OT_Newton(COMTMesh& mesh, CPolygon2D& D, const double step_length, const bool nearest)
	{
		mesh._update_direction(D);
		_OT_damping<COMTMesh>(mesh, D, step_length, nearest);
		mesh._error();

		/*! debug */
		double epsilon = 1e-10;
		for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
		{
			typename COMTMesh::CVertex* pv = *viter;
			if (mesh._infinite(pv)) continue;
			bool legal = false;
			for (typename COMTMesh::VertexEdgeIterator veiter(pv); !veiter.end(); veiter++)
			{
				typename COMTMesh::CEdge* pe = *veiter;
				if (fabs(pe->weight()) > epsilon)
				{
					legal = true;
					break;
				}
			}
			if (!legal)
			{
				std::cout << "Illegal Vertex " << pv->id() << std::endl;
				break;
			}
		}
	};

	//initalize for computing OT map
	template<class COMTMesh>
	void OT_Initialize(COMTMesh& mesh, CPolygon2D& D, const bool nearest)
	{
		int idx = 0;
		for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
		{
			typename COMTMesh::CVertex* pv = *viter;
			pv->index() = idx++;
		}

		COperator<COMTMesh> pS(&mesh);
		pS._calculate_face_vertex_area();
	
		double s = 0;
		for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
		{
			typename COMTMesh::CVertex* pv = *viter;
			s += pv->target_area();
		}
		double ds = D.area();

		for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
		{
			typename COMTMesh::CVertex* pv = *viter;
			pv->target_area() *= (ds / s);
		}

		mesh._add_infinity_vertex();

		//set the inital height funciton
		for (typename COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
		{
			typename COMTMesh::CVertex* pv = *viter;
			//specially handle the inifinity vertex
			if (mesh._infinite(pv))
			{
				pv->uv() = CPoint2(0, 0);
				pv->phi() = CPoint(pv->uv()[0], pv->uv()[1], (nearest)? -2.0:2.0);
				continue;
			}
			double height = pv->uv().norm2() / 2.0;
			pv->phi() = CPoint(pv->uv()[0], pv->uv()[1], (nearest)?height: -height);
		}
		mesh._Lawson_edge_swap(nearest);
		mesh._Legendre_transform(nearest, D);
		mesh._power_cell_clip(D);
	}
}

#undef T_PD_MESH
#endif  //! OPTIMAL_TRANSPORT_DYNAMIC_MESH_H_
