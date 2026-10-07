#ifndef _VIEWER_COMTMesh_H_
#define _VIEWER_COMTMesh_H_

/*!
*      \file COMTMeshViewer.h
*      \brief Optimal Mass Transportation Map Mesh Viewer
*	   \author David Gu
*      \date 05/01/2020
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>

#ifdef MAC_OS
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif // MAC_OS

#include "viewer/Arcball.h"     /* Arc Ball  Interface */
#include "bmp/RgbImage.h"
#include "Operator/Operator.h"
#include "OptimalMassTransportationMesh.h"
#include "Geometry/Segment2D.h"

namespace MeshLib
{
	class COMTMeshViewerConfiguration
	{
	public:
		//compute optimal transportation map or worest transportation map
		bool m_nearest;
		//show prime or dual, or both
		int	 m_prime_dual_view;
		//show 3D surface or planar projection
		bool m_show_uv;
		//show edges
		bool m_show_edges;
		//show backface
		bool m_show_backface;

		bool m_show_boundary;
		bool m_show_sharp_edges;
		bool m_light_edit;
		bool m_show_normal_map;
		bool m_show_cell_center;
		bool m_show_dual_cell_boundary;

		//shade mode: 1 vertex normal; 0 face normal
		int  m_shade_mode;
		//texture name
		GLuint m_texture_name;
		//texture mode
		int    m_texture_mode;

	public:
		COMTMeshViewerConfiguration()
		{
			m_nearest = true;
			m_show_uv = false;
			m_show_edges = false;
			m_show_sharp_edges = false;
			m_show_backface = true;
			m_show_boundary = false;
			m_light_edit = false;
			m_prime_dual_view = 0;
			m_show_normal_map = false;
			m_show_cell_center = false;
			m_show_dual_cell_boundary = false;
			m_shade_mode = 1;
			m_texture_mode = 2;
		}
	};

	void drawString(double x, double y, double z, const char* string, double scale) {

		// Save the current matrix
		glPushMatrix();
		// Translate to the appropriate starting point
		glTranslated(x, y, z);

		// Render the characters
		glScaled(scale, scale, scale);
		for (const char* c = string; *c != '\0'; c++) {
			glutStrokeCharacter(GLUT_STROKE_ROMAN, *c);
		}
		// Another useful function
		//    int glutStrokeWidth(void *font, int character);
		// Retrieve the original matrix
		glPopMatrix();
	};

	void draw_prime_edges(COMTMesh* pMesh, COMTMeshViewerConfiguration& configure)
	{
		glDisable(GL_LIGHTING);
		glColor3f(0, 0, 0);
		glBegin(GL_LINES);

		for (COMTMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); ++eiter)
		{
			COMTMesh::CEdge* pE = *eiter;
			if (pMesh->_infinite(pE)) continue;
			COMTMesh::CVertex* v1 = pMesh->edgeVertex1(pE);
			COMTMesh::CVertex* v2 = pMesh->edgeVertex2(pE);

			CPoint p1 = v1->phi();
			CPoint p2 = v2->phi();
			if (!configure.m_show_uv)
			{
				glVertex3d(p1[0], p1[1], p1[2]);
				glVertex3d(p2[0], p2[1], p2[2]);
			}
			else
			{
				glVertex3d(p1[0], p1[1], 0);
				glVertex3d(p2[0], p2[1], 0);
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}


	void draw_cell_center(COMTMesh *pMesh)
	{

		glDisable(GL_LIGHTING);
		glColor3f(0, 0, 1);
		for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
		{
			COMTMesh::CVertex* pv = *viter;
			glPushMatrix();
			glTranslated(pv->uv()[0], pv->uv()[1], 0);
			glutSolidSphere(0.002, 8, 8);
			glPopMatrix();
		}
		glEnable(GL_LIGHTING);

	}
	void draw_dual_cell_boundary(COMTMesh* pMesh, COMTMeshViewerConfiguration & configure)
	{
		glDisable(GL_LIGHTING);
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);

		for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
		{
			COMTMesh::CVertex* pv = *viter;
			if (pv->dual_cell3D().edges().empty()) continue;
			for (size_t i = 0; i < pv->dual_cell3D().edges().size(); i++)
			{
				CSegment3D& s = pv->dual_cell3D().edges()[i];
				if (configure.m_show_uv)
				{
					glVertex3d(s.start()[0], s.start()[1], 0);
					glVertex3d(s.end()[0], s.end()[1], 0);
				}
				else
				{
					glVertex3d(s.start()[0], s.start()[1], s.start()[2]);
					glVertex3d(s.end()[0], s.end()[1], s.end()[2]);
				}
			}
		}

		if (configure.m_show_uv)
		{
			glColor3f(1, 0, 0);
			for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
			{
				COMTMesh::CVertex* pv = *viter;
				if (pv->dual_cell2D().edges().empty()) continue;
				for (size_t i = 0; i < pv->dual_cell2D().edges().size(); i++)
				{
					CSegment2D& s = pv->dual_cell2D().edges()[i];
					glVertex3d(s.start()[0], s.start()[1], 0);
					glVertex3d(s.end()[0], s.end()[1], 0);
				}
			}
			glColor3f(0, 0, 1);
		}
		glEnd();

/*
		for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
		{
			COMTMesh::CVertex* pv = *viter;
			char line[1024];
			sprintf_s(line, "%d", pv->id());
			drawString(pv->uv()[0], pv->uv()[1], 0, line, 1e-4);
		}
*/
		glEnable(GL_LIGHTING);
	}

	void draw_sharp_edges(COMTMesh* pMesh)
	{
		glDisable(GL_LIGHTING);
		glColor3f(1, 0, 0);
		glLineWidth(2.0);
		glBegin(GL_LINES);
		for (COMTMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); ++eiter)
		{
			COMTMesh::CEdge* pE = *eiter;
			if (!pE->sharp()) continue;

			COMTMesh::CVertex* pV0 = pMesh->edgeVertex1(pE);
			COMTMesh::CVertex* pV1 = pMesh->edgeVertex2(pE);

			CPoint p0 = pV0->phi();
			CPoint p1 = pV1->phi();

			glVertex3d(p0[0], p0[1], p0[2]);
			glVertex3d(p1[0], p1[1], p1[2]);
		}
		glEnd();
		glLineWidth(1.0);
		glEnable(GL_LIGHTING);
	}

	void draw_boundary(COMTMesh* pMesh, bool show_uv)
	{
		glDisable(GL_LIGHTING);
		glColor3f(0, 0, 1);
		glLineWidth(3.0);
		glBegin(GL_LINES);
		for (COMTMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); ++eiter)
		{
			COMTMesh::CEdge* pE = *eiter;
			if (!pE->boundary() && !pMesh->_boundary(pE)) continue;

			COMTMesh::CVertex* pV0 = pMesh->edgeVertex1(pE);
			COMTMesh::CVertex* pV1 = pMesh->edgeVertex2(pE);

			CPoint p0 = pV0->phi();
			CPoint p1 = pV1->phi();

			if (!show_uv)
			{
				glVertex3d(p0[0], p0[1], p0[2]);
				glVertex3d(p1[0], p1[1], p1[2]);
			}
			else
			{
				glVertex3d(p0[0], p0[1], 0);
				glVertex3d(p1[0], p1[1], 0);
			}
		}
		glEnd();
		glLineWidth(1.0);
		glEnable(GL_LIGHTING);
	}

	void draw_dual_edges(COMTMesh* pMesh, COMTMeshViewerConfiguration & configure )
	{
		glDisable(GL_LIGHTING);
		glLineWidth(1.0);
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);
		for (COMTMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); ++eiter)
		{
			COMTMesh::CEdge* pE = *eiter;
			COMTMesh::CFace* f1 = pMesh->edgeFace1(pE);
			COMTMesh::CFace* f2 = pMesh->edgeFace2(pE);

			CPoint p1 = f1->dual_point();
			CPoint p2 = f2->dual_point();

			if (!configure.m_show_uv)
			{
				glVertex3d(p1[0], p1[1], p1[2]);
				glVertex3d(p2[0], p2[1], p2[2]);
			}
			else
			{
				glVertex3d(p1[0], p1[1], 0);
				glVertex3d(p2[0], p2[1], 0);
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}

	/*! draw mesh */
	void _draw_dual_cell_face_2D(COMTMesh * pMesh, COMTMeshViewerConfiguration & configure)
	{

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glBegin(GL_TRIANGLES);
		for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
		{
			COMTMesh::CVertex* pv = *viter;
			if (pMesh->_infinite(pv)) continue;
			if (pv->dual_cell2D().edges().empty()) continue;

			CPoint n = (configure.m_show_normal_map)?pv->normal_map():pv->normal();

			glNormal3d(n[0], n[1], n[2]);
			glColor3d(pv->rgb()[0], pv->rgb()[1], pv->rgb()[2]);

			std::vector<CPoint2> ps;
			
			for (size_t i = 0; i < pv->dual_cell2D().edges().size(); i++)
			{
				ps.push_back(pv->dual_cell2D().edges()[i].start());
			}
			for (size_t i = 0; i < ps.size() - 2; i++)
			{
				CPoint2 ts[3];
				ts[0] = ps[0]; ts[1] = ps[i + 1]; ts[2] = ps[i + 2];
				for (int j = 0; j < 3; j++)
				{
					glVertex3d(ts[j][0], ts[j][1], 0);
				}
			}
		}
		glEnd();

		if (configure.m_show_backface)
		{
			glBegin(GL_TRIANGLES);
			for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
			{
				COMTMesh::CVertex* pv = *viter;
				if (pMesh->_infinite(pv)) continue;
				if (pv->dual_cell2D().edges().empty()) continue;

				CPoint n = (configure.m_show_normal_map) ? pv->normal_map() : pv->normal();

				glNormal3d(-n[0], -n[1], -n[2]);
				glColor3d(pv->rgb()[0], pv->rgb()[1], pv->rgb()[2]);

				std::vector<CPoint2> ps;

				for (size_t i = 0; i < pv->dual_cell2D().edges().size(); i++)
				{
					ps.push_back(pv->dual_cell2D().edges()[i].start());
				}
				for (size_t i = 0; i < ps.size() - 2; i++)
				{
					CPoint2 ts[3];
					ts[0] = ps[0]; ts[1] = ps[i + 1]; ts[2] = ps[i + 2];
					for (int j = 2; j >= 0; j--)
					{
						glVertex3d(ts[j][0], ts[j][1], 0);
					}
				}
			}
			glEnd();

		}
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	/*! draw dual mesh */
	/*
	void _draw_dual_cell_face_3D(COMTMesh* pMesh, COMTMeshViewerConfiguration & configure)
	{

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glBegin(GL_TRIANGLES);
		for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
		{
			COMTMesh::CVertex* pv = *viter;
			std::vector<CPoint> ps;

			if (pMesh->_infinite(pv)) continue;

			bool cell_inside = true;
			for (COMTMesh::VertexFaceIterator vfiter(pv); !vfiter.end(); ++vfiter)
			{
				COMTMesh::CFace* pf = *vfiter;
				CPoint pt = pf->dual_point();
				ps.push_back(pt);
				if (!pf->inside()) cell_inside = false;
			}

			CPoint n = (configure.m_show_normal_map)?pv->normal_map():pv->normal();
			glNormal3d(n[0], n[1], n[2]);

			if (cell_inside)
				glColor3d(pv->rgb()[0], pv->rgb()[1], pv->rgb()[2]);
			else
				glColor3d(1, 0, 0);

			for (size_t i = 0; i < ps.size() - 2; i++)
			{
				CPoint ts[3];
				ts[0] = ps[0]; ts[1] = ps[i + 1]; ts[2] = ps[i + 2];
				for (int j = 0; j < 3; j++)
				{
					glVertex3d(ts[j][0], ts[j][1], ts[j][2]);
				}
			}
		}

		if (configure.m_show_backface)
		{
			glBegin(GL_TRIANGLES);

			for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
			{
				COMTMesh::CVertex* pv = *viter;
				if (pMesh->_infinite(pv)) continue;
				
				std::vector<CPoint> ps;

				for( COMTMesh::VertexFaceIterator vfiter(pv); !vfiter.end(); ++vfiter )
				{
					COMTMesh::CFace* pf = *vfiter;
					CPoint pt = pf->dual_point();
					ps.push_back(pt);
				}

				CPoint n = (configure.m_show_normal_map)?pv->normal_map():pv->normal();
				glNormal3d(-n[0], -n[1], -n[2]);
				glColor3d(pv->rgb()[0], pv->rgb()[1], pv->rgb()[2]);

				for (size_t i = 0; i < ps.size() - 2; i++)
				{
					CPoint ts[3];
					ts[0] = ps[0]; ts[1] = ps[i + 1]; ts[2] = ps[i + 2];
					for (int j = 2; j >= 0; j--)
					{
						glVertex3d(ts[j][0], ts[j][1], ts[j][2]);
					}
				}
			}
			glEnd();
		}

		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	*/


	/*! draw dual mesh */
	void _draw_dual_cell_face_3D(COMTMesh* pMesh, COMTMeshViewerConfiguration & configure)
	{

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glBegin(GL_TRIANGLES);
		for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
		{
			COMTMesh::CVertex* pv = *viter;
			if (pMesh->_infinite(pv)) continue;
			if (pv->dual_cell3D().edges().empty()) continue;

			std::vector<CPoint> ps;

			for (size_t i = 0; i < pv->dual_cell3D().edges().size(); i++)
			{
				ps.push_back(pv->dual_cell3D().edges()[i].start());
			}

			CPoint n = (configure.m_show_normal_map)?pv->normal_map():pv->normal();
			glNormal3d(n[0], n[1], n[2]);

			glColor3d(pv->rgb()[0], pv->rgb()[1], pv->rgb()[2]);
			
			for (size_t i = 0; i < ps.size() - 2; i++)
			{
				CPoint ts[3];
				ts[0] = ps[0]; ts[1] = ps[i + 1]; ts[2] = ps[i + 2];
				for (int j = 0; j < 3; j++)
				{
					glVertex3d(ts[j][0], ts[j][1], ts[j][2]);
				}
			}

		}
		glEnd();


		if (configure.m_show_backface)
		{
			glBegin(GL_TRIANGLES);

			for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
			{
				COMTMesh::CVertex* pv = *viter;

				if (pMesh->_infinite(pv)) continue;
				if (pv->dual_cell3D().edges().empty()) continue;

				std::vector<CPoint> ps;

				for (size_t i = 0; i < pv->dual_cell3D().edges().size(); i++)
				{
					ps.push_back(pv->dual_cell3D().edges()[i].start());
				}

				CPoint n = (configure.m_show_normal_map)?pv->normal_map():pv->normal();
				glNormal3d(-n[0], -n[1], -n[2]);

				glColor3d(pv->rgb()[0], pv->rgb()[1], pv->rgb()[2]);

				for (size_t i = 0; i < ps.size() - 2; i++)
				{
					CPoint ts[3];
					ts[0] = ps[0]; ts[1] = ps[i + 1]; ts[2] = ps[i + 2];
					for (int j = 2; j >= 0; j--)
					{
						glVertex3d(ts[j][0], ts[j][1], ts[j][2]);
					}
				}
			}
			glEnd();
		}

		glDisable(GL_POLYGON_OFFSET_FILL);
	}


	/*! draw priem mesh */
	void draw_prime_mesh(COMTMesh* pMesh, COMTMeshViewerConfiguration & configure )
	{
		//glDisable(GL_LIGHTING);

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glBindTexture(GL_TEXTURE_2D, configure.m_texture_name);
		glBegin(GL_TRIANGLES);

		for (COMTMesh::MeshFaceIterator fiter(pMesh); !fiter.end(); ++fiter)
		{
			COMTMesh::CFace* pf = *fiter;

			if (pMesh->_infinite(pf)) continue;

			for (COMTMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
			{
				COMTMesh::CVertex* v = *fviter;
				CPoint  pt = v->phi();
				CPoint2 uv = v->uv();

				CPoint n = (configure.m_shade_mode) ? v->normal() : pf->normal();
				if (configure.m_show_normal_map)
				{
					n = v->normal_map();
				}
				glNormal3d(n[0], n[1], n[2]);
				glColor3d(v->rgb()[0], v->rgb()[1], v->rgb()[2]);
				glTexCoord2d(uv[0], uv[1]);

				if (!configure.m_show_uv)
				{
					glVertex3d(pt[0], pt[1], pt[2]);
					continue;
				}
				glVertex3d(pt[0], pt[1], 0);
			}
		}

		if (configure.m_show_backface)
		{
			for (COMTMesh::MeshFaceIterator fiter(pMesh); !fiter.end(); ++fiter)
			{
				COMTMesh::CFace* pf = *fiter;
				if (pMesh->_infinite(pf)) continue;

				std::vector<COMTMesh::CVertex*> vs;

				for (COMTMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
				{
					COMTMesh::CVertex* v = *fviter;
					vs.push_back(v);
				}

				for (int i = 2; i >= 0; i--)
				{
					COMTMesh::CVertex* v = vs[i];
					CPoint pt = v->phi();
					CPoint n = (configure.m_shade_mode) ? v->normal() : pf->normal();
					if (configure.m_show_normal_map)
					{
						n = v->normal_map();
					}
					CPoint2 uv = v->uv();
					glNormal3d(-n[0], -n[1], -n[2]);
					glTexCoord2d(uv[0], uv[1]);
					glColor3d(v->rgb()[0], v->rgb()[1], v->rgb()[2]);
					if (!configure.m_show_uv)
					{
						glVertex3d(pt[0], pt[1], pt[2]);
						continue;
					}
					glVertex3d(pt[0], pt[1], 0);
				}
			}
		}
		glEnd();

		/* Draw here */

		glDisable(GL_POLYGON_OFFSET_FILL);

		if (configure.m_show_edges)
			draw_prime_edges(pMesh, configure);
		if (configure.m_show_boundary)
			draw_boundary(pMesh, configure.m_show_uv);
		if (configure.m_show_sharp_edges)
			draw_sharp_edges(pMesh);
	}


	/*! draw mesh */
	void draw_dual_mesh(COMTMesh* pMesh, COMTMeshViewerConfiguration & configure )
	{
		//draw_vertex(pVertexMaxError, configure.m_show_uv);

		/* Draw here */

		if (configure.m_show_uv)
			_draw_dual_cell_face_2D(pMesh, configure);
		else
			_draw_dual_cell_face_3D(pMesh, configure);

		if (configure.m_show_dual_cell_boundary)
			draw_dual_cell_boundary(pMesh, configure);
		if (configure.m_show_edges)
			draw_dual_edges(pMesh, configure);
		if (configure.m_show_boundary)
			draw_boundary(pMesh, configure.m_show_uv);
		if (configure.m_show_sharp_edges)
			draw_sharp_edges(pMesh);
	}

}//meshlib namespace

#endif
