#ifndef _TRANSFORMATION_H_
#define _TRANSFORMATION_H_

#include <vector>
#include "Geometry/Point.h"
#include "Geometry/Point2.h"

namespace MeshLib
{
	class CTransformation
	{
	public:
		CTransformation() {};
		~CTransformation() {};
		CTransformation(CPoint A, CPoint B, CPoint C)
		{
			CPoint e[3];
			e[2] = (B - A) ^ (C - A);
			e[2] /= e[2].norm();
			e[0] = (B - A);
			e[0] /= e[0].norm();
			e[1] = e[2] ^ e[0];

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					m_a[i][j] = e[i][j];
				}
			m_t = A;
		};

		CPoint operator*(CPoint  v)
		{
			CPoint p = v - m_t;
			CPoint q(0, 0, 0);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					q[i] += m_a[i][j] * p[j];
				}
			return q;
		}
		CPoint operator*(CPoint2 v)
		{
			CPoint p(v[0], v[1], 0);
			return (*this) * p;
		}
		CTransformation inverse()
		{
			CTransformation tr;

			CPoint q(0, 0, 0);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					q[i] += m_a[i][j] * m_t[j];
				}
			tr.m_t = -q;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					tr.m_a[j][i] = m_a[i][j];
				}
			return tr;
		}

	protected:
		CPoint m_t;
		double m_a[3][3];
	};


	class CTransformation2D
	{
	public:
		CTransformation2D() {};
		~CTransformation2D() {};
		//find a transformation from p[i] to q[i]
		
		CTransformation2D(std::vector<CPoint2>& p, std::vector<CPoint2>& q)
		{
			double a[2][2];
			CPoint2 t = p[1] - p[0];
			a[0][0] = t[0]; a[1][0] = t[1];
			t = p[2] - p[0];
			a[0][1] = t[0]; a[1][1] = t[1];
			double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

			double b[2][2];

			b[0][0] = a[1][1] / det; b[1][1] = a[0][0] / det;
			b[0][1] = -a[0][1] / det; b[1][0] = -a[1][0] / det;

			t = q[1] - q[0];
			a[0][0] = t[0]; a[1][0] = t[1];
			t = q[2] - q[0];
			a[0][1] = t[0]; a[1][1] = t[1];

			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
				{
					m_a[i][j] = 0;
					for (int k = 0; k < 2; k++)
						m_a[i][j] += a[i][k] * b[k][j];
				}

			m_t[0] = q[0][0] - (m_a[0][0] * p[0][0] + m_a[0][1] * p[0][1]);
			m_t[1] = q[0][1] - (m_a[1][0] * p[0][0] + m_a[1][1] * p[0][1]);
		};

		CTransformation2D(CPoint2 p[3], CPoint2 q[3])
		{
			double a[2][2];
			CPoint2 t = p[1] - p[0];
			a[0][0] = t[0]; a[1][0] = t[1];
			t = p[2] - p[0];
			a[0][1] = t[0]; a[1][1] = t[1];
			double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

			double b[2][2];

			b[0][0] = a[1][1] / det; b[1][1] = a[0][0] / det;
			b[0][1] = -a[0][1] / det; b[1][0] = -a[1][0] / det;

			t = q[1] - q[0];
			a[0][0] = t[0]; a[1][0] = t[1];
			t = q[2] - q[0];
			a[0][1] = t[0]; a[1][1] = t[1];

			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
				{
					m_a[i][j] = 0;
					for (int k = 0; k < 2; k++)
						m_a[i][j] += a[i][k] * b[k][j];
				}

			m_t[0] = q[0][0] - m_a[0][0] * p[0][0] + m_a[0][1] * p[0][1];
			m_t[1] = q[0][1] - m_a[1][0] * p[0][0] + m_a[1][1] * p[0][1];
		};

		CPoint2 operator*(CPoint2 p)
		{
			CPoint2 q;

			q[0] = m_a[0][0] * p[0] + m_a[0][1] * p[1] + m_t[0];
			q[1] = m_a[1][0] * p[0] + m_a[1][1] * p[1] + m_t[1];

			return q;
		};

		CTransformation2D inverse()
		{
			CTransformation2D tr;

			double det = m_a[0][0] * m_a[1][1] - m_a[1][0] * m_a[0][1];
			
			tr.m_a[0][0] = m_a[1][1] / det; tr.m_a[1][1] = m_a[0][0] / det;
			tr.m_a[0][1] = -m_a[0][1] / det; tr.m_a[1][0] = -m_a[1][0] / det;

			tr.m_t[0] = -tr.m_a[0][0] * m_t[0] - tr.m_a[0][1] * m_t[1];
			tr.m_t[1] = -tr.m_a[1][0] * m_t[0] - tr.m_a[1][1] * m_t[1];

			return tr;
		};

	protected:
		CPoint2 m_t;
		double  m_a[2][2];
	};



}
#endif
