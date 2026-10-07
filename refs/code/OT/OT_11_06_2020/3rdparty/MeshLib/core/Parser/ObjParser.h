#ifndef _OBJ_PARSER_H_
#define _OBJ_PARSER_H_

#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>

#include "../Geometry/Point.h"
#include "../Geometry/Point2.h"
#include "../Parser/StrUtil.h"

namespace MeshLib
{


	namespace ObjectParser
	{
		class CMaterial
		{
		public:

			std::string m_name;
			std::string m_map_Kd;

			CPoint  m_Ka;
			CPoint  m_Kd;
			CPoint  m_Ks;
			float   m_Tr;
			float   m_illum;
			float   m_Ns;
			float   m_d;
			
		};

		class CCorner
		{
		public:
			int m_vid;
			int m_tid;
			int m_nid;
			int m_mtl_id;
		};

		class CFace
		{
		public:
			std::vector<CCorner> m_corners;
		};

		class CObjParser
		{
		public:
			CObjParser() {};
			~CObjParser()
			{
				for (size_t i = 0; i < m_mtls.size(); i++)
				{
					delete m_mtls[i];
				}
			}

			//read an object file
			void read_obj(const char* file_name);
			//read material
			void read_mtl(const char* file_name);

		public:

			std::vector<CPoint>  m_vs;
			std::vector<CPoint2> m_vts;
			std::vector<CPoint>  m_vns;
			std::vector<CPoint>  m_rgbs;

			std::vector<CFace>     m_faces;
			std::vector<CMaterial*> m_mtls;

		};

		void CObjParser::read_obj(const char * filename)
		{

			std::fstream f(filename, std::fstream::in);
			if (f.fail())
			{
				std::cerr << "Error in reading file " << filename << std::endl;
				return;
			}

			char cmd[1024];


			bool with_uv = false;
			bool with_normal = false;
			bool with_rgb = false;

			int  mtl_id = 0;

			while (f.getline(cmd, 1024))
			{
				std::string line(cmd);
				line = strutil::trim(line);

				strutil::Tokenizer stokenizer(line, " \t\r\n");

				stokenizer.nextToken();
				std::string token = stokenizer.getToken();

				//mtllib . / sculpt.obj.mtl
				if (token == "mtllib")
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					read_mtl(token.c_str());
				}
				//usemtl material_0

				if (token == "usemtl")
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					for (size_t i = 0; i < m_mtls.size(); i++)
					{
						if (m_mtls[i]->m_name == token)
						{
							mtl_id = (int)i;
							break;
						}
					}
				}
				if (token == "v")
				{
					CPoint p;
					for (int i = 0; i < 3; i++)
					{
						stokenizer.nextToken();
						token = stokenizer.getToken();
						p[i] = strutil::parseString<float>(token);
					}
					m_vs.push_back(p);


					if (stokenizer.nextToken())
					{
						CPoint rgb;
						for (int i = 0; i < 3; i++)
						{
							token = stokenizer.getToken();
							rgb[i] = strutil::parseString<float>(token);
							stokenizer.nextToken();
						}

						m_rgbs.push_back(rgb);
						with_rgb = true;
					}
					continue;
				}


				if (token == "vt")
				{
					with_uv = true;
					CPoint2 uv;
					for (int i = 0; i < 2; i++)
					{
						stokenizer.nextToken();
						token = stokenizer.getToken();
						uv[i] = strutil::parseString<float>(token);
					}
					m_vts.push_back(uv);
					continue;
				}


				if (token == "vn")
				{
					with_normal = true;

					CPoint n;
					for (int i = 0; i < 3; i++)
					{
						stokenizer.nextToken();
						token = stokenizer.getToken();
						n[i] = strutil::parseString<float>(token);
					}
					m_vns.push_back(n);
					continue;
				}


				if (token == "f")
				{
					CFace face;

					for (int i = 0; i < 3; i++)
					{
						CCorner corner;

						stokenizer.nextToken();
						token = stokenizer.getToken();


						strutil::Tokenizer tokenizer(token, " /\t\r\n");

						int ids[3];
						int k = 0;
						while (tokenizer.nextToken())
						{
							std::string token = tokenizer.getToken();
							ids[k] = strutil::parseString<int>(token);
							k++;
						}

						corner.m_mtl_id = mtl_id;

						corner.m_vid = ids[0];
						if (with_uv)
							corner.m_tid = ids[1];
						if (with_normal)
							corner.m_nid = ids[2];

						face.m_corners.push_back(corner);
					}
					m_faces.push_back(face);
				}
			}

			f.close();
		};


		void CObjParser::read_mtl(const char * filename)
		{

			std::fstream f(filename, std::fstream::in);
			if (f.fail())
			{
				std::cerr << "Error in reading file " << filename << std::endl;
				return;
			}

			char cmd[1024];


			CMaterial * pM = NULL;

			while (f.getline(cmd, 1024))
			{

				std::string line(cmd);
				line = strutil::trim(line);

				strutil::Tokenizer stokenizer(line, " \t\r\n");

				stokenizer.nextToken();
				std::string token = stokenizer.getToken();

				

				//newmtl material_0
				if (token == "newmtl")
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					pM = new CMaterial;
					assert(pM != NULL);
					pM->m_name = token;
					m_mtls.push_back(pM);
					continue;
				}

				if (token == "Ka")
				{
					CPoint p;
					for (int i = 0; i < 3; i++)
					{
						stokenizer.nextToken();
						token = stokenizer.getToken();
						p[i] = strutil::parseString<float>(token);
					}
					pM->m_Ka = p;
					continue;
				}

				if (token == "Kd")
				{
					CPoint p;
					for (int i = 0; i < 3; i++)
					{
						stokenizer.nextToken();
						token = stokenizer.getToken();
						p[i] = strutil::parseString<float>(token);
					}
					pM->m_Kd = p;
					continue;
				}

				if (token == "Ks")
				{
					CPoint p;
					for (int i = 0; i < 3; i++)
					{
						stokenizer.nextToken();
						token = stokenizer.getToken();
						p[i] = strutil::parseString<float>(token);
					}
					pM->m_Ks = p;
					continue;
				}

				if (token == "map_Kd")
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					pM->m_map_Kd = token;
					continue;
				}

				if(token == "Tr")
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					pM->m_Tr = strutil::parseString<float>(token);
					continue;
				}


				if( token == "illum" )
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					pM->m_illum = strutil::parseString<float>(token);
					continue;
				}

				if( token == "Ns")
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					pM->m_Ns = strutil::parseString<float>(token);
					continue;
				}
			}

			f.close();
		};

	}//namespace ObjectParser

} //namespace MeshLib

#endif
