/*****************************************************************************************
*      Optimal Transportation Map based on Gu-Luo-Sun-Yau theorem
*
*
*    Purpose:
* 
*       Compute the Optimal Transportation Map between probability measures defined 
*		on planar domains
*
*       David Xianfeng Gu October 12, 2020
*
*	    Computer Science Department
*	    Stony Brook University 
*		Stony Brook, New York 11794
*		gu@cs.stonybrook.edu
* 
*		"Variational Principles for Minkowski Type Problems, Discrete Optimal Transport, 
*		and Discrete Monge-Ampere Equations"
*		Xianfeng Gu, Feng Luo, Jian Sun and Shing-Tung Yau, AJM, 20(2), 383-398 
* 
*		This code is only for education and research purpose. For any potential commercial
*		application, please contact the author.
* 
*		All rights are preserved.
*	
*****************************************************************************************/


/*! \file omt_viewer.cpp
*   \brief Brenier/Alexandrov/Gu-Yau Theorem
*   \author David Gu
*   \date documented on 10/12/2020
*
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <Eigen/Eigen>

#ifdef MAC_OS
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif // MAC_OS

#include "viewer/Arcball.h"         /*  Arc Ball  Interface         */
#include "bmp/RgbImage.h"
#include "Operator/Operator.h"
#include "OptimalMassTransportationMesh.h"
#include "Geometry/Segment2D.h"
#include "COMTMeshViewer.h"

using namespace MeshLib;

/* window width and height */
int win_width, win_height;
int gButton;
int startx,starty;

COMTMeshViewerConfiguration configure;

/* rotation quaternion and translation vector for the object */
CQrot       ObjRot(0,0,1,0);
CPoint      ObjTrans(0,0,0);
CQrot       LightRot(1,0,0,0);
/* light position */
CPoint LightPosition;
/* arcball object */
CArcball arcball;
/* texture image */
RgbImage image;

/* the source domain */
CPolygon2D  D;

/* global target mesh */
COMTMesh mesh;

/* global step length */
double step_length = 5e-2;

//copy frame buffer to an image
/*! save frame buffer to an image "snap_k.bmp"
*/
void read_frame_buffer()
{
	static int id = 0;

	GLfloat * buffer = new GLfloat[win_width * win_height*3];
	assert( buffer );
	glReadBuffer(GL_FRONT_LEFT);
	glReadPixels( 0, 0, win_width, win_height, GL_RGB, GL_FLOAT, buffer);

	RgbImage  image(win_height,win_width);

	for( int i = 0; i < win_height; i ++ )
	for( int j = 0; j < win_width;j ++ )
	{
		float r = buffer[(i*win_width + j)*3+0];
		float g = buffer[(i*win_width + j)*3+1];
		float b = buffer[(i*win_width + j)*3+2];

		image.SetRgbPixelf( i,j, r,g,b);
	}
	delete []buffer;

	char name[256];
	std::ostringstream os(name);
	os << "snape_"<< id++ << ".bmp";
	image.WriteBmpFile( os.str().c_str() );

}


/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
    double rot[16];

    glTranslated( ObjTrans[0], ObjTrans[1], ObjTrans[2]);
    ObjRot.convert( rot );
    glMultMatrixd(( GLdouble *)  rot );	
}

/*! the eye is always fixed at world z = +5 */

void setupEye(void){
  glLoadIdentity();
  gluLookAt( 0,0, 5,0,0,0,0,1,0);
  
}

/*! setup light */
void setupLight()
{
	CQrot v( 0, 0,0, -1 );
	CQrot CL( LightRot.m_w, -LightRot.m_x, -LightRot.m_y, -LightRot.m_z);
	CL = LightRot*v*CL;
	LightPosition = CPoint( CL.m_x, CL.m_y, CL.m_z );
	GLfloat lightOnePosition[4]={(GLfloat)LightPosition[0], (GLfloat)LightPosition[1], (GLfloat)LightPosition[2], 0.0f};
	glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);

}

/*! display call back function
*/
void display()
{
	/* clear frame buffer */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	  
	setupLight();
	/* transform from the eye coordinate system to the world system */
	setupEye();
	glPushMatrix();
	/* transform from the world to the ojbect coordinate system */
	setupObject();

	//draw_morph_prime_mesh(&mesh, configure);
	
  /* draw the mesh */

	COMTMesh* pMesh = &mesh;


	switch (configure.m_prime_dual_view)
	{
	case 0:
		if (pMesh != NULL)
		{
			draw_prime_mesh(pMesh, configure);
		}
		break;
	case 1:
		if (pMesh != NULL)
		{
			draw_dual_mesh(pMesh, configure );
		}
		break;
	case 2:
		if (pMesh != NULL)
		{
			draw_prime_mesh(pMesh, configure);
			draw_dual_mesh(pMesh,configure);
		}
		break;
	}

	if (configure.m_show_cell_center)
	{
		draw_cell_center(pMesh);
	}

	glPopMatrix();
	glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */

void reshape(int w, int h)
{
  float ar;

  win_width=w;
  win_height=h;

  ar = (float)(w)/h;
  glViewport(0, 0, w, h);               /* Set Viewport */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); 

  // magic imageing commands
  gluPerspective( 40.0, /* field of view in degrees */
		  ar, /* aspect ratio */
		  0.0001, /* Z near */ 
		  100.0 /* Z far */);

  glMatrixMode(GL_MODELVIEW);
  
  glutPostRedisplay();
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
	printf("-----------------------------------------------------------------\n");
	printf(" Demo for Alexandrov Theorem \n");
	printf(" Input  : Genus zero surface with a single boundary, with UV coordinages on the unit disk\n");
	printf(" Output : Optimal transportation map between the surface and the plan area element \n");
	printf(" Command: Viewer.exe -target target_mesh -source source_mesh\n");
	printf(" Debug  : Press '!' for one step Newton's method\n");
	printf("		  Press 'm' for Calculating the dual cell center\n");
	printf("          Press 'W' for output to 'output_phi.m' and 'output_ot.m' with Legendre dual and ot map\n");
	printf(" Display: prime view and dual view\n");
	printf("-----------------------------------------------------------------\n");
	printf(" d  -  Configure: show prime mesh, dual mesh or both\n");
	printf(" L  -  Configure: edit light source position\n");
	printf(" l  -  Configure: show dual cell edges\n");
	printf(" B  -  Configure: Show back faces\n");
	printf(" c  -  Configure: show cell center\n");
	printf(" b  -  Configure: Show boundary edges\n");
	printf(" e  -  Configure: Show prime edges and dual edges\n");
	printf(" n  -  Configure: Show original surface normal map\n");
	printf(" m  -  Calculate: dual cell center\n");
	printf("------------------------------------------------------------------\n");
	printf(" w  -  Wireframe Display\n");
	printf(" f  -  Flat Shading \n");
	printf(" s  -  Smooth Shading\n");
	printf(" t  -  Toggle texture Mapping\n");
	printf(" o  -  Save frame buffer to snap_n.bmp\n");
	printf(" ?  -  Help Information\n");
	printf(" esc - quit\n");
	printf("------------------------------------------------------------------\n");
}

void key_process( unsigned char key )
{
	switch( key )
	{

	case 'f':
		//Flat Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		configure.m_shade_mode = 0;
		break;
	case 's':
    //Smooth Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		configure.m_shade_mode = 1;
		break;
	case 'w':
	  //Wireframe mode
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
  case 't':
    configure.m_texture_mode = (configure.m_texture_mode +1 )%3;
    switch(configure.m_texture_mode)
    {
    case 0:
        glDisable(GL_TEXTURE_2D);
        break;
    case 1:
        glEnable(GL_TEXTURE_2D);
		    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        break;
    case 2:
       glEnable(GL_TEXTURE_2D);
		    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
       break;
    }
    break;
  case '?':
    help();
    break;

  case 'o':
	  read_frame_buffer();
	  break;

  case 27:
		exit(0);
		break;

  case 'e':
	configure.m_show_edges = ! configure.m_show_edges;
	break;

  case 'b':
	configure.m_show_boundary = ! configure.m_show_boundary;
	break;

  case 'S':
	configure.m_show_sharp_edges = ! configure.m_show_sharp_edges;
	break;

  case 'B':
	  configure.m_show_backface = ! configure.m_show_backface;
	  break;
  case 'g':
	  configure.m_show_uv = !configure.m_show_uv;
	  break;
  case 'n':
	  configure.m_show_normal_map = !configure.m_show_normal_map;
	  break;
  case 'd':
	  configure.m_prime_dual_view = (configure.m_prime_dual_view+1)%3;
	  break;
	  break;
  case 'c':
	  configure.m_show_cell_center = !configure.m_show_cell_center;
	  break;
  case 'L':
	  configure.m_light_edit = !configure.m_light_edit;
	  break;
  case 'l':
	  configure.m_show_dual_cell_boundary = !configure.m_show_dual_cell_boundary;
	  break;
  case '!':
		OT_Newton(mesh, D, step_length, configure.m_nearest);
	  break;
  case 'm':
	  mesh._power_cell_center();
	  break;
  case 'W':
	{	  
		std::cout << "Output to 'output_ot.m' and 'output_phi.m' " << std::endl;
		mesh.write_m("output_phi.m");
		mesh.write_m("output_ot.m");
	}
	 break;
	}
}
/*! Keyboard call back function */

void keyBoard(unsigned char key, int x, int y) 
{
	key_process(key );
	glutPostRedisplay();  
}

/*! setup GL states */
void setupGLstate(){

  GLfloat lightOneColor[] = {1, 1, 1, 1};
  GLfloat globalAmb[] = {.1f, .1f, .1f, 1.0f};
  GLfloat lightOnePosition[] = {.0f,  .0f, 1.0f, 0.0f};

  GLfloat lightTwoColor[]    = { 1, 1, 1, 1 };
  GLfloat lightTwoPosition[] = { .0f,  .0f, 1.0f, 0.0f };

  glEnable(GL_CULL_FACE);
  glFrontFace(GL_CCW);      
  //glFrontFace(GL_CW);      
  glEnable(GL_DEPTH_TEST);
  //glClearColor(0,0,0,0);
  glClearColor(1.0,1.0,1.0,1.0);
  glShadeModel(GL_SMOOTH);


  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);

  glEnable(GL_LIGHT2);
  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, lightTwoColor);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
  glDisable(GL_LIGHT2);
  
	const GLfloat specular[]= {1.0f,1.0f,1.0f,1.0f};

	glLightfv(GL_LIGHT1, GL_SPECULAR, specular );
	
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 64.0f );

	GLfloat mat_ambient[] = { 0.0f , 0.0f , 0.0f , 1.0f };
	GLfloat mat_diffuse[] = { 0.01f , 0.01f , 0.01f , 1.0f };
	GLfloat mat_specular[] = {0.5f , 0.5f , 0.5f, 1.0f };
	GLfloat mat_shininess[] = { 32 };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

}

/*! mouse click call back function */
void  mouseClick(int button , int state, int x, int y){

  
  /* set up an arcball around the Eye's center
	  switch y coordinates to right handed system  */

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) 
  { 
      gButton = GLUT_LEFT_BUTTON;  
	  arcball = CArcball( win_width, win_height,  x-win_width/2,  win_height-y-win_height/2); 
  }
  
  if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
      startx = x;
      starty = y;
      gButton = GLUT_MIDDLE_BUTTON;
   }
  
  if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
      startx = x;
      starty = y;
      gButton = GLUT_RIGHT_BUTTON;
   }
  return ;
}

/*! mouse motion call back function */

void mouseMove(int x, int y)
{
  CPoint trans;
  CQrot       rot;
  
  /* rotation, call arcball */
  if (gButton == GLUT_LEFT_BUTTON ) 
  {
      rot = arcball.update( x-win_width/2, win_height-y-win_height/2);
	  if( configure.m_light_edit )
		  LightRot = rot * LightRot;
	  else
	      ObjRot =  rot * ObjRot;
      glutPostRedisplay();
  }
  
  /*xy translation */
  if (gButton == GLUT_MIDDLE_BUTTON) 
  {
	  double scale = 10./win_height;
      trans =  CPoint(scale*(x-startx), scale*(starty-y), 0  );
	    startx = x;
	    starty = y;
      ObjTrans = ObjTrans + trans;
      glutPostRedisplay();
  }
  
  /* zoom in and out */
  if (gButton == GLUT_RIGHT_BUTTON ) {
      double scale = 10./win_height;
      trans =  CPoint(0,0, scale*(starty-y)   );
	    startx = x;
	    starty = y;
      ObjTrans = ObjTrans+trans;
      glutPostRedisplay();
  }
  
}


/*! initialize bitmap image texture */

void initialize_bmp_texture()
{
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &configure.m_texture_name);
	glBindTexture(GL_TEXTURE_2D, configure.m_texture_name);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,   GL_LINEAR);

	int ImageWidth  = image.GetNumCols();
	int ImageHeight = image.GetNumRows();
	GLubyte * ptr   = (GLubyte * )image.ImageData();

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 
				        ImageWidth,
						ImageHeight, 
				        0, 
				        GL_RGB, 
				        GL_UNSIGNED_BYTE,
						ptr);

    if(configure.m_texture_mode == 1)
		  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    else if(configure.m_texture_mode == 2)
		  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_2D);
}

void process_key( char * keys )
{
	for( size_t i = 0; i < strlen( keys ); i ++ )
	{
		unsigned char ch = keys[i];
		key_process( ch );
	}
}


/*! main function for viewer
*/
int main( int argc, char * argv[] )
{

	if( argc < 3 ) 
	{
		printf("Usage: %s -nearest -source mesh_name -target target_mesh -step_length lambda\n", argv[0]);
		return -1;
	}

	std::string name;

	for( int i = 1; i < argc; i ++ )
	{
		if( strcmp( argv[i], "-target" ) == 0 )
		{

			std::string mesh_name( argv[++i] );
			name = mesh_name;
			
			if( strutil::endsWith( mesh_name, ".obj" ) )
			{
				mesh.read_obj( mesh_name.c_str() );
				for (COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
				{
					COMTMesh::CVertex* pv = *viter;
					COMTMesh::CHalfEdge* ph = mesh.vertexHalfedge(pv);
					pv->uv() = ph->uv();
				}
			}

			if( strutil::endsWith( mesh_name, ".m" ) )
			{
				mesh.read_m( mesh_name.c_str() );
			}
		
			COperator<COMTMesh> pS( &mesh );
			pS._normalize();
			pS._calculate_face_vertex_normal();
	
			//copy original surface vertex normal to normal map
			for (COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
			{
				COMTMesh::CVertex* pV = *viter;
				pV->normal_map() = pV->normal();
			}

			//process color
			bool has_rgb = false;
			for( COMTMesh::MeshVertexIterator viter( &mesh ); !viter.end(); viter ++ )
			{
				COMTMesh::CVertex * pV = *viter;
				CPoint dp = pV->rgb() - CPoint(1,1,1);
				if( dp.norm() > 0 ) has_rgb = true;
			}
			if( !has_rgb )
			{
				for( COMTMesh::MeshVertexIterator viter( &mesh ); !viter.end(); viter ++ )
				{
					COMTMesh::CVertex * pV = *viter;
					pV->rgb() = CPoint( 229.0/255.0, 162.0/255.0, 141.0/255.0 );
				}
			}
		}

		if (strcmp(argv[i], "-source") == 0)
		{
			COMTMesh source;
			std::string mesh_name(argv[++i]);
			name = mesh_name;

			if (strutil::endsWith(mesh_name, ".obj"))
			{
				source.read_obj(mesh_name.c_str());
				for (COMTMesh::MeshVertexIterator viter(&source); !viter.end(); viter++)
				{
					COMTMesh::CVertex* pv = *viter;
					COMTMesh::CHalfEdge* ph = source.vertexHalfedge(pv);
					pv->uv() = ph->uv();
				}
			}

			if (strutil::endsWith(mesh_name, ".m"))
			{
				source.read_m(mesh_name.c_str());
			}

			COMTMesh::CBoundary bnd(&source);
			COMTMesh::CLoop * pL = bnd.loops().front();
			std::vector<CPoint2> pts;

			for (std::list<COMTMesh::CHalfEdge*>::iterator iter = pL->halfedges().begin(); iter != pL->halfedges().end(); iter++)
			{
				COMTMesh::CHalfEdge* ph = *iter;
				COMTMesh::CVertex* pv = mesh.halfedgeTarget(ph);
				pts.push_back(CPoint2(pv->point()[0], pv->point()[1])*1.2);
			}
			for (int i = 0; i < pts.size(); i++)
			{
				CSegment2D s(pts[i], pts[(i + 1) % pts.size()]);
				D.add(s);
			}
		}


		if (strcmp(argv[i], "-mesh2obj") == 0)
		{
			COMTMesh mesh;
			mesh.read_m(argv[++i]);
			
			COperator<COMTMesh> pS(&mesh);
			pS._normalize();
			pS._calculate_face_vertex_normal();

			mesh.write_obj(argv[++i]);
			return 0;
		}

		if (strcmp(argv[i], "-step_length") == 0)
		{
			step_length = atof(argv[++i]);
		}

		if (strcmp(argv[i], "-nearest") == 0)
		{
			configure.m_nearest = true;
			std::cout << "Optimal Transportation Plan" << std::endl;
		}
		
		if (strcmp(argv[i], "-farthest") == 0)
		{
			configure.m_nearest = false;
			std::cout << "Worst Transportation Plan" << std::endl;
		}

		if( strcmp( argv[i], "-key" ) == 0 )
		{
			process_key( argv[i+1] );
		}
	}

	//initialize the OT Map
	OT_Initialize<COMTMesh>(mesh, D, configure.m_nearest);

	/* glut stuff */
	glutInit(&argc, argv);                /* Initialize GLUT */
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize(600, 600);
	std::string title("Optimal Transportation Map Viewer - ");
	title = title + name;
	glutCreateWindow(title.c_str());        /* Create window with given title */
	glViewport(0,0,800,800 );

	glutDisplayFunc(display);             /* Set-up callback functions */
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyBoard);
	setupGLstate();
	initialize_bmp_texture();

	glutMainLoop();                       /* Start GLUT event-processing loop */

	return 0;
}