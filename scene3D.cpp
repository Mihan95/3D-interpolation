#include <QtGui>
#include <QGL>
#include <math.h>
#include "scene3D.h"
#include <stdio.h>
#include <iostream>

using namespace std;

void Scene3D::get_triangle( triangle *a )
{
    t = a;
}

void Scene3D::get_nTriang( int n )
{
    nTriang = n;
}

Scene3D::Scene3D(QWidget* parent) : QGLWidget(parent) 
{ 
   xRot=-120; yRot=0; zRot=-45; zTra=0; nSca=1; nTriang = 0;
   t = NULL; drawMode = 1;
} 

Scene3D::~Scene3D()
{
    delete [] t;
}


void Scene3D::initializeGL()
{
   qglClearColor(Qt::white); 
   glEnable(GL_DEPTH_TEST);

   glEnableClientState(GL_VERTEX_ARRAY);
}

void Scene3D::resizeGL(int nWidth, int nHeight)
{ 
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
 
   GLfloat ratio=(GLfloat)nHeight/(GLfloat)nWidth;
   
   if (nWidth>=nHeight)
      glOrtho(-1.5/ratio, 1.5/ratio, -1.5, 1.5, -15.0, 1.5);
   else
      glOrtho(-1.5, 1.5, -1.5*ratio, 1.5*ratio, -15.0, 1.5);
   glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
}

void Scene3D::paintGL()
{ 
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

   glMatrixMode(GL_MODELVIEW); 
   glLoadIdentity();           

   glScalef(nSca, nSca, nSca);
   glTranslatef(0.0f, zTra, 0.0f);    
   glRotatef(xRot, 1.0f, 0.0f, 0.0f);     
   glRotatef(yRot, 0.0f, 1.0f, 0.0f);
   glRotatef(zRot, 0.0f, 0.0f, 1.0f);
 
   drawAxis();
   drawFigure();
}  

void Scene3D::mousePressEvent(QMouseEvent* pe)
{
   ptrMousePosition = pe->pos();                     
} 

void Scene3D::mouseReleaseEvent(QMouseEvent* pe)
{ 
        
} 

void Scene3D::mouseMoveEvent(QMouseEvent* pe) 
{   
   xRot += 180/nSca*(GLfloat)(pe->y()-ptrMousePosition.y())/height(); 
   zRot += 180/nSca*(GLfloat)(pe->x()-ptrMousePosition.x())/width(); 
   
   ptrMousePosition = pe->pos();
   
   updateGL();
}

void Scene3D::wheelEvent(QWheelEvent* pe) 
{ 
   if ((pe->delta())>0) scale_plus(); else if ((pe->delta())<0) scale_minus();   
    
   updateGL();       
}

void Scene3D::keyPressEvent(QKeyEvent* pe) 
{  
   switch (pe->key())
   {
      case Qt::Key_Q:
        emit press_Mult();
      break;

      case Qt::Key_W:
        emit press_Div();
      break;

      case Qt::Key_1:
        drawMode = 1;
      break;

      case Qt::Key_2:
        drawMode = 2;
      break;

      case Qt::Key_3:
        drawMode = 3;
      break;

      case Qt::Key_Plus:
         scale_plus();
      break;
         
      case Qt::Key_Equal:  
         scale_plus(); 
      break;
         
      case Qt::Key_Minus: 
         scale_minus();
      break;

      case Qt::Key_Up:  
         rotate_up();
      break;
         
      case Qt::Key_Down:  
         rotate_down();
      break;         
         
      case Qt::Key_Left:  
        rotate_left();
      break;
         
      case Qt::Key_Right:  
         rotate_right();
      break;                           
         
      case Qt::Key_Z:
         translate_down();
      break;  
         
      case Qt::Key_X:
         translate_up();
      break; 
         
      case Qt::Key_Space:
         defaultScene();
      break;
      
      case Qt::Key_Escape:
         this->close();
         emit finished();
      break;                                                           
   }
   
   updateGL();
}

void Scene3D::scale_plus()
{
   nSca = nSca*1.1;
}

void Scene3D::scale_minus()
{
   nSca = nSca/1.1;
}

void Scene3D::rotate_up()
{
   xRot += 1.0;
}

void Scene3D::rotate_down()
{
   xRot -= 1.0;
}

void Scene3D::rotate_left()
{
   zRot += 5.0;
}

void Scene3D::rotate_right()
{
   zRot -= 5.0;
}

void Scene3D::translate_down()
{
   zTra -= 0.05;
}

void Scene3D::translate_up()
{
   zTra += 0.05;
}

void Scene3D::defaultScene()
{
   xRot=-120; yRot=0; zRot=-45; zTra=0; nSca=1;
}

void Scene3D::drawAxis()
{
   glLineWidth(3.0f); 
   
   glColor4f(1.00f, 0.00f, 0.00f, 1.0f); 
   glBegin(GL_LINES);
      glVertex3f( 2000.0f / nSca,  0.0f,  0.0f);
      glVertex3f(-2000.0f / nSca,  0.0f,  0.0f);
   glEnd();  
   
   QColor halfGreen(0, 128, 0, 255);
   qglColor(halfGreen);
   glBegin(GL_LINES);
      glVertex3f( 0.0f,  2000.0f / nSca,  0.0f);
      glVertex3f( 0.0f, -2000.0f / nSca,  0.0f);
  
      glColor4f(0.00f, 0.00f, 1.00f, 1.0f);
      glVertex3f( 0.0f,  0.0f,  2000.0f / nSca);
      glVertex3f( 0.0f,  0.0f, -2000.0f / nSca);
   glEnd();
}

void Scene3D::drawFigure()
{
    glLineWidth(1);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    if( drawMode == 1 )
    {
        glColor3d(1,0,0);
        glVertexPointer(3, GL_DOUBLE, 0, t);
        glDrawArrays(GL_TRIANGLES , 0, 3*nTriang);
    }
    if( drawMode == 2 )
    {
        glColor3d(1,0,0);
        glVertexPointer(3, GL_DOUBLE, 0, t);
        glDrawArrays(GL_TRIANGLES , 0, 3*nTriang);

        glBegin(GL_TRIANGLES);
          glColor3d(0,1,0);
         for( int i = 0; i < nTriang; i += 4 )
        {
               glVertex3d( t[i].a.x, t[i].a.y, f( t[i].a.x, t[i].a.y ));
               glVertex3d( t[i].b.x, t[i].b.y, f( t[i].b.x, t[i].b.y ));
               glVertex3d( t[i].c.x, t[i].c.y, f( t[i].c.x, t[i].c.y ));


               glVertex3d( t[i + 1].a.x, t[i + 1].a.y, f( t[i + 1].a.x, t[i + 1].a.y ));
               glVertex3d( t[i + 1].b.x, t[i + 1].b.y, f( t[i + 1].b.x, t[i + 1].b.y ));
               glVertex3d( t[i + 1].c.x, t[i + 1].c.y, f( t[i + 1].c.x, t[i + 1].c.y ) );

               glVertex3d( t[i + 2].a.x, t[i + 2].a.y, f( t[i + 2].a.x, t[i + 2].a.y ));
               glVertex3d( t[i + 2].b.x, t[i + 2].b.y, f( t[i + 2].b.x, t[i + 2].b.y ));
               glVertex3d( t[i + 2].c.x, t[i + 2].c.y, f( t[i + 2].c.x, t[i + 2].c.y ));

               glVertex3d( t[i + 3].a.x, t[i + 3].a.y, f( t[i + 3].a.x, t[i + 3].a.y ));
               glVertex3d( t[i + 3].b.x, t[i + 3].b.y, f( t[i + 3].b.x, t[i + 3].b.y ));
               glVertex3d( t[i + 3].c.x, t[i + 3].c.y, f( t[i + 3].c.x, t[i + 3].c.y ));
         }

         glEnd();
    }
    if( drawMode == 3 )
    {
        glBegin(GL_TRIANGLES);
          glColor3d(0,1,0);
         for( int i = 0; i < nTriang; i += 4 )
        {
               glVertex3d( t[i].a.x, t[i].a.y, t[i].a.z );
               glVertex3d( t[i].b.x, t[i].b.y, t[i].b.z );
               glVertex3d( t[i].c.x, t[i].c.y, t[i].c.z );


               glVertex3d( t[i + 1].a.x, t[i + 1].a.y, t[i + 1].a.z );
               glVertex3d( t[i + 1].b.x, t[i + 1].b.y, t[i + 1].b.z );
               glVertex3d( t[i + 1].c.x, t[i + 1].c.y, t[i + 1].c.z );

               glVertex3d( t[i + 2].a.x, t[i + 2].a.y, t[i + 2].a.z );
               glVertex3d( t[i + 2].b.x, t[i + 2].b.y, t[i + 2].b.z );
               glVertex3d( t[i + 2].c.x, t[i + 2].c.y, t[i + 2].c.z );

               glVertex3d( t[i + 3].a.x, t[i + 3].a.y, t[i + 3].a.z );
               glVertex3d( t[i + 3].b.x, t[i + 3].b.y, t[i + 3].b.z );
               glVertex3d( t[i + 3].c.x, t[i + 3].c.y, t[i + 3].c.z );
         }

         glEnd();
    }

}

void Scene3D::need_paint( triangle *q, int nT , double h )
{
    delete [] t;
    t = q;
    nTriang = nT;
    nSca = 1/h;

    updateGL();
}
