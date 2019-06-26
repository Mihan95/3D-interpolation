#ifndef SCENE3D_H
#define SCENE3D_H

#include <QGLWidget>
#include <QGLBuffer>
#include "funcs.h"

class Scene3D : public QGLWidget
{ 
    Q_OBJECT
   private: 
      GLfloat xRot;
      GLfloat yRot; 
      GLfloat zRot;
      GLfloat zTra;
      GLfloat nSca;
      
      QPoint ptrMousePosition;

      triangle *t;
      triangle *temp_t;
      int nTriang;

      double residual;
      double matrix_res;
      int it;

      int drawMode;
      
      void scale_plus();
      void scale_minus();
      void rotate_up();
      void rotate_down();
      void rotate_left();
      void rotate_right();

      void translate_down();
      void translate_up();     
      void defaultScene();  

      void drawAxis();
        
      void drawFigure();
        
   protected:
      void initializeGL();         
      void resizeGL(int nWidth, int nHeight);  
      void paintGL();                          
      void mousePressEvent(QMouseEvent* pe); 
      void mouseMoveEvent(QMouseEvent* pe);  
      void mouseReleaseEvent(QMouseEvent* pe); 
      void wheelEvent(QWheelEvent* pe);  
      void keyPressEvent(QKeyEvent* pe);       
      
   public: 
      Scene3D(QWidget* parent = 0);
      ~Scene3D();
      void get_triangle( triangle *a );
      void get_nTriang( int n );

public slots:
      void need_paint(triangle *q, int nT, double h);

  signals:
      void press_Mult();
      void press_Div();
      void finished();
}; 
#endif
