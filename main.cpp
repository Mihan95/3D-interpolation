#include <QApplication>
#include <iostream>
#include <math.h>
#include "pthread.h"
#include "scene3D.h"
#include "solver.h"
#include "solver_thread.h"

using namespace std;

pthread_barrier_t barr;

int main(int argc, char* argv[])
{
    QApplication app( argc, argv );

    Scene3D       *scene1 = new Scene3D;
    Solver       *solver1 = new Solver;
    Solve_thread  *thread = new Solve_thread;

    QObject::connect(scene1,  SIGNAL(press_Mult()), solver1, SLOT(more_points()));
    QObject::connect(scene1,  SIGNAL(press_Div()),  solver1, SLOT(less_points()));
    QObject::connect(solver1, SIGNAL(solve_finded(triangle*, int, double)), scene1,
                     SLOT(need_paint(triangle*, int, double)));
    QObject::connect(scene1,  SIGNAL(finished()), thread, SLOT(quit()));

    solver1->get_command_line_prm( argc, argv );
    pthread_barrier_init( &barr, NULL, solver1->p );

    solver1->moveToThread( thread );
    thread->start();
    solver1->get_solution();

    scene1->resize( 500, 500 );
    scene1->show();

    int r = app.exec();
    delete scene1;
    delete solver1;
    delete thread;
    return r;
}



















