#ifndef SOLVER_H
#define SOLVER_H

#include <QWidget>
#include <QMainWindow>
#include <QEvent>
#include "funcs.h"

class Solver : public QObject
{
    Q_OBJECT

private:

    pthread_t tid;

    int N;
    int m;
    int ms;
    int max_it;
    int it;

    int n;

    double h;
    double J;
    double eps;
    double residual;
    double matrix_res;

    int *I;

    double *a;
    double *b;
    double *x;
    double *workspace;

    args_t *args;

    bool isFirstStart;

    QKeyEvent *pe;

public:

    int nTrian;
    int p;
    triangle *t;

    Solver( QObject* parent = 0 );
    ~Solver();
    void get_command_line_prm( int argc, char* argv[] );
    int get_solution();
    void print_msr_matrix();
    void print_vector();

public slots:
    void more_points();
    void less_points();

signals:
    void solve_finded(triangle *q, int nT, double h);
};

#endif
