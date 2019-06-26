#ifndef FUNCS_H
#define FUNCS_H

#include "pthread.h"

class args_t{
public:
    int k;
    int p;
    int n;
    double *workspace;
    double *a;
    double *b;
    int *I;
    int max_it;
    double eps;
    double *x;
    int it;
    double matrix_res;


    args_t() : k(0), p(0), n(0), workspace(NULL), a(NULL), b(NULL), I(NULL), max_it(0), eps(0), x(NULL), it(0),
        matrix_res(0) {}
};

class point{
public:
    double x;
    double y;

    point() : x(0), y(0) {}
    point( double a, double b ) : x(a), y(b) {}

    friend const point operator+( const point& a, const point& b )
    {
        return point (a.x + b.x, a.y + b.y);
    }
    friend const point operator/( const point& a, const double c )
    {
        return point (a.x / c, a.y / c);
    }
    friend const point operator-( const point& a, const point& b )
    {
        return point (a.x - b.x, a.y - b.y);
    }
    friend const point operator*( const point& a, const double c )
    {
        return point (a.x * c, a.y * c);
    }
};

class point3D{
public:
    double x;
    double y;
    double z;

    point3D() : x(0), y(0), z(0) {}
    point3D( double a, double b, double c ) : x(a), y(b), z(c) {}
    point3D& operator=(const point& a){
        x = a.x;
        y = a.y;
        return *this;
    }

};

class triangle{
public:
    point3D a;
    point3D b;
    point3D c;
    triangle() : a(0, 0, 0), b(0, 0, 0), c(0, 0, 0) {}
};


void  calc_param(char* argv[], int &n, int &p);
void *msr_solve(void *ptr_arg);
void  fill_nz(int m, double *a, int *I, int N , double J);
int   matrix_msr_size(int m, int N );
int   sum( int m, int p, int t );
point polar_to_dec( double r, double fi );
int   n_matrix_rows( int m, int N );
void  fill_vector( double *b, double J, int m, int N, double h/*длина стороны*/ );
double approximate( double *x, int m, int N, double h, triangle *t );
double f( double x, double y );

extern pthread_barrier_t barr;

#endif
