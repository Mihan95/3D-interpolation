#include <QtGui>
#include "solver.h"
#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

using namespace std;

Solver::Solver(QObject *parent) : QObject(parent)
{
    max_it = 1000;
    eps = 1e-6;
    isFirstStart = true;
    residual = 0;
}

Solver::~Solver()
{
     delete [] a;
     delete [] b;
     delete [] x;
     delete [] I;
     delete [] workspace;
     delete [] args;
}

void Solver::get_command_line_prm( int argc, char* argv[] )
{
    if( argc < 5 )
    {
        cout << "Wrong input parameters\n";
        exit(0);
    }
    N = atoi( argv[1] );
    m = atoi( argv[2] );
    h = atoi( argv[3] );
    p = atoi( argv[4] );
}

void Solver::print_msr_matrix()
{
    for( int i = 0; i < ms; i++ )
    {
        cout << i << " " << I[i] << " " << a[i] << endl;
    }
}

void Solver::print_vector()
{
    for( int i = 0; i < n; i++ )
    {
        cout << i << " " << b[i] << endl;
    }
}

void Jacobi( double *a, int *I, int n, double *y, double *f, int k, int p )
{
    (void) I;
    int i1, i2, i;

    i1  = k * n;
    i1 /= p;
    i2  = (k + 1) * n;
    i2 /= p;

    for( i = i1; i < i2; i++ )
    {
        y[i] = f[i] / a[i];
    }

    pthread_barrier_wait(&barr);
}

double multiply( const double *a, const int *I, const double *x, int i )
{
    double s;
    int len, l, j;
    s   = a[i] * x[i];
    len = I[i + 1] - I[i];
    l   = I[i];
    for( j = 0; j < len; j++ )
        s += a[l + j] * x[I[l + j]];

    return s;
}

void matrix_mult_vector( const double *a, const int *I, int n, const double *x, double *b, int k, int p )
{
    int i, i1, i2;

    i1  = k * n;
    i1 /= p;
    i2  = (k + 1) * n;
    i2 /= p;


    for( i = i1; i < i2; i++ )
    {
        b[i] = multiply( a, I, x, i );
    }

    pthread_barrier_wait( &barr );
}

double scalar_product( double *x, double *y, double *part_s,  int n, int k, int p )
{
    int i1, i2, i;

    part_s[k] = 0;

    i1  = k * n;
    i1 /= p;
    i2  = (k + 1) * n;
    i2 /= p;

    double s = 0;
    for( i =  i1; i < i2; i++ )
        part_s[k] += x[i] * y[i];

    pthread_barrier_wait( &barr );

    for( i = 0; i < p; i++ )
        s += part_s[i];
    pthread_barrier_wait( &barr );
    return s;
}

void linear_comb( double *y, double *z, double c, int n, int k, int p )
{
    int i1, i2, i;

    i1  = k * n;
    i1 /= p;
    i2  = (k + 1) * n;
    i2 /= p;


    for( i = i1; i < i2; i++ )
        y[i] += c * z[i];

    pthread_barrier_wait( &barr );
}

void msr_residual( double *a, int *I, int n, double *x, double *b, double *r, int k, int p )
{
    matrix_mult_vector( a, I, n, x, r, k, p );

    int i1, i2, i;

    i1  = k * n;
    i1 /= p;
    i2  = (k + 1) * n;
    i2 /= p;


    for( i = i1; i < i2; i++ )
        r[i] = b[i] - r[i];

    pthread_barrier_wait( &barr );
}


void *msr_solve(void *ptr_arg)
{
    args_t *ptr = (args_t*) ptr_arg;

    double eps = ptr->eps;
    int max_it = ptr->max_it;
    double *a = ptr->a;
    int *I = ptr->I;
    int n = ptr->n;
    double *b = ptr->b;
    double *x = ptr->x;
    double *workspace = ptr->workspace;
    int k = ptr->k;
    int p = ptr->p;

    double *u = workspace;

    double *v = u + n;
    double *r = v + n;

    double *part_s = r + n;


    int it;
    double c1 = 0.;
    double c2 = 0.;
    double tau;

    msr_residual( a, I, n, x, b, r, k, p );

    for( it = 0; it < max_it; it++ )
    {
        Jacobi( a, I, n, u, r, k, p );

        matrix_mult_vector( a, I, n, u, v, k, p );

        part_s[k] = 0;
        c1 = scalar_product( v, r, part_s, n, k, p );
        c2 = scalar_product( v, v, part_s, n, k, p );
        if( c1 < eps * eps || c2 < eps * eps )
            break;

        tau = c1 / c2;
        linear_comb( x, u,  tau, n, k, p );
        linear_comb( r, v, -tau, n, k, p );
    }

    if( c1 < 0 || c2 < 0 || it >= max_it )
    {
        cout << "c1 == " << c1 << endl;
        cout << "c2 == " << c2 << endl;
        cout << "max_it == " << max_it << endl;
        exit(0);
    }

    if( k == 0 )
        cout << "it == " << it << endl;
    if( k == 0)
    {
        double max = 0;
        for( int i = 0; i < n; i++ )
        {
            if( fabs( r[i] ) > max )
                max = fabs( r[i] );
        }
        cout << "Matrix residual is " << max << endl;
    }

    return NULL;
}

int Solver::get_solution()
{
    ms = matrix_msr_size( m, N );
    J  = h * h * cos(M_PI / N) / (2 * m * m * sin(M_PI / N));
    n  = n_matrix_rows( m, N );

    nTrian = 4 * N * ( m * m );

    I = new int    [ms];
    a = new double [ms];
    b = new double [n];
    x = new double [n];

    t = new triangle[nTrian];

    workspace = new double [4 * n + p];
    args = new args_t [p];

    memset( b, 0, n  * sizeof( double ));
    memset( a, 0, ms * sizeof( double ));
    memset( I, 0, ms * sizeof( int ));
    memset( x, 0, n  * sizeof( double ));
    memset( t, 0,      sizeof( triangle ));
    memset( workspace, 0, (4 * n + p) * sizeof( double ));

    fill_nz( m, a, I, N, J );
    fill_vector( b, J, m, N, h );

    for( int i = 0; i < p; i++ )
    {
        args[i].a = a;
        args[i].b = b;
        args[i].eps = eps;
        args[i].I = I;
        args[i].k = i;
        args[i].max_it = max_it;
        args[i].n = n;
        args[i].p = p;
        args[i].workspace = workspace;
        args[i].x = x;
    }

    for( int k = 1; k < p; k++ )
    {
        if (pthread_create(&tid, 0, &msr_solve, args+k)){
            printf ("Error %d\n", k);
            return 100;
        }
    }
    msr_solve( (void*) args );

    residual =  approximate( x, m, N, h, t );

    cout << "Residual is " << residual << endl;

    residual = 0;

    emit solve_finded( t, nTrian, h );

    return 0;
}

void Solver::more_points()
{
    h *= 2;
    m *= 2;
    t = 0;
    delete [] a;
    delete [] b;
    delete [] x;
    delete [] I;
    delete [] workspace;
    delete [] args;

    cout << "\n\nNumber of points " << m << endl;

    get_solution();
}

void Solver::less_points()
{
    if( m != 1 && h > 1e-15 )
    {
        m /= 2;
        h /= 2;
        t = 0;
        delete [] a;
        delete [] b;
        delete [] x;
        delete [] I;
        delete [] workspace;
        delete [] args;

        cout << "\n\nNumber of points " << m << endl;

        get_solution();
    }
}





















