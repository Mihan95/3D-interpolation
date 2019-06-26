#include "funcs.h"
#include <math.h>
#include <iostream>

using namespace std;


void calculate_triangle( triangle *t, double *x, int &q, int i, int j, int k,
                       point pA, point p0, point p1, point pA0, point pA1, point p01, double  &residual )
{
    double temp = 0.0;

    t[q].a.x = pA.x;  t[q].a.y = pA.y;  t[q].a.z =  x[i];
    t[q].b.x = pA0.x; t[q].b.y = pA0.y; t[q].b.z = (x[i] + x[j]) / 2;
    t[q].c.x = pA1.x; t[q].c.y = pA1.y; t[q].c.z = (x[i] + x[k]) / 2;


    temp = fabs( f( t[q].a.x, t[q].a.y ) - t[q].a.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].b.x, t[q].b.y ) - t[q].b.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].c.x, t[q].c.y ) - t[q].c.z );
    if( temp > residual )
        residual = temp;
    q++;

    t[q].a.x = p0.x;  t[q].a.y = p0.y;  t[q].a.z =  x[j];
    t[q].b.x = pA0.x; t[q].b.y = pA0.y; t[q].b.z = (x[i] + x[j]) / 2;
    t[q].c.x = p01.x; t[q].c.y = p01.y; t[q].c.z = (x[j] + x[k]) / 2;


    temp = fabs( f( t[q].a.x, t[q].a.y ) - t[q].a.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].b.x, t[q].b.y ) - t[q].b.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].c.x, t[q].c.y ) - t[q].c.z );
    if( temp > residual )
        residual = temp;
    q++;

    t[q].a.x = p1.x;  t[q].a.y = p1.y;  t[q].a.z =  x[k];
    t[q].b.x = p01.x; t[q].b.y = p01.y; t[q].b.z = (x[j] + x[k]) / 2;
    t[q].c.x = pA1.x; t[q].c.y = pA1.y; t[q].c.z = (x[i] + x[k]) / 2;


    temp = fabs( f( t[q].a.x, t[q].a.y ) - t[q].a.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].b.x, t[q].b.y ) - t[q].b.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].c.x, t[q].c.y ) - t[q].c.z );
    if( temp > residual )
        residual = temp;
    q++;

    t[q].a.x = pA0.x; t[q].a.y = pA0.y; t[q].a.z = (x[i] + x[j]) / 2;
    t[q].b.x = p01.x; t[q].b.y = p01.y; t[q].b.z = (x[j] + x[k]) / 2;
    t[q].c.x = pA1.x; t[q].c.y = pA1.y; t[q].c.z = (x[i] + x[k]) / 2;
    temp = fabs( f( t[q].a.x, t[q].a.y ) - t[q].a.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].b.x, t[q].b.y ) - t[q].b.z );
    if( temp > residual )
        residual = temp;
    temp = fabs( f( t[q].c.x, t[q].c.y ) - t[q].c.z );
    if( temp > residual )
        residual = temp;
    q++;


}

void fill_triangles_III( double *x, int m, int N, double r, double alf, int l, triangle *t, point shift_x, point shift_y,
                         int &q, double &residual )
{
    int i, j = 2;
    int i1 = 2 + l * m * ( m + 1 ) / 2;
    int i2 = m - 1 + l * m * ( m + 1 ) / 2;
    int nElems = N * m * ( m + 1 ) / 2;

    point pA, p0, p1, pA0, pA1, p01;
    point c(0, 0);
    point a = polar_to_dec( r,     M_PI - l * alf );
    point d = polar_to_dec( j * r, M_PI - l * alf );

    int i3 = 1 + l * m * ( m + 1 ) / 2;
    int k;
    pA = a;

    p0 = c;     p1 = c + shift_y;
    pA0  = ( pA + p0 ) / 2;    pA1 = ( pA + p1 ) / 2;    p01 = ( p0 + p1 ) / 2;
    k = ( l + 1 ) * m * ( m + 1 ) / 2 + 1;
    if ( k > nElems ) k -= nElems;
    calculate_triangle( t, x, q, i3, 0, k, pA, p0, p1, pA0, pA1, p01, residual );

    p0 = p1;    p1 = p0 + shift_x;
    pA0  = ( pA + p0 ) / 2;    pA1 = ( pA + p1 ) / 2;    p01 = ( p0 + p1 ) / 2;
    calculate_triangle( t, x, q, i3, k, i3 + m, pA, p0, p1, pA0, pA1, p01, residual );

    j++;
    c = a; a = d; d = polar_to_dec( j * r, M_PI - l * alf );

    for( i = i1; i <= i2; i++ )
    {
        pA  = a;

        p0   = c;        p1  = c + shift_y;
        pA0  = ( pA + p0 ) / 2;    pA1 = ( pA + p1 ) / 2;    p01 = ( p0 + p1 ) / 2;
        calculate_triangle( t, x, q, i, i - 1, i + m - 1, pA, p0, p1, pA0, pA1, p01, residual );

        p0   = p1;       p1 = p0 + shift_x;
        pA0  = ( pA + p0 ) / 2;    pA1 = ( pA + p1 ) / 2;    p01 = ( p0 + p1 ) / 2;
        calculate_triangle( t, x, q, i, i + m - 1, i + m, pA, p0, p1, pA0, pA1, p01, residual );

        j++;
        c = a; a = d; d = polar_to_dec( j * r, M_PI - l * alf );
    }
}

void fill_triangles_V( double *x, int m, double r, double alf, int l, triangle *t, point shift_x, point shift_y, int &q, double &residual )
{
    int i = m + l * m * ( m + 1 ) / 2;
    point pA, p0, p1, pA0, pA1, p01;
    point c = polar_to_dec( r * (m - 1), M_PI - l * alf );
    point a = polar_to_dec( r *  m,      M_PI - l * alf );

    (void) shift_x;

    pA = a;
    p0 = c;       p1 = c + shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    calculate_triangle( t, x, q, i, i - 1, i + m - 1, pA, p0, p1, pA0, pA1, p01, residual );
}

void fill_triangles_VI( double *x, int m, int N, int p, double r, double alf, int l, triangle *t, point &a,
                        point shift_x, point shift_y, int &q, double &residual )
{
    // p = m; p >= 3; p--
    int i = 1 + sum( m, p, 0 ) + l * m * ( m + 1 ) / 2;
    int j, k;

    point pA, p0, p1, pA0, pA1, p01;

    a = polar_to_dec( r, M_PI - l * alf );
    a = a + shift_y * (m - p + 1);

    pA = a;
    p0 = pA - shift_x;    p1 = p0 + shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    if( l != N - 1 )
    {

        j = m - p + 1 + (l+1) * m * ( m + 1 ) / 2;
        k = j + 1;

    }
    else
    {
        j =  m - p + 1;
        k = j + 1;
    }
    calculate_triangle( t, x, q, i, j, k, pA, p0, p1, pA0, pA1, p01, residual );

    p0 = p1;    p1 = p1 + shift_x;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    calculate_triangle( t, x, q, i, k, i + p - 1, pA, p0, p1, pA0, pA1, p01, residual );

    a = a + shift_x;
}

void fill_triangles_II( double *x, int m, int p, int l, triangle *t, point &a, point shift_x, point shift_y, int &q, double &residual )
{   //p = m; p >= 4; p--

    int i, i1, i2;
    point pA, p0, p1, pA0, pA1, p01;

    i1 = 2 + sum( m, p, 0 ) + l * m * ( m + 1 ) / 2;
    i2 = i1 + p - 4;

    for( i = i1; i <= i2; i++ )
    {
        pA = a;
        p0 = pA - shift_x;    p1 = p0 + shift_y;
        pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
        calculate_triangle( t, x, q, i, i - 1, i + p - 2, pA, p0, p1, pA0, pA1, p01, residual );

        p0 = p1;    p1 = p1 + shift_x;
        pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
        calculate_triangle( t, x, q, i, i + p - 2, i + p - 1, pA, p0, p1, pA0, pA1, p01, residual );

        a = a + shift_x;
    }
}

void fill_triangles_IV( double *x, int m, int p, int l, triangle *t, point &a, point shift_x, point shift_y, int &q, double &residual )
{
    //p = m; p >= 3; p--
    int i = m + sum( m, p, 1 ) + l * m * ( m + 1 ) / 2;
    point pA, p0, p1, pA0, pA1, p01;

    pA = a;
    p0 = pA - shift_x;    p1 = pA - shift_x + shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    calculate_triangle( t, x, q, i, i - 1, i + p - 2, pA, p0, p1, pA0, pA1, p01, residual );
}

void fill_triangles_VII( double *x, int m, int N, double r, double alf, int l, triangle *t,
                         point shift_x, point shift_y, int &q, double &residual )
{
    int i = ( l + 1 ) * m * ( m + 1 ) / 2;
    int j, k;
    point pA, p0, p1, pA0, pA1, p01;

    point a = polar_to_dec( r, M_PI - l * alf );
    a = a + shift_y * (m - 1);

    pA = a;
    p0 = pA - shift_x;    p1 = pA - shift_x + shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    if( l != N - 1 )
    {
        j = i + m - 1;
        k = j + 1;
    }
    else
    {
        j = i + m - 1 - N * m * ( m + 1 ) / 2;
        k = j + 1;
    }
    calculate_triangle( t, x, q, i, j, k, pA, p0, p1, pA0, pA1, p01, residual );
}

void fill_triangles_m_1( double *x, int N, int l, triangle *t, point shift_x, point shift_y, int &q, double &residual )
{
    int i = l + 1;
    int j, k;
    point pA, p0, p1, pA0, pA1, p01;

    point a = shift_x;

    pA = a;
    p0 = pA - shift_x;  p1 = shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    j = 0;
    k = i + 1;
    if( l == N-1 ) k -= N;
    calculate_triangle( t, x, q, i, j, k, pA, p0, p1, pA0, pA1, p01, residual );
}

double approximate( double *x, int m, int N, double h, triangle *t )
{
    int l;
    double alf = 2 * M_PI / N;
    double r = h / (2 * m * ( sin( alf / 2 )));
    point shift_x;
    point shift_y;
    point a;
    int q = 0;

    double residual = 0.0;

    if( m > 2 )
    {
        for( l = 0; l < N; l++ )
        {
            shift_x = polar_to_dec( r, M_PI -  l      * alf );
            shift_y = polar_to_dec( r, M_PI - (l + 1) * alf );

            fill_triangles_III( x, m, N, r, alf, l, t, shift_x, shift_y, q, residual );
            fill_triangles_V  ( x, m, r,    alf, l, t, shift_x, shift_y, q, residual );
            for( int p = m; p >= 4; p-- )
            {
                fill_triangles_VI( x, m, N, p, r, alf, l, t, a, shift_x, shift_y, q, residual );
                fill_triangles_II( x, m, p, l, t, a, shift_x, shift_y, q, residual );
                fill_triangles_IV( x, m, p, l, t, a, shift_x, shift_y, q, residual );
            }
            fill_triangles_VI ( x, m, N, 3, r, alf, l, t, a, shift_x, shift_y, q, residual );
            fill_triangles_IV ( x, m, 3, l, t, a, shift_x, shift_y, q, residual );
            fill_triangles_VII( x, m, N, r, alf, l, t, shift_x, shift_y, q, residual );
        }
    }
    else if( m == 2 )
          {
            for( l = 0; l < N; l++ )
            {
                shift_x = polar_to_dec( r, M_PI -  l      * alf );
                shift_y = polar_to_dec( r, M_PI - (l + 1) * alf );

                fill_triangles_III( x, m, N, r, alf, l, t, shift_x, shift_y, q, residual );

                fill_triangles_V  ( x, m, r,    alf, l, t, shift_x, shift_y, q, residual );
                fill_triangles_VII( x, m, N, r, alf, l, t, shift_x, shift_y, q, residual );
            }
          }
          else if( m == 1 )
                {
                    for( l = 0; l < N; l++ )
                    {
                        shift_x = polar_to_dec( r, M_PI -  l      * alf );
                        shift_y = polar_to_dec( r, M_PI - (l + 1) * alf );

                        fill_triangles_m_1( x, N, l, t, shift_x, shift_y, q, residual );
                    }
                }

    return residual;
}
