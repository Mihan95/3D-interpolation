#include <math.h>
#include "funcs.h"
#include <iostream>

using namespace std;

double f( double x, double y )
{
    return x+y;
}

point polar_to_dec( double r, double fi )
{
    return point(r * sin(fi), r * cos(fi));
}

double add_to_b( double J, point pA, point p0, point p1, point pA0, point pA1, point p01 )
{
    return J * ( f( pA.x, pA.y ) / 32 + ( f( p0.x, p0.y ) + f( p1.x, p1.y )) / 192
                 +(f( pA0.x, pA0.y ) + f( pA1.x, pA1.y ) ) * 5 / 96 + f( p01.x, p01.y ) / 48 );
}

void calculate_b( double *b, double J, point pA, point p0, point p1, point &pA0, point &pA1, point &p01, int i )
{
    pA0  = pA1;    pA1 = ( pA + p1 ) / 2;    p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b(J, pA, p0, p1, pA0, pA1, p01);
}

void fill_vector_I( double *b, int N, double J, double r, double alf )
{
    point a(0, 0);
    point pA, p0, p1, pA0, pA1, p01;

    pA = a;

    for( int l = 0; l < N; l++ )
    {
        p0  = polar_to_dec( r, M_PI -  l      * alf );
        p1  = polar_to_dec( r, M_PI - (l + 1) * alf );

        pA0 = ( pA + p0 ) / 2;
        pA1 = ( pA + p1 ) / 2;
        p01 = ( p0 + p1 ) / 2;

        b[0] = b[0] + add_to_b(J, pA, p0, p1, pA0, pA1, p01);
    }
}

void fill_vector_III( double *b, int m, double J, double r, double alf, int l, point shift_x, point shift_y, point shift_z )
{
    int i, j = 2;
    int i1 = 1 + l * m * ( m + 1 ) / 2;
    int i2 = m - 1 + l * m * ( m + 1 ) / 2;

    point pA, p0, p1, pA0, pA1, p01;
    point c(0, 0);
    point a = polar_to_dec( r,     M_PI - l * alf );
    point d = polar_to_dec( j * r, M_PI - l * alf );

    for( i = i1; i <= i2; i++ )
    {
        pA  = a;

        p0   = c;        p1  = c + shift_y;
        pA0  = ( pA + p0 ) / 2;    pA1 = ( pA + p1 ) / 2;    p01 = ( p0 + p1 ) / 2;
        b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

        p0   = p1;       p1 = p0 + shift_x;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0   = p1;       p1 = d;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0   = p1;       p1 = a + shift_z;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0   = p1;       p1 = c + shift_z;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0   = p1;       p1 = c;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        j++;
        c = a; a = d; d = polar_to_dec( j * r, M_PI - l * alf );
    }
}

void fill_vector_V( double *b, int m, double J, double r, double alf, int l, point shift_x, point shift_y, point shift_z )
{
    int i = m + l * m * ( m + 1 ) / 2;
    point pA, p0, p1, pA0, pA1, p01;
    point c = polar_to_dec( r * (m - 1), M_PI - l * alf );
    point a = polar_to_dec( r *  m,      M_PI - l * alf );

    (void) shift_x;

    pA = a;
    p0 = c;       p1 = c + shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

    p1 = p0;      p0 = c + shift_z;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );
}

void fill_vector_VI( double *b, int m, double J, int l, int p, double r, double alf, point &a, point shift_x, point shift_y,
                     point shift_z )
{
    // p = m; p >= 3; p--
    int i = 1 + sum( m, p, 0 ) + l * m * ( m + 1 ) / 2;

    point pA, p0, p1, pA0, pA1, p01;
    (void) shift_z;

    a = polar_to_dec( r, M_PI - l * alf );
    a = a + shift_y * (m - p + 1);

    pA = a;
    p0 = pA - shift_x;    p1 = p0 + shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

    p0 = p1;    p1 = p1 + shift_x;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    p0 = p1;    p1 = pA + shift_x;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    p0 = p1;    p1 = p1 - shift_y;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    p0 = p1;    p1 = p1 - shift_x;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    p0 = p1;    p1 = pA - shift_x;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    a = a + shift_x;
}

void fill_vector_II( double *b, int m, double J, int l, int p, point &a, point shift_x, point shift_y, point shift_z )
{   //p = m; p >= 4; p--

    int i, i1, i2;
    point pA, p0, p1, pA0, pA1, p01;
    (void) shift_z;

    i1 = 2 + sum( m, p, 0 ) + l * m * ( m + 1 ) / 2;
    i2 = i1 + p - 4;


    for( i = i1; i <= i2; i++ )
    {
        pA = a;
        p0 = pA - shift_x;    p1 = p0 + shift_y;
        pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
        b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

        p0 = p1;    p1 = p1 + shift_x;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0 = p1;    p1 = pA + shift_x;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0 = p1;    p1 = p1 - shift_y;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0 = p1;    p1 = p1 - shift_x;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        p0 = p1;    p1 = pA - shift_x;
        calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

        a = a + shift_x;
    }
}

void fill_vector_IV( double *b, int m, double J, int l, int p, point &a, point shift_x, point shift_y, point shift_z )
{
    //p = m; p >= 3; p--
    int i = m + sum( m, p, 1 ) + l * m * ( m + 1 ) / 2;
    point pA, p0, p1, pA0, pA1, p01;

    (void) shift_z;

    pA = a;
    p0 = pA - shift_y; p0 = p0 + shift_x;    p1 = pA - shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

    p0 = p1;    p1 = pA - shift_x;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    p0 = p1;    p1 = p1 + shift_y;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );
}

void fill_vector_VII( double *b, int m, double J, int l, double r, double alf, point shift_x, point shift_y, point shift_z )
{
    int i = ( l + 1 ) * m * ( m + 1 ) / 2;
    point pA, p0, p1, pA0, pA1, p01;
    (void) shift_z;

    point a = polar_to_dec( r, M_PI - l * alf );
    a = a + shift_y * (m - 1);

    pA = a;
    p0 = pA - shift_y; p0 = p0 + shift_x;    p1 = pA - shift_y;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

    p0 = p1;    p1 = pA - shift_x;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );

    p0 = p1;    p1 = p1 + shift_y;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );
}

void fill_vector_m_1( double *b, double J, int l, point shift_x, point shift_y, point shift_z )
{
    int i = 1 + l;
    point pA, p0, p1, pA0, pA1, p01;

    point a = shift_x;

    p0 = shift_z; p1 = a - shift_x;
    pA0 = ( pA + p0 ) / 2;  pA1 = ( pA + p1 ) / 2;  p01 = ( p0 + p1 ) / 2;
    b[i] = b[i] + add_to_b( J, pA, p0, p1, pA0, pA1, p01 );

    p0 = p1; p1 = shift_y;
    calculate_b( b, J, pA, p0, p1, pA0, pA1, p01, i );
}

void fill_vector( double *b, double J, int m, int N, double h/*длина стороны*/ )
{
    int l;
    double alf = 2 * M_PI / N;
    double r = h / (2 * m * ( sin( alf / 2 )));
    point shift_x;
    point shift_y;
    point shift_z;
    point a;

    if( m > 2 )
    {
        fill_vector_I( b, N, J, r, alf );
        for( l = 0; l < N; l++ )
        {
            shift_x = polar_to_dec( r, M_PI -  l      * alf );
            shift_y = polar_to_dec( r, M_PI - (l + 1) * alf );
            shift_z = polar_to_dec( r, M_PI - (l - 1) * alf );

            fill_vector_III( b, m, J, r, alf, l, shift_x, shift_y, shift_z );
            fill_vector_V  ( b, m, J, r, alf, l, shift_x, shift_y, shift_z );
            for( int p = m; p >= 4; p-- )
            {
                fill_vector_VI( b, m, J, l, p, r, alf, a, shift_x, shift_y, shift_z );
                fill_vector_II( b, m, J, l, p, a, shift_x, shift_y, shift_z );
                fill_vector_IV( b, m, J, l, p, a, shift_x, shift_y, shift_z );
            }
            fill_vector_VI ( b, m, J, l, 3, r, alf, a, shift_x, shift_y, shift_z );
            fill_vector_IV ( b, m, J, l, 3, a,   shift_x, shift_y, shift_z );
            fill_vector_VII( b, m, J, l, r, alf, shift_x, shift_y, shift_z );
        }
    }
    else if(  m == 2 )
            {
                fill_vector_I( b, N, J, r, alf );
                for( l = 0; l < N; l++ )
                {
                    shift_x = polar_to_dec( r, M_PI -  l      * alf );
                    shift_y = polar_to_dec( r, M_PI - (l + 1) * alf );
                    shift_z = polar_to_dec( r, M_PI - (l - 1) * alf );

                    fill_vector_III( b, m, J, r, alf, l, shift_x, shift_y, shift_z );
                    fill_vector_V  ( b, m, J, r, alf, l, shift_x, shift_y, shift_z );
                    fill_vector_VII( b, m, J, l, r, alf, shift_x, shift_y, shift_z );
                }
            }
            else if( m == 1 )
                 {
                    fill_vector_I( b, N, J, r, alf );

                    for( l = 0; l < N; l++ )
                    {
                        shift_x = polar_to_dec( r, M_PI -  l      * alf );
                        shift_y = polar_to_dec( r, M_PI - (l + 1) * alf );
                        shift_z = polar_to_dec( r, M_PI - (l - 1) * alf );
                        fill_vector_m_1( b, J, l, shift_x, shift_y, shift_z );
                    }

                 }

}
