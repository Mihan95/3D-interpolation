#include "funcs.h"
#include <iostream>

using namespace std;

int sum( int m, int p, int t )
{
    int k, s = 0;
    for( k = 0; k <= m - p; k++ )
    {
        s = s + (m - k - t);
    }

    return s;
}

int matrix_msr_size( int m, int N )
{
    if( m == 1)
    {
        return N + 1 + N * 4 +1;
    }
    int n = 0;
    int i1, i2, r, i;
    int l = 0;
    for( r = m; r >= 4; r--)
    {
        i1 = 2 + sum( m, r, 0 ) + l * m * ( m + 1 ) / 2;
        i2 = i1 + r - 4;

        for( i = i1; i <= i2; i++ )
        {
            n += 7;
        }
    }

    i1 = 1 + l * m * ( m + 1 ) / 2;
    i2 = m - 1 + l * m * ( m + 1 ) / 2;
    for( i = i1; i <= i2; i++ )
    {
        n += 7;
    }

    for( r = m; r >= 3; r-- )
    {
        n += 5;
    }

    n += 4;

    for( r = m; r >= 3; r -- )
    {
        n += 7;
    }

    n += 5;

    n *= N;

    n = n + 1 + N + 1;
    return n;
}

void fill_nz_type_I( int n, double *a, int *I, int N, double J, int m )
{//должен быть вынесен за цикл по l
    int j;

    I[0] = n + 1;
    I[1] = I[0] + N;

    a[0] = N * J / 12;
    for( j = 0; j < N; j++ )
    {
        a[ I[0] + j ] = J / 12;
        I[ I[0] + j ] = 1 + j * m * ( m + 1 ) / 2;
    }

    if( m == 1 )
    {
        for( j = 1; j <= N; j++ )
        {
            I[ I[0] + j-1 ] = j;
        }
    }
}

void fill_nz_type_II( double *a, int *I, double J, int l, int m, int r )
{   //r = m; r >= 4; r--
    int i, j, i1, i2;

    i1 = 2 + sum( m, r, 0 ) + l * m * ( m + 1 ) / 2;
    i2 = i1 + r - 4;

    for( i = i1 + 1; i <= i2 + 1; i++ )
    {
        I[i] = I[i - 1] + 6;
    }
    for( i = i1; i <= i2; i++ )
    {
        a[i] = J / 2;
        for( j = 0; j < 6; j++ )
            a[I[i] + j] = J / 12;

        I[ I[i] ]     = i - r;
        I[ I[i] + 1 ] = i - r + 1;
        I[ I[i] + 2 ] = i - 1;
        I[ I[i] + 3 ] = i + 1;
        I[ I[i] + 4 ] = i + r - 2;
        I[ I[i] + 5 ] = i + r - 1;
    }
}

void fill_I_3( int *I, int i1, int i2, int m, int n_it )
{
    int i;
    for( i = i1; i <= i2; i++ )
    {
        I[ I[i] ]     = i - ( n_it + sum( m, n_it + 2, n_it ) + m - n_it + 1 );
        I[ I[i] + 1 ] = i - ( n_it + sum( m, n_it + 2, n_it ) );
        I[ I[i] + 2 ] = i - 1;
        I[ I[i] + 3 ] = i + 1;
        I[ I[i] + 4 ] = i + m - 1;
        I[ I[i] + 5 ] = i + m;
        n_it++;
    }
}

void fill_I_3( int *I, int i1, int i2, int m, int n_it, int q )
{
    int i;
    for( i = i1; i <= i2; i++ )
    {
        I[ I[i] ]     = i - 1;
        I[ I[i] + 1 ] = i + 1;
        I[ I[i] + 2 ] = i + m - 1;
        I[ I[i] + 3 ] = i + m;
        I[ I[i] + 4 ] = i - ( n_it + sum( m, n_it + 2, n_it ) + m - n_it + 1 ) + q;
        I[ I[i] + 5 ] = i - ( n_it + sum( m, n_it + 2, n_it )) + q;
        n_it++;
    }
}

void fill_nz_type_III( double *a, int *I, double J, int l, int m, int N )
{
    int i, j, i3;
    int i1 = 2 + l * m * ( m + 1 ) / 2;
    int i2 = m - 1 + l * m * ( m + 1 ) / 2;

    int q = N * m * ( m + 1 ) / 2;
    int n_it;

    //обработаем соединенный с нулем отдельно
    i3 = 1 + l * m * ( m + 1 ) / 2;
    I[i3 + 1] = I[i3] + 6;
    a[i3] = J / 2;
    for( j = 0; j < 6; j++ )
        a[I[i3] + j] = J / 12;
    if( l != 0 )
    {
        n_it = 1;
        if( l != N-1 )
        {
            I[ I[i3] ]     = 0;
            I[ I[i3] + 1 ] = i3 - (n_it + sum( m, n_it + 2, n_it ) + m - n_it + 1);
            I[ I[i3] + 2 ] = i3 - (n_it + sum( m, n_it + 2, n_it ));
            I[ I[i3] + 3 ] = i3 + 1;
            I[ I[i3] + 4 ] = i3 + m;
            I[ I[i3] + 5 ] = ( l + 1 ) * m * ( m + 1 ) / 2 + n_it;

            if( I[ I[i3] + 5 ] > q ) I[ I[i3] + 5 ] -= q;
        }
        else
        {
            I[ I[i3] ]     = 0;
            I[ I[i3] + 2 ] = i3 - (n_it + sum( m, n_it + 2, n_it ) + m - n_it + 1);
            I[ I[i3] + 3 ] = i3 - (n_it + sum( m, n_it + 2, n_it ));
            I[ I[i3] + 4 ] = i3 + 1;
            I[ I[i3] + 5 ] = i3 + m;
            I[ I[i3] + 1 ] = ( l + 1 ) * m * ( m + 1 ) / 2 + n_it;

            if( I[ I[i3] + 1 ] > q ) I[ I[i3] + 1 ] -= q;
        }
    }
    else
    {
        n_it = 1;
        I[ I[i3] ]     = 0;
        I[ I[i3] + 1 ] = i3 + 1;
        I[ I[i3] + 2 ] = i3 + m;
        I[ I[i3] + 3 ] = m * ( m + 1 ) / 2 + n_it;
        I[ I[i3] + 4 ] = i3 - (n_it + sum( m, n_it + 2, n_it ) + m - n_it + 1) + q;
        I[ I[i3] + 5 ] = i3 - (n_it + sum( m, n_it + 2, n_it )) + q;

        //if( I[ I[i3] + 3 ] > q ) I[ I[i3] + 3 ] -= q;
    }
/**************************************************/
    n_it = 2;
    for( i = i1 + 1; i <= i2 + 1; i++ )
    {
        I[i] = I[i - 1] + 6;
    }
    for( i = i1; i <= i2; i++ )
    {
        a[i] = J / 2;
        for( j = 0; j < 6; j++ )
            a[I[i] + j] = J / 12;
     }

    if( l != 0 )
         fill_I_3( I, i1, i2, m, n_it );
    else fill_I_3( I, i1, i2, m, n_it, q );
}

void fill_nz_type_IV( double *a, int *I, double J, int l, int m, int r )
{   //r = m; r >= 3; r--
    int i;

    i = m + sum( m, r, 1 ) + l * m * ( m + 1 ) / 2;
    I[i + 1] = I[i] + 4;

    a[i] = J / 4;

    a[ I[i] ]     = J / 12;
    a[ I[i] + 1 ] = J / 24;
    a[ I[i] + 2 ] = J / 12;
    a[ I[i] + 3 ] = J / 24;

    I[ I[i]     ] = i - r;
    I[ I[i] + 1 ] = i - r + 1;
    I[ I[i] + 2 ] = i - 1;
    I[ I[i] + 3 ] = i + r - 2;
}

void fill_I_5( int *I, double *a, int i, double J, int m )
{
    (void) J;
    a[ I[i] ]     = J / 24;
    a[ I[i] + 1 ] = J / 12;
    a[ I[i] + 2 ] = J / 24;

    I[ I[i] ]    = i - m;
    I[ I[i] + 1] = i - 1;
    I[ I[i] + 2] = i + m - 1;
}

void fill_I_5( int *I, double *a, int i, double J, int m, int q )
{
    (void) J;
    a[ I[i] ]    = J / 12;
    a[ I[i] + 1] = J / 24;
    a[ I[i] + 2] = J / 24;

    I[ I[i] ]    = i - 1;
    I[ I[i] + 1] = i + m - 1;
    I[ I[i] + 2] = i - m + q;
}

void fill_nz_type_V( double *a, int *I, double J, int l, int m, int N )
{
    int i = m + l * m * ( m + 1 ) / 2;
    int q =     N * m * ( m + 1 ) / 2;

    I[i + 1] = I[i] + 3;

    a[i] = J / 6;

    if( l != 0 )
        fill_I_5( I, a, i, J, m );
    else fill_I_5( I, a, i, J, m, q );
}

void fill_I_6( int *I, int i, int r, int m, int l )
{
    I[ I[i] ]     = i - r;
    I[ I[i] + 1 ] = i - r + 1;
    I[ I[i] + 2 ] = i + 1;
    I[ I[i] + 3 ] = i + r - 1;
    I[ I[i] + 4 ] = m - r + 1 + (l+1) * m * ( m + 1 ) / 2;
    I[ I[i] + 5 ] = m - r + 2 + (l+1) * m * ( m + 1 ) / 2;
}

void fill_I_6( int *I, int i, int r, int m, int l, int q )
{
    I[ I[i] ]     = m - r + 1 + (l+1) * m * ( m + 1 ) / 2 - q;
    I[ I[i] + 1 ] = m - r + 2 + (l+1) * m * ( m + 1 ) / 2 - q;
    I[ I[i] + 2 ] = i - r;
    I[ I[i] + 3 ] = i - r + 1;
    I[ I[i] + 4 ] = i + 1;
    I[ I[i] + 5 ] = i + r - 1;

    /*if( I[ I[i] + 0 ] < 0 ) I[ I[i] + 0 ] += q;
    if( I[ I[i] + 1 ] < 0 ) I[ I[i] + 1 ] += q;*/
}

void fill_nz_type_VI( double *a, int *I, double J, int l, int m, int N, int r, int q )
{   // r = m; r >= 3; r--
    int j;

    int i = 1 + sum( m, r, 0 ) + l * m * ( m + 1 ) / 2;
    I[i + 1] = I[i] + 6;

    a[i] = J / 2;
    for( j = 0; j < 6; j++ )
        a[I[i] + j] = J / 12;

    if( l != N - 1)
        fill_I_6( I, i, r, m, l );
    else fill_I_6( I, i, r, m, l, q );

}

void fill_nz_type_VII( double *a, int *I, double J, int l, int m, int N )
{
    int i;
    int q = N * m * ( m + 1 ) / 2;

    i = ( l + 1 ) * m * ( m + 1 ) / 2;
    I[i + 1] = I[i] + 4;

    a[i] = J / 4;

    a[ I[i] ]     = J / 12;
    a[ I[i] + 1 ] = J / 24;
    a[ I[i] + 2 ] = J / 12;
    a[ I[i] + 3 ] = J / 24;

    if( l != (N - 1) )
    {
        I[ I[i] ]     = i - 2;
        I[ I[i] + 1 ] = i - 1;
        I[ I[i] + 2 ] = i + m - 1;
        I[ I[i] + 3 ] = i + m;
    }
    else
    {
        I[ I[i] ]     = i + m - 1 - q;
        I[ I[i] + 1 ] = i + m - q;
        I[ I[i] + 2 ] = i - 2;
        I[ I[i] + 3 ] = i - 1;
    }
}

int n_matrix_rows( int m, int N )
{
    return m * (m + 1) * N / 2 + 1;
}

void fill_nz_m_1(double *a, int *I, double J, int l, int N )
{
    int i = l+1;
    I[i + 1] = I[i] + 3;
    a[i] = J / 6;

    a[ I[i] ]    = J / 12;
    a[ I[i] + 1] = J / 24;
    a[ I[i] + 2] = J / 24;

    if( l == (N - 1) )
    {
        I[ I[i] ]     = 0;
        I[ I[i] + 1 ] = 1;
        I[ I[i] + 2 ] = 3;
    }
    if( l == 0 )
    {
        I[ I[i] ]     = 0;
        I[ I[i] + 1 ] = 2;
        I[ I[i] + 2 ] = 4;
    }
    if( l != 0 && l != N-1 )
    {
        I[ I[i] ]     = 0;
        I[ I[i] + 1 ] = i - 1;
        I[ I[i] + 2 ] = i + 1;
    }
}

void fill_nz( int m, double *a, int *I, int N, double J )
{
    int n_rows = n_matrix_rows( m, N );
    int l = 0;
    int q = N * m * ( m + 1 ) / 2;

    if( m > 2)
    {
        fill_nz_type_I( n_rows, a, I, N, J, m );
        for( l = 0; l < N; l++ )
        {
            fill_nz_type_III( a, I, J, l, m, N);
            fill_nz_type_V( a, I, J, l, m, N );

            for( int r = m; r >= 4; r-- )
            {
                fill_nz_type_VI( a, I, J, l, m, N, r, q );
                fill_nz_type_II( a, I, J, l, m, r );
                fill_nz_type_IV( a, I, J, l, m, r );
            }
            fill_nz_type_VI( a, I, J, l, m, N, 3, q );
            fill_nz_type_IV( a, I, J, l, m, 3 );
            fill_nz_type_VII( a, I, J, l, m, N );
        }
    } else if( m == 2 )
            {
                fill_nz_type_I( n_rows, a, I, N, J, m );
                for( l = 0; l < N; l++ )
                {
                    fill_nz_type_III( a, I, J, l, m, N);
                    fill_nz_type_V( a, I, J, l, m, N );
                    fill_nz_type_VII( a, I, J, l, m, N );
                }
            }
            else if( m == 1 )
                  {
                        fill_nz_type_I( n_rows, a, I, N, J, m );
                        for( l = 0; l < N; l++ )
                        {
                            fill_nz_m_1( a, I, J, l, N );
                        }
                  }
}
