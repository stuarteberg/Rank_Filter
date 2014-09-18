#ifndef RANK_FILTER
#define RANK_FILTER

#include <vigra/multi_array.hxx>
#include <set>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <utility>

#include "sorted_vector.hxx"

namespace vigra
{

template<unsigned int N,
        class T1, class S1,
        class T2, class S2,
        typename std::enable_if<(N == 1), int>::type = 0>
inline void lineRankOrderFilterND(const vigra::MultiArrayView <N, T1, S1> &src,
        vigra::MultiArrayView <N, T2, S2> dest,
        unsigned int half_length, float rank, unsigned int axis = 0)
{
    // Rank must be in the range 0 to 1
    assert((0 <= rank) && (rank <= 1));

    const int rank_pos = round(rank * (2 * half_length));

    // The position of the
    typename vigra::MultiArrayView<N, T1, S1>::difference_type_1 window_begin(0);
    SortedVector<T1> sorted_window;

    // Get the initial sorted window.
    // Include the reflection.
    for (typename vigra::MultiArrayView<N, T1, S1>::difference_type_1 j(half_length); j > 0; j--)
    {
        sorted_window.insert(src[window_begin + j]);
    }

    for (typename vigra::MultiArrayView<N, T1, S1>::difference_type_1 j(0); j <= half_length; j++)
    {
        sorted_window.insert(src[window_begin + j]);
    }

    T1 prev_value;
    T1 next_value;
    while ( window_begin < src.size() )
    {
        dest[window_begin] = sorted_window[rank_pos];
        if ( window_begin < half_length )
        {
            prev_value = src[half_length - window_begin];
        }
        else
        {
            prev_value = src[window_begin - half_length];
        }

        window_begin++;

        if ( window_begin < (src.size() - half_length) )
        {
            next_value = src[window_begin + half_length];
        }
        else
        {
            next_value = src[2 * src.size() - (window_begin + half_length + 2)];
        }

        sorted_window.erase(prev_value);
        sorted_window.insert(next_value);
    }
    dest[window_begin] = sorted_window[rank_pos];
}

template<unsigned int N,
        class T1, class S1,
        class T2, class S2,
        typename std::enable_if<(N > 1), int>::type = 0>
inline void lineRankOrderFilterND(const vigra::MultiArrayView <N, T1, S1> &src,
        vigra::MultiArrayView <N, T2, S2> dest,
        unsigned int half_length, float rank, unsigned int axis = 0)
{
    typename vigra::MultiArrayView<N, T1, S1>::difference_type transposed_axes;

    for (unsigned int i = 0; i < N; i++)
    {
        transposed_axes[i] = i;
    }

    std::swap(transposed_axes[0], transposed_axes[axis]);

    vigra::MultiArray<N, T1> src_transposed(src.transpose(transposed_axes));

    vigra::MultiArrayView<N, T2, S2> dest_transposed(dest.transpose(transposed_axes));


    typename vigra::MultiArrayView<N - 1, T1, S1>::difference_type pos;
    pos = 0;

    bool done = false;
    bool carry = true;
    while (!done)
    {
        lineRankOrderFilterND(src_transposed.bindOuter(pos), dest_transposed.bindOuter(pos), half_length, rank, 0);

        carry = true;
        for (unsigned int i = 0; ( carry && (i < (N - 1)) ); i++)
        {
            if ( (++pos[i]) < src_transposed.shape(i + 1) )
            {
                carry = false;
            }
            else
            {
                pos[i] = 0;
                carry = true;
            }
        }

        done = carry;
    }
}

template<unsigned int N,
        class T1, class S1,
        class T2, class S2>
inline void lineRankOrderFilter(const vigra::MultiArrayView <N, T1, S1> &src,
        vigra::MultiArrayView <N, T2, S2> dest,
        unsigned int half_length, float rank, unsigned int axis = 0)
{
    lineRankOrderFilterND(src, dest, half_length, rank, axis);
}

}
#endif
