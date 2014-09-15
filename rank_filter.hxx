#include <list>
#include <vigra/multi_array.hxx>
#include <vigra/linear_algebra.hxx>
#include <set>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iterator>

#define DISPLAY(X)  std::cout << #X << " = " << X << std::endl

template <class T, class C>
std::ostream& operator<<(std::ostream& out, const vigra::MultiArrayView<1, T, C>& that)
{
    std::ios::fmtflags flags = out.setf(std::ios::right |std::ios::fixed,std::ios::adjustfield |std::ios::floatfield);
    for(vigra::MultiArrayIndex i = 0; i < that.size(); ++i)
    {
        out << that(i) << " ";
    }
    out << std::endl;
    out.setf(flags);
    return out;
}

template <class T1, class S1,
          class T2, class S2>
inline void lineRankOrderFilter(const vigra::MultiArrayView<1, T1, S1> & src,
                                vigra::MultiArrayView<1, T2, S2> dest,
                                unsigned long half_length, double rank)
{
    // Will ignore boundaries initially.
    // Then will try adding reflection.

    // Rank must be in the range 0 to 1
    assert((0 <= rank) && (rank <= 1));

    const int rank_pos = round(rank * (2*half_length));

    // The position of the
    typename vigra::MultiArrayView<1, T1, S1>::difference_type_1 window_begin(0);
    std::multiset<T1> sorted_window;
    std::list< typename std::multiset<T1>::iterator > window_iters;

    // Get the initial sorted window.
    // Include the reflection.
    for (typename vigra::MultiArrayView<1, T1, S1>::difference_type_1 j(half_length); j > 0; j--)
    {
        window_iters.push_back(sorted_window.insert(src[window_begin + j]));
    }
    for (typename vigra::MultiArrayView<1, T1, S1>::difference_type_1 j(0); j <= half_length; j++)
    {
        window_iters.push_back(sorted_window.insert(src[window_begin + j]));
    }

    typename std::multiset<T1>::iterator rank_point = sorted_window.begin();

    for (int i = 0; i < rank_pos; i++)
    {
        rank_point++;
    }

    while ( window_begin < src.size() )
    {
        std::cout << "sorted_window = ";
        for (auto sorted_window_iter = sorted_window.begin(); sorted_window_iter != sorted_window.end(); sorted_window_iter++)
        {
            std::cout << *sorted_window_iter << " " ;
        }
        std::cout << std::endl;

        dest[window_begin] = *rank_point;

        typename std::multiset<T1>::iterator prev_iter(window_iters.front());
        T1 prev_value = *prev_iter;
        window_iters.pop_front();

        window_begin++;

        T1 next_value;
        if ( window_begin < (src.size() - half_length) )
        {
            next_value = src[window_begin + half_length];
        }
        else
        {
            next_value = src[2 * src.size() - (window_begin + half_length + 2)];
        }

        if ( ( *rank_point < *prev_iter ) && ( *rank_point <= next_value ) )
        {
            sorted_window.erase(prev_iter);
            window_iters.push_back(sorted_window.insert(next_value));
        }
        else if ( ( *rank_point >= *prev_iter ) && ( *rank_point > next_value ) )
        {
            if ( rank_point == prev_iter )
            {
                window_iters.push_back(sorted_window.insert(next_value));
                rank_point--;

                sorted_window.erase(prev_iter);
            }
            else
            {
                sorted_window.erase(prev_iter);
                window_iters.push_back(sorted_window.insert(next_value));
            }
        }
        else if ( ( *rank_point < *prev_iter ) && ( *rank_point > next_value ) )
        {
            sorted_window.erase(prev_iter);
            window_iters.push_back(sorted_window.insert(next_value));

            rank_point--;
        }
        else if ( ( *rank_point >= *prev_iter ) && ( *rank_point <= next_value ) )
        {
            if (rank_point == prev_iter)
            {
                window_iters.push_back(sorted_window.insert(next_value));
                rank_point++;

                sorted_window.erase(prev_iter);
            }
            else
            {
                sorted_window.erase(prev_iter);
                window_iters.push_back(sorted_window.insert(next_value));

                rank_point++;
            }
        }
    }
}
