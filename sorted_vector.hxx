#ifndef SORTED_VECTOR
#define SORTED_VECTOR

#include <deque>
#include <algorithm>

template <typename T>
class SortedVector
{
public:
    std::deque<T> m_vec ;

public:
    SortedVector()
    {
    }
    
    ~SortedVector()
    {
    }

    auto begin() -> decltype( m_vec.begin() )
    {
        return m_vec.begin();
    }
    
    auto end() -> decltype( m_vec.end() )
    {
        return m_vec.end();
    }
    
    void insert(T const & t)
    {
        auto it = std::upper_bound(m_vec.begin(), m_vec.end(), t);
        m_vec.insert(it, t);
    }
    
    void erase(T const & t)
    {
        auto it = std::lower_bound(m_vec.begin(), m_vec.end(), t);
        assert(*it == t && "Can't erase non-existent element!");
        m_vec.erase(it);
    }

    T operator[] (unsigned long i) const 
    {
        return m_vec[i];
    }
    
    auto size() -> decltype( m_vec.size() )
    {
        return m_vec.size();
    }
};

#endif
