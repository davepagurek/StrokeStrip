/*--
    VectorC.h  

    This file is part of the Cornucopia curve sketching library.
    Copyright (C) 2010 Ilya Baran (baran37@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CORNUCOPIA_VECTORC_H_INCLUDED
#define CORNUCOPIA_VECTORC_H_INCLUDED

#include "defs.h"
#include <vector>
#include "Eigen/StdVector"

NAMESPACE_Cornu

enum CircularType
{
    NOT_CIRCULAR,
    CIRCULAR
};

template<typename T> struct default_allocator_traits
{
    typedef std::allocator<T> allocator;
};

template<> struct default_allocator_traits<Eigen::Vector2d>
{
    typedef Eigen::aligned_allocator<Eigen::Vector2d> allocator;
};

//VectorC represents a vector with possibly circular access.
//Warning: inherits class with non-virtual destructor--do not use polymorphically!
template<typename T, typename Alloc = typename default_allocator_traits<T>::allocator >
class VectorC : public std::vector<T, Alloc>
{
public:
    typedef std::vector<T, Alloc> Base;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;

    VectorC() : _circular(NOT_CIRCULAR) {}
    VectorC(int size, CircularType circular) : Base(size), _circular(circular) {}
    VectorC(const Base &base, CircularType circular) : Base(base), _circular(circular) {}
    VectorC(const VectorC &other) : Base(other), _circular(other._circular) {}

    int size() const { return (int)Base::size(); }

    reference operator[](int idx) { return Base::operator[](toLinearIdx(idx)); }
    const_reference operator[](int idx) const { return Base::operator[](toLinearIdx(idx)); }

    //noncircular
    reference flatAt(int idx) { return Base::operator[](idx); }
    const_reference flatAt(int idx) const { return Base::operator[](idx); }

    CircularType circular() const { return _circular; }
    void setCircular(CircularType circular) { _circular = circular; }

    //returns the size for iteration where at each iteration we access elements i, i+1, ..., i+offset
    int endIdx(int offset) const { return _circular ? (int)Base::size() : std::max(0, (int)size() - offset); }

    class Circulator
    {
    public:
        Circulator(const VectorC<T> *ptr, int idx) : _ptr(ptr), _idx(idx), _startIdx(idx) {}

        const_reference operator*() const { return (*_ptr)[_idx]; }
        bool operator==(const Circulator &other) const { return _ptr == other._ptr && _ptr->toLinearIdx(_idx) == _ptr->toLinearIdx(other._idx); }
        bool operator!=(const Circulator &other) const { return !(*this == other); }

        Circulator &operator++() { ++_idx; return *this; }
        Circulator &operator--() { --_idx; return *this; }
        Circulator &operator+=(int x) { _idx += x; return *this; }
        Circulator &operator-=(int x) { _idx -= x; return *this; }
        Circulator operator+(int x) const { return Circulator(_ptr, _idx + x, _startIdx); }
        Circulator operator-(int x) const { return Circulator(_ptr, _idx - x, _startIdx); }

        bool done() const { if(_ptr->circular()) return abs(_idx - _startIdx) >= _ptr->size(); else return _idx < 0 || _idx >= _ptr->size(); }
        int index() const { return _ptr->toLinearIdx(_idx); }

    private:
        Circulator(const VectorC<T> *ptr, int idx, int startIdx) : _ptr(ptr), _idx(idx), _startIdx(startIdx) {}

        const VectorC<T> *_ptr;
        int _idx;
        int _startIdx;
    };

    Circulator beginCirculator() const { return Circulator(this, 0); }
    Circulator circulator(int idx) const { return Circulator(this, toLinearIdx(idx)); }

    int toLinearIdx(int idx) const
    {
        if(!_circular)
            return idx;
        int out = idx % size();
        if(out >= 0)
            return out;
        return out + size();
    }

    int numElems(int from, int to) const //returns the number of elements between from (inclusive) and to (exclusive)
    {
        if(from > to)
            to += size();
        return to - from;
    }

private:
    CircularType _circular;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_VECTORC_H_INCLUDED
