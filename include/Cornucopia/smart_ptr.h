/*--
    smart_ptr.h  

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

#ifndef CORNUCOPIA_SMART_PTR_H_INCLUDED
#define CORNUCOPIA_SMART_PTR_H_INCLUDED

#include "defs.h"
#include <algorithm>

NAMESPACE_Cornu

template<typename T> class smart_ptr;

class smart_base
{
private:
    mutable int _refCount;

public:
    smart_base() : _refCount(0) {}
    smart_base(const smart_base &) : _refCount(0) {}
    virtual ~smart_base() {}

    //assigning to a smart_base should not change the reference count
    smart_base &operator=(const smart_base &) { return *this; }

protected:
    template<class U> friend class smart_ptr;

    virtual void addRef() const
    {
        ++_refCount;
    }
    virtual void releaseRef() const
    {
        bool free = false;
        free = (--_refCount <= 0);
        if (free)
            const_cast<smart_base *>(this)->freeRef();
    }
    virtual void freeRef() { delete this; }

    int getRefCount() const { return _refCount; }
};

#define CORNU_SMART_TYPEDEFS(classname) \
    typedef Cornu::smart_ptr<classname> classname##Ptr; \
    typedef Cornu::smart_ptr<const classname> classname##ConstPtr;

#define CORNU_SMART_FORW_DECL(classname) \
    class classname; \
    CORNU_SMART_TYPEDEFS(classname)

template<typename T>
class smart_ptr
{
public:
    typedef T element_type;

    smart_ptr() : ptr(0), typedPtr(0) {}

    smart_ptr(T *inPtr) : ptr(inPtr), typedPtr(inPtr)
    {
        if(ptr)
            ptr->addRef();
    }

    template<class U>
    smart_ptr(const smart_ptr<U> &smartPtr) : ptr(smartPtr.ptr), typedPtr(smartPtr.typedPtr)
    {
        if(ptr)
            ptr->addRef(); 
    }

    smart_ptr(const smart_ptr &smartPtr) : ptr(smartPtr.ptr), typedPtr(smartPtr.typedPtr)
    {
        if(ptr)
            ptr->addRef(); 
    }

    ~smart_ptr()
    {
        if(ptr)
            ptr->releaseRef();
    }


    T *operator->() const { return get(); }
    T &operator*() const { return *get(); }

    template<class U>
    smart_ptr<T> &operator=(const smart_ptr<U> &other)
    {
        smart_ptr<T>(other).swap(*this);
        return *this;
    }

    smart_ptr<T> &operator=(const smart_ptr &other)
    {
        smart_ptr<T>(other).swap(*this);
        return *this;
    }

    smart_ptr<T> &operator=(T *other)
    {
        smart_ptr<T>(other).swap(*this);
        return *this;
    }

    template<class U> bool operator<(const smart_ptr<U> &other) const { return ptr < other.ptr; }

    template<class U> bool operator==(const smart_ptr<U> &other) const { return ptr == other.ptr; }
    bool operator==(T *other) const { return get() == other; }
    template<class U> bool operator!=(const smart_ptr<U> &other) const { return ptr != other.ptr; }
    bool operator!=(T *other) const { return get() != other; }

    typedef const smart_base *smart_ptr<T>::*unspecified_bool_type;
    operator unspecified_bool_type () const
    {
        return ptr == 0 ? 0 : &smart_ptr<T>::ptr;
    }

    T *get() const { return typedPtr; }
    void reset() 
    {
        if(ptr)
            ptr->releaseRef();
        ptr = 0;
        typedPtr = 0;
    }
    void swap(smart_ptr<T> &other) 
    { 
        std::swap(ptr, other.ptr);
        std::swap(typedPtr, other.typedPtr); 
    }

    template<class U> friend class smart_ptr; 

private:
    const smart_base *ptr;
    T *typedPtr;
};

template<class T>
bool operator==(T *ptr1, const smart_ptr<T> &ptr2) { return ptr1 == ptr2.get(); }
template<class T>
bool operator!=(T *ptr1, const smart_ptr<T> &ptr2) { return ptr1 != ptr2.get(); }

template<class T, class U>
smart_ptr<T> dynamic_pointer_cast(const smart_ptr<U> &p) {
    return dynamic_cast<T*>(p.get());
}
template<class T, class U>
smart_ptr<T> static_pointer_cast(const smart_ptr<U> &p) {
    return static_cast<T*>(p.get());
}
template<class T, class U>
smart_ptr<T> const_pointer_cast(const smart_ptr<U> &p) {
    return const_cast<T*>(p.get());
}

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_SMART_PTR_H_INCLUDED
