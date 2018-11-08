/*
 * Copyright (c) 2013 - Idiap Research Institute - http://www.idiap.ch/
 * 
 * See file COPYING for the licence associated to this software
 * 
 * This file is a part of the TopicModels4Video software version: academic-1.0
*/

#ifndef ENSURE_PROPER_TIME_UNIT
#define ENSURE_PROPER_TIME_UNIT(t) // Suppose the time units are ok
#endif

template <class T> class SmartPointer {
public:
    SmartPointer(T* pT) : pT(pT) {}
    ~SmartPointer() {
        delete pT;
    }
    T* operator*() {return pT;}
    T* operator->() {return pT;}
private:
    T* pT;
};
template <class T> class SmartArray {
public:
    SmartArray(T* p) : p(p) {}
    SmartArray(size_t size) : p(new T[size]) {}
    ~SmartArray() {
        delete[] p;
    }
    T& operator[](int index) {return p[index];}
    T* operator*() {return p;}
    T* operator->() {return p;}
private:
    T* p;
};
 
