#ifndef _FUNOPT_MACROS_H_
#define _FUNOPT_MACROS_H_
#include <cstdio>
#include <process.h>

inline void m_assert(bool b, char* msg, char* file, int line) {
    if(!b) {
        printf("[assert] %s in %s line %d\n", msg, file, line);
        abort();
    }
}

#define massert(b, msg) m_assert(b, msg, __FILE__, __LINE__) 

#endif
