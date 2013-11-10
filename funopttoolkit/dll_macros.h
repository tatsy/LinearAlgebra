#ifndef _DLL_MACROS_H_
#define _DLL_MACROS_H_

#ifdef __EXPORT__
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif

#endif
