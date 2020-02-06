#ifndef RUNIONS_GLOBAL_H
#define RUNIONS_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(RUNIONS_LIBRARY)
#  define RUNIONSSHARED_EXPORT Q_DECL_EXPORT
#else
#  define RUNIONSSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // RUNIONS_GLOBAL_H
