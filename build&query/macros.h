//
// Created by  on 2020/7/24.
//

#ifndef P2H_MACROS_H
#define P2H_MACROS_H

#include <x86intrin.h>
#include <cassert>

#define CEIL(X,Y) (((X)-1) / (Y) + 1)

#define POPCNT32(X) (_mm_popcnt_u32(X))

#define E(X) ((X) & ((X)-1))

#define P(X) ((X) ^ -(X))

#endif //P2H_MACROS_H
