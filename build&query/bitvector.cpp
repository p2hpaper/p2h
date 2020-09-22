//
//


#include "bitvector.h"
#include "macros.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
//#include <x86intrin.h>
#include <nmmintrin.h>

BitVector::BitVector(int num):
num_(num), num_bytes_(CEIL(num, 8)), num_words_(CEIL(num, 32))
{
    int bytes_needed = CEIL(num, 32) * 4;
    data_ = (uint8_t*) malloc(sizeof(uint8_t) * bytes_needed);
    SetZeros(); // The bits should be unset initially.
}

BitVector::~BitVector() {
    free(data_);
}

bool BitVector::GetBit(int pos) {
    int byte_id = pos / 8;
    int offset = pos % 8;
    uint8_t mask = 1 << offset;
    return (data_[byte_id] & mask);
}

void BitVector::SetBit(int pos) {
    int byte_id = pos / 8;
    int offset = pos % 8;
    uint8_t mask = 1 << offset;
    data_[byte_id] |= mask;
}

void BitVector::UnsetBit(int pos) {
    int byte_id = pos / 8;
    int offset = pos % 8;
    uint8_t mask = 1 << offset;
    data_[byte_id] &= ~mask;
}

void BitVector::SetOnes() {
    memset(data_, 0xff, CEIL(num_, 32) * 4);
}


void BitVector::SetZeros() {
    memset(data_, 0x0, CEIL(num_, 32) * 4);
}

/*
inline int BitVector::CountOnes() {
    int count = 0;
    const uint32_t * words = (uint32_t*)data_;
    //int num_words = CEIL(num_, 32);

    for(int i = 0; i < num_words_; i++){
        count += _mm_popcnt_u32(words[i]);
    }
    return count;
}*/

void BitVector::And(const BitVector* bv){

    uint32_t* words = (uint32_t*) data_;
    for(int i = 0; i < num_words_; i++){
        words[i] &= bv->GetWord(i);
    }
}

void BitVector::Or(const BitVector* bv){
    uint32_t* words = (uint32_t*) data_;
    for(int i = 0; i < num_words_; i++){
        words[i] |= bv->GetWord(i);
    }
}

uint32_t BitVector::GetWord(int pos) const {
    uint32_t* words = (uint32_t*) data_;
    return words[pos];
}





