//
// Created by on 2020/7/24.
//

#ifndef P2H_BITVECTOR_H
#define P2H_BITVECTOR_H

#include <cstdint>
/**
 * Warning:
 * Pay attention to the bit sequence:
 * Every 8 IDs are placed in an uint8_t unit.
 * Within the uint8_t unit, smaller id is placed
 * at LOWER (less significant) bits.
 */

class BitVector{
public:
    BitVector(int num);
    ~BitVector();

    void SetOnes();
    void SetZeros();

    int CountOnes();

    void And(const BitVector* bv);
    void Or(const BitVector* bv);


    // bit manupliation
    bool GetBit(int pos);
    void SetBit(int pos);
    void UnsetBit(int pos);

    void SetByte(uint8_t byte, int pos);
    uint8_t GetByte(int pos) const;
    uint32_t GetWord(int pos) const;
    int Num() const;
    int NumByte() const;

private:

    uint8_t* data_ = nullptr;
    int num_;
    int num_bytes_;
    int num_words_;
};

inline void BitVector::SetByte(uint8_t byte, int pos) {
    data_[pos] = byte;
}

inline uint8_t BitVector::GetByte(int pos) const {
    return data_[pos];
}

inline int BitVector::Num() const {
    return num_;
}

inline int BitVector::NumByte() const {
    return num_bytes_;
}


#endif //P2H_BITVECTOR_H
