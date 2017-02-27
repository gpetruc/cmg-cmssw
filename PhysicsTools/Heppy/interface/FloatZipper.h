#ifndef PhysicsTools_Heppy_FloatZipper_h
#define PhysicsTools_Heppy_FloatZipper_h

#include <cstdint>

namespace heppy {
    class FloatZipper {
        public:
            FloatZipper() : shift(0),mask(0xFFFFFFFF),test(0),maxn(0) {}
            FloatZipper(int bits) :
                shift((23-bits)),    // bits I throw away
                mask((0xFFFFFFFF >> (shift)) << (shift)), // mask for truncation
                test(1 << (shift-1)), // most significant bit I throw away
                maxn((1<<bits)-2) {}  // max number I can increase before overflowing

            void zip(unsigned int n, float *arr) {
                for (unsigned int i = 0; i < n; ++i) {
                    reduceMantissaToNbitsRounding(arr[i]);
                }
            }

        private:
            const int      shift; // bits I throw away
            const uint32_t mask;  // mask for truncation
            const uint32_t test;  // most significant bit I throw away
            const uint32_t maxn;  // max number I can increase before overflowing

            inline void reduceMantissaToNbitsRounding(float &f) {
                constexpr uint32_t low23 = (0x007FFFFF); // mask to keep lowest 23 bits = mantissa
                constexpr uint32_t  hi9  = (0xFF800000); // mask to keep highest 9 bits = the rest
                union { float flt; uint32_t i32; } conv;
                conv.flt=f;
                if (conv.i32 & test) { // need to round
                    uint32_t mantissa = (conv.i32 & low23) >> shift;
                    if (mantissa < maxn) mantissa++;
                    conv.i32 = (conv.i32 & hi9) | (mantissa << shift);
                } else {
                    conv.i32 &= mask;
                }
                f = conv.flt;
            }
    };
}

#endif
