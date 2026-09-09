#include "comparaison_masse_native.cpp"
static_assert(sizeof(Pivot)==96 && sizeof(Summary)==64,
              "ABI des preuves incompatible avec le pont Rust");
