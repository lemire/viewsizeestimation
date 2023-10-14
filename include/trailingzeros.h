#ifndef NGRAMESTIMATOR_H
#define NGRAMESTIMATOR_H

// Cette fonction renvoi le nombre de bits les moins significatifs mis à zéro du
// nombre k
inline int trailingzeros(uint64 k, int L = 19) {
  int ans = 0;
  while (((k & 1) == 0) && (ans < L - 1)) {
    ans += 1;
    k >>= 1;
  }
  return ans;
}

// Cette fonction renvoi le nombre de bits les plus significatifs mis à zéro du
// nombre k
inline int lefttrailingzeros(uint64 k, int L = 19) {
  int ans = 0;
  uint64 mask = 1;
  mask <<= L - 1; // << (int) pow(2, (double) (L-1));
  while (((k & mask) == 0) && (ans < L - 1)) {
    ans += 1;
    k <<= 1;
  }
  return ans;
}

#endif
