#include <bits/stdc++.h>
using namespace std;
uint64_t bin_pow(
    uint64_t n, uint64_t p,
    uint64_t mod) { /**    n*m = 1 (mod p)  =>  m = n**(p-2) (mod p)    **/
  uint64_t res = 1;
  while (p) {
    if (p & 1)
      res = ((__uint128_t)res * n) % mod;
    n = ((__uint128_t)n * n) % mod;
    p >>= 1;
  }
  return res;
}

struct montgomery {
  uint64_t n, nr;

  constexpr montgomery(uint64_t n) : n(n), nr(1) {
    // log(2^64) = 6
    for (int i = 0; i < 6; i++)
      nr *= 2 - n * nr;
  }

  [[nodiscard]] uint64_t reduce(__uint128_t x) const {
    uint64_t q = __uint128_t(x) * nr;
    uint64_t m = ((__uint128_t)q * n) >> 64;
    uint64_t res = (x >> 64) + n - m;
    if (res >= n)
      res -= n;
    return res;
  }

  [[nodiscard]] uint64_t multiply(uint64_t x, uint64_t y) const {
    return reduce((__uint128_t)x * y);
  }

  [[nodiscard]] uint64_t transform(uint64_t x) const {
    return (__uint128_t(x) << 64) % n;
  }
};

vector<int> bit_sort(int n) {
  int h = -1;
  vector<int> rev(n, 0);
  int skip = __lg(n) - 1;
  for (int i = 1; i < n; ++i) {
    if (!(i & (i - 1)))
      ++h;
    rev[i] = rev[i ^ (1 << h)] | (1 << (skip - h));
  }
  return rev;
}

const uint64_t mod = 2524775926340780033, gen = 3;
//const uint64_t mod = 998244353, gen = 3;
void ntt(vector<uint64_t>& a, vector<int>& rev, montgomery& red, uint64_t inv_n,
         uint64_t root, uint64_t inv_root, bool invert) {
  int n = (int)a.size();

  for (int i = 0; i < n; ++i)
    if (i < rev[i])
      swap(a[i], a[rev[i]]);

  uint64_t w = invert ? inv_root : root;
  vector<uint64_t> W(n >> 1, red.transform(1));
  for (int i = 1; i < (n >> 1); ++i)
    W[i] = red.multiply(W[i - 1], w);

  int lim = __lg(n);
  for (int i = 0; i < lim; ++i)
    for (int j = 0; j < n; ++j)
      if (!(j & (1 << i))) {
        uint64_t t = red.multiply(a[j ^ (1 << i)],
                                  W[(j & ((1 << i) - 1)) * (n >> (i + 1))]);
        a[j ^ (1 << i)] = a[j] >= t ? a[j] - t : a[j] + mod - t;
        a[j] = a[j] + t < mod ? a[j] + t : a[j] + t - mod;
      }

  if (invert)
    for (int i = 0; i < n; i++)
      a[i] = red.multiply(a[i], inv_n);
}

void mul(vector<uint64_t>& a, vector<uint64_t>& b) {
  montgomery red(mod);
  for (auto& x : a)
    x = red.transform(x);
  for (auto& x : b)
    x = red.transform(x);

  int n = 1;
  while (n < a.size() || n < b.size())
    n <<= 1;
  n <<= 1;
  a.resize(n);
  b.resize(n);

  uint64_t inv_n = red.transform(bin_pow(n, mod - 2, mod));
  uint64_t root = red.transform(bin_pow(gen, (mod - 1) / n, mod));
  uint64_t inv_root = red.transform(bin_pow(red.reduce(root), mod - 2, mod));
  auto rev = bit_sort(n);

  ntt(a, rev, red, inv_n, root, inv_root, false);
  ntt(b, rev, red, inv_n, root, inv_root, false);

  for (int i = 0; i < n; i++)
    a[i] = red.multiply(a[i], b[i]);
  ntt(a, rev, red, inv_n, root, inv_root, true);

  for (auto& x : a)
    x = red.reduce(x);
}

int main() {

  montgomery red(mod);
  std::cout << red.nr;
}