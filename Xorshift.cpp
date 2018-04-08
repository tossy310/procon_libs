template<typename T=uint32_t>
class Xorshift128 {
private:
  static constexpr T MASK = 0x7FFFFFFF;
  T x = 123456789, y = 987654321, z = 1000000007, w;
public:
  T rnd(){
    T t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w & MASK;
  }
  Xorshift128(const T seed = 1000000009) : w(seed) {}
};
