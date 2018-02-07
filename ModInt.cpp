template<int mod=1000000007>
class ModInt {
  int x;
public:
  ModInt() : x(0){}
  ModInt(int64_t y){ x = y % mod; if(x < 0) x += mod; }
  ModInt &operator += (const ModInt &p){ x += p.x; if(x >= mod) x -= mod; return *this; }
  ModInt &operator -= (const ModInt &p){ x -= p.x; if(x < 0) x += mod; return *this; }
  ModInt &operator *= (const ModInt &p){ x = (int) (1LL * x * p.x % mod); return *this; }
  ModInt &operator /= (const ModInt &p){ *this *= p.inverse(); return *this; }
  ModInt &operator += (const int64_t y){ x = (x + y)%mod; if(x < 0) x += mod; return *this; }
  ModInt &operator -= (const int64_t y){ x = (x - y)%mod; if(x < 0) x += mod; return *this; }
  ModInt &operator *= (const int64_t y){ x = (int) (x * y % mod); return *this; }
  ModInt &operator /= (const int64_t y){ *this *= ModInt(y).inverse(); return *this; }
  ModInt operator -() const { return ModInt(-x); }
  ModInt operator + (const ModInt &p) const { return ModInt(*this) += p; }
  ModInt operator - (const ModInt &p) const { return ModInt(*this) -= p; }
  ModInt operator * (const ModInt &p) const { return ModInt(*this) *= p; }
  ModInt operator / (const ModInt &p) const { return ModInt(*this) /= p; }
  bool operator == (const ModInt &p) const { return x == p.x; }
  bool operator != (const ModInt &p) const { return x != p.x; }
  ModInt operator + (const int64_t y) const { return ModInt(*this) += y; }
  ModInt operator - (const int64_t y) const { return ModInt(*this) -= y; }
  ModInt operator * (const int64_t y) const { return ModInt(*this) *= y; }
  ModInt operator / (const int64_t y) const { return ModInt(*this) /= y; }
  bool operator == (const int64_t y) const { return x == (mod + y%mod)%mod; }
  bool operator != (const int64_t y) const { return x != (mod + y%mod)%mod; }
  ModInt operator = (const int64_t y) { return *this = ModInt(y); }
  ModInt inverse() const {
    int a = x, b = mod, u = 1, v = 0, t;
    while(b > 0){
      t = a/b; a -= t*b; swap(a, b);
      u -= t*v; swap(u, v);
    }
    return ModInt(u);
  }
  ModInt pow(int64_t y) const {
    if(x==0) return ModInt(0);
    int64_t r = 1, t = x;
    while(y > 0){
      if(y&1) r = r*t%mod;
      t = t*t%mod; y >>= 1;
    }
    return ModInt(r);
  }
  friend ostream &operator << (ostream &os, const ModInt<mod> &p) { return os<<p.x; }
  friend istream &operator >> (istream &is, ModInt<mod> &a) { int64_t x; is>>x; a = ModInt<mod>(x); return is; }
};
using Int = ModInt<>;
