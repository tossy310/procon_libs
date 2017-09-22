// TODO template使ってLENの定義周りを綺麗にしたい

#define LEN 200000
class RollingHash {
private:
  static const long B1 = 10007, B2 = 10009;
  static const long MD1 = 1000000007, MD2 = 1000000009;
  static long pow1[], pow2[];
public:
  static void initHash(){
    pow1[0]=1;
    rep(i,LEN+4) pow1[i+1] = pow1[i]*B1 % MD1;
    pow2[0]=1;
    rep(i,LEN+4) pow2[i+1] = pow2[i]*B2 % MD2;
  }
  long h1[LEN+5], h2[LEN+5];
  int n;
  string &s;
  RollingHash(string &str) : s(str){
    if(pow1[0] == 0 || pow2[0]==0) initHash();
    n = s.size();
    assert(n<=LEN);
    h1[0] = h2[0] = s[0];
    rep(i,1,n){
      h1[i] = (h1[i-1]*B1 + s[i])%MD1;
      h2[i] = (h2[i-1]*B2 + s[i])%MD2;
    }
  }
  pair<long,long> get(int from, int l){
    if(l == 0) return mp(0L,0L);
    if(from == 0) return mp(h1[l-1], h2[l-1]);
    long ret1 = (h1[from+l-1] - h1[from-1]*pow1[l])%MD1;
    if(ret1<0) ret1 += MD1;
    long ret2 = (h2[from+l-1] - h2[from-1]*pow2[l])%MD2;
    if(ret2<0) ret2 += MD2;
    return mp(ret1, ret2);
  }
};
long RollingHash::pow1[LEN+5], RollingHash::pow2[LEN+5];
