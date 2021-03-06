// TODO template使ってLENの定義周りを綺麗にしたい

#define LEN 200000
class RollingHash {
private:
  static const int64_t B1 = 10007, B2 = 10009;
  static const int64_t MD1 = 1000000007, MD2 = 1000000009;
  static int64_t pow1[], pow2[];
  static void initHash(){
    pow1[0]=1;
    rep(i,LEN+4) pow1[i+1] = pow1[i]*B1 % MD1;
    pow2[0]=1;
    rep(i,LEN+4) pow2[i+1] = pow2[i]*B2 % MD2;
  }
public:
  int64_t h1[LEN+5], h2[LEN+5];
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
  pair<int,int> get(int from, int l) const {
    if(l == 0) return {0,0};
    if(from == 0) return {h1[l-1], h2[l-1]};
    int ret1 = (h1[from+l-1] - h1[from-1]*pow1[l])%MD1;
    if(ret1<0) ret1 += MD1;
    int ret2 = (h2[from+l-1] - h2[from-1]*pow2[l])%MD2;
    if(ret2<0) ret2 += MD2;
    return {ret1, ret2};
  }
};
int64_t RollingHash::pow1[LEN+5] = {0}, RollingHash::pow2[LEN+5] = {0};

// Validation ARC055 C ABCAC
