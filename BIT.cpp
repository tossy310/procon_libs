// Binary Indexed Tree
template<typename T>
class BIT {
private:
  vector<T> bit;
  int n;
public:
  BIT(int _n) : n(_n) {
    bit = vector<T>(n+1, 0);
  }
  void add(int v, T a){ //0-indexed
    for(int x=v+1; x<=n; x += x&(-x)) bit[x] += a;
  }
  T sum(int v){ //0-indexed
    T ret=0;
    for(int x=v+1; x>0; x -= x&(-x)) ret += bit[x];
    return ret;
  }
  int lower_bound(T w){ //w以上となる最小のsumの位置(0-indexed)
    if(w<=0) return 0;
    int x=0, d=0;
    while(n > (1<<d)) d++;
    for(int k=(1<<(d-1)); k>0; k/=2){
      if(x+k<=n && bit[x+k]<w){
        w -= bit[x+k];
        x += k;
      }
    }
    return x;
  }
}; // END class BIT


// 2-dimention Binary Indexed Tree
template<typename T>
class BIT2D {
private:
  vector<vector<T> > bit;
  int x,y;
public:
  BIT2D(int _x, int _y) : x(_x), y(_y) {
    bit = vector<vector<T> >(x+1, vector<T>(y+1, 0));
  }
  void add(int _a, int _b, T w){ //0-indexed
    for(int a=_a+1; a<=x; a += a&(-a)) {
      for(int b=_b+1; b<=y; b += b&(-b)) {
        bit[a][b] += w;
      }
    }
  }
  T sum(int _a, int _b){ //0-indexed
    T ret=0;
    for(int a=_a+1; a>0; a -= a&(-a)) {
      for(int b=_b+1; b>0; b -= b&(-b)) {
        ret += bit[a][b];
      }
    }
    return ret;
  }
}; // END class BIT2D
