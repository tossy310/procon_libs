// Convex Hull Trick for min query
template<typename T, typename U=__int128>
class ConvexHullTrick {
private:
  deque<pair<T,T>> data;
  bool check(const pair<T,T> &p3) const {
    auto &p1 = data[data.size()-2];
    auto &p2 = data[data.size()-1];
    return (U)(p2.fi-p1.fi)*(p3.se-p2.se) < (U)(p2.se-p1.se)*(p3.fi-p2.fi);
  }
public:
  ConvexHullTrick(){ }
  // add line ax + b
  // additions of 'a' shold be given in NON-INCREASING order
  void add(T a, T b){
    auto p = make_pair(a, b);
    while(data.size()>=2 && !check(p)) data.pop_back();
    data.push_back(p);
  }
  // get value of i-th data in deque with x
  inline T val(int i, T x) const {
    return data[i].fi*x + data[i].se;
  }
  // get minimun value of query x
  // queries shold be given in NON-DECREASING order
  T query(T x){
    while(data.size()>=2 && val(0, x) >= val(1, x)) data.pop_front();
    return val(0, x);
  }
  inline bool empty(){
    return data.empty();
  }
};

// no constraint for x (use vector instead of deque, deque is too slow!!)
T query(T x) const {
  int sz = data.size();
  int l = 0, r = sz;
  while(r-l>1){
    int m = (l+r)/2;
    if(val(m, x) - val(m-1, x) < 0) l = m;
    else r = m;
  }
  return val(l, x);
}



// Dynamic Convex Hull Trick
// TODO understand and verify
// from https://github.com/niklasb/contest-algos/blob/master/convex_hull/dynamic.cpp
const long is_query = -(1LL<<62);
struct CHTLine {
  long m, b;
  mutable function<const CHTLine*()> succ;
  bool operator<(const CHTLine& rhs) const {
    if(rhs.b != is_query) return m < rhs.m;
    const CHTLine* s = succ();
    if(!s) return 0;
    long x = rhs.m;
    return b - s->b < (s->m - m) * x;
  }
};
struct DynamicCHT : public multiset<CHTLine> {
  // will maintain upper hull for maximum
  bool bad(iterator y){
    auto z = next(y);
    if(y == begin()) {
      if(z == end()) return false;
      return y->m == z->m && y->b <= z->b;
    }
    auto x = prev(y);
    if(z == end()) return y->m == x->m && y->b <= x->b;
    return (x->b - y->b)*(z->m - y->m) >= (y->b - z->b)*(y->m - x->m);
  }
  void insert_line(long m, long b){
    auto y = insert({m, b});
    y->succ = [=] { return next(y) == end() ? 0 : &*next(y); };
    if(bad(y)){ erase(y); return; }
    while(next(y) != end() && bad(next(y))) erase(next(y));
    while(y != begin() && bad(prev(y))) erase(prev(y));
  }
  long eval(long x){
    auto l = *lower_bound((CHTLine){x, is_query});
    return l.m * x + l.b;
  }
};
