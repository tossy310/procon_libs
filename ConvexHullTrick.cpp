// Convex Hull Trick for min query
template<typename T, typename U=__int128>
class ConvexHullTrick {
private:
  deque<pair<T,T>> data;
  bool check(const pair<T,T> &p3) const {
    auto &p1 = data[data.size()-2];
    auto &p2 = data[data.size()-1];
    return (U)(p2.first-p1.first)*(p3.second-p2.second) < (U)(p2.second-p1.second)*(p3.first-p2.first);
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
    return data[i].first*x + data[i].second;
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
    // care about overflow in next line
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
    if(size() == 0) { /* do something or return LLONG_MIN */ }
    auto l = *lower_bound((CHTLine){x, is_query});
    return l.m * x + l.b;
  }
};


template<class T>
class LiChaoTree {
private:
  struct Line {
    T a, b;
    Line(T _a, T _b) : a(_a), b(_b) {}
    inline T get(T x) const { return a*x + b; }
    inline bool over(const Line &lb, const T &x) const {
      return get(x) < lb.get(x);
    }
  };
  void update(Line &x, int k, int l, int r) {
    int mid = (l + r)/2;
    auto lover = x.over(seg[k], xs[l]), rover = x.over(seg[k], xs[mid]);
    if(rover) swap(seg[k], x);
    if(l+1 >= r) return;
    else if(lover != rover) update(x, 2*k+1, l, mid);
    else update(x, 2*k+2, mid, r);
  }
  T query(T x, int k) const {
    k += sz - 1;
    T ret = seg[k].get(x);
    while(k > 0) {
      k = (k - 1) / 2;
      ret = min(ret, seg[k].get(x));
    }
    return ret;
  }
public:
  vector<T> xs;
  vector<Line> seg;
  int sz = 0;
  LiChaoTree() {}
  LiChaoTree(const vector<T> &x, T _INF) : xs(x) {
    sz = 1;
    while(sz < (int)xs.size()) sz *= 2;
    while((int)xs.size() < sz) xs.push_back(xs.back() + 1);
    seg.assign(2*sz - 1, Line(0, _INF));
  }
  void update(T a, T b) { // ax+b
    if(sz == 0) return;
    Line l(a, b);
    update(l, 0, 0, sz);
  }
  T query(int k) const { // min-query to k-th x (xs[k])
    return query(xs[k], k);
  }
  T query_x(T x) const { // min-query x (actual value)
    int k = lower_bound(begin(xs), end(xs), x) - begin(xs);
    return query(x, k);
  }
};
