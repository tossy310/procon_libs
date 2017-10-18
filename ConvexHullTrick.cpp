// Convex Hull Trick for min query
template<typename T, typename U=__int128>
class ConvexHullTrick {
private:
  deque<pair<T,T>> data;
  bool check(const pair<T,T> &p3){
    auto &p1 = data[data.size()-2];
    auto &p2 = data[data.size()-1];
    return (U)(p2.fi-p1.fi)*(p3.se-p2.se) < (U)(p2.se-p1.se)*(p3.fi-p2.fi);
  }
public:
  ConvexHullTrick(){ }
  // additions of 'a' shold be NON-INCREASING
  void add(T a, T b){
    auto p = mp(a,b);
    while(data.size()>=2 && !check(p)) data.pop_back();
    data.push_back(p);
  }
  inline T val(int i, T x) {
    return data[i].fi*x + data[i].se;
  }
  // queries shold be NON-DECREASING
  T query(T x){
    while(data.size()>=2 && val(0, x) >= val(1, x)) data.pop_front();
    return val(0, x);
  }
  inline bool empty(){
    return data.empty();
  }
};

// TODO dynamic convex hull Trick
// cf. https://github.com/niklasb/contest-algos/blob/master/convex_hull/dynamic.cpp
