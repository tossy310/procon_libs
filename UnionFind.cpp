// if par is negative, it is a root, of which tree size is the abs-value.
// rank represents depth
class UnionFind {
public:
  vector<int> par, rank;
  UnionFind(int sz) : par(sz, -1), rank(sz, 0){}
  int find(int x){
    if(par[x]<0) return x;
    else return par[x] = find(par[x]);
  }
  void unite(int x, int y){
    x=find(x); y=find(y);
    if(x==y) return;
    if(rank[x] < rank[y]) swap(x,y);
    par[x] += par[y];
    par[y] = x;
    if(rank[x]==rank[y]) rank[x]++;
  }
  inline bool same(int x, int y){ return find(x) == find(y); }
  inline int size(int x){ return -par[find(x)]; }
};


template <class T>
class WeightedUnionFind {
  int n;
  vector<int> par,sz;
  vector<T> w;
  // weight[i] : iがpar[i]よりどれだけ重いか
public:
  WeightedUnionFind(){}
  WeightedUnionFind(int _n) : n(_n) {
    par.resize(n, -1);
    sz.resize(n, 1);
    w.resize(n, 0);
  }
  int find(int i){
    if(par[i] < 0) return i;
    int p = find(par[i]);
    w[i] += w[par[i]];
    return par[i] = p;
  }
  T weight(int x){ find(x); return w[x]; }
  bool same(int x, int y){ return find(x) == find(y); }
  // weight[i] +w の位置に j を配置
  void unite(int i, int j, T nw){
    nw += weight(i) - weight(j);
    int x = find(i), y = find(j);
    if(x == y) return;
    if(sz[x] < sz[y]){
      swap(x, y);
      nw = -nw;
    }
    sz[x] += sz[y];
    par[y] = x;
    w[y] = nw;
  }
  // x からみた y の相対位置
  T diff(int x, int y){ return weight(y) - weight(x); }
};


class RewindableUnionFind {
public:
  vector<int> par, rank;
  vector<pair<int,pair<int,int>>> history; // <idx, <par, rank>>
  RewindableUnionFind(int sz) : par(sz, -1), rank(sz, 0){}
  int find(int x){
    if(par[x]<0) return x;
    else return find(par[x]);
  }
  void record(int x){
    history.pb(mp(x,mp(par[x],rank[x])));
  }
  void rewind(){
    assert(history.size()>=2);
    rep(_,2){
      auto &pp = history.back();
      par[pp.fi] = pp.se.fi;
      rank[pp.fi] = pp.se.se;
      history.pop_back();
    }
  }
  void unite(int x, int y){
    x = find(x); y = find(y);
    record(x); record(y);
    if(x==y) return;
    if(rank[x] < rank[y]) swap(x,y);
    par[x] += par[y];
    par[y] = x;
    if(rank[x]==rank[y]) rank[x]++;
  }
  inline bool same(int x, int y){ return find(x) == find(y); }
  inline int size(int x){ return -par[find(x)]; }
};
