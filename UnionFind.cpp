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
