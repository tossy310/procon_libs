class UnionFind {
public:
  vector<int> par, rank; // parent(negative := its root and abs-value is its size), depth
  UnionFind(int sz) : par(sz, -1), rank(sz, 0){}
  int find(int x){
    if(par[x]<0) return x;
    else return par[x] = find(par[x]);
  }
  void unite(int x, int y){
    x=find(x); y=find(y);
    if(x==y) return;  // already belong to same tree
    if(rank[x] < rank[y]){  // y becomes parent node
      par[y] += par[x];
      par[x] = y;
    } else {  // x becomes parent node
      par[x] += par[y];
      par[y] = x;
      if(rank[x]==rank[y]) rank[x]++;
    }
  }
  inline bool same(int x, int y){ return find(x) == find(y); }
  inline int size(int x){ return -par[find(x)]; }
};
