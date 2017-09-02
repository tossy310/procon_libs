// Bellman Ford O(VE)
// 閉路検出時はfalseを返す
// 入力は隣接リストのグラフ表現
template<class T>
bool bellmanFord(const vector<vector<pair<int,T>>> &vec, vector<T> &d, int from=0){
  fill(all(d), INF);
  d[from] = 0;
  int n = d.size();
  rep(_,n){
    bool update = false;
    rep(i,n){
      for(auto &p : vec[i]){
        int to = p.first;
        T nd = d[i] + p.second;
        if(d[to] > nd){
          d[to] = nd;
          update = true;
        }
      }
    }
    if(!update) return true;
    if(d[from]<0) return false;
  }
  return false;
}


// 辺集合を入力に受け取る
struct edge {int from, to; long cost;};
bool bellmanFord(int from, vector<long> &d, vector<edge> &es){
  fill(all(d), INF);
  d[from]=0;
  int n=d.size();
  rep(i,n){
    bool update=false;
    for(int j=0; j<es.size(); j++){
      edge &e = es[j];
      long target = d[e.from]+e.cost;
      if(d[e.from]!=INF && d[e.to]>target){
        d[e.to] = target;
        update=true;
      }
    }
    if(!update) return true;
    if(d[from]<0) return false;
  }
  return false;
}
