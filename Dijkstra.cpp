// weighted
template<class T>
void dijkstra(const vector<vector<pair<int,T>>> &vec, vector<T> &d, int from=0){
  using P = pair<T,int>;
  fill(all(d), INF);
  priority_queue<P, vector<P>, greater<P>> pq;
  d[from] = 0;
  pq.push(mp(0,from));
  while(!pq.empty()){
    auto p = pq.top(); pq.pop();
    int v = p.second;
    T dd = p.first;
    if(d[v] < dd) continue;
    for(auto to : vec[v]){
      T nd = dd + to.second;
      int ni = to.first;
      if(d[ni] > nd){
        d[ni] = nd;
        pq.push(mp(nd, ni));
      }
    }
  }
}

// TODO unweighted 01 dijkstra
