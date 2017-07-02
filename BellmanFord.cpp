struct edge {int from, to; long cost;};

vector<long> dist;
vector<edge> EG;

// Bellman-Ford Shortest Path O(V * E)
// 負の閉路を検出したらfalseを返す
bool bellmanFord(int from, vector<long> &d, vector<edge> &es){
  fill(all(d), INF);
  d[from]=0;
  int n=d.size();
  rep(i,n){
    bool update=false;
    for(int j=0; j<es.size(); j++){
      edge e = es[j];
      long target = d[e.from]+e.cost;
      if(d[e.from]!=INF && d[e.to]>target){
        d[e.to] = target;
        update=true;
      }
    }
    if(!update) return true;
  }
  return false;
}
