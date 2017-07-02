struct edge {int to; long cost;};

vector<vector<edge> > EG;
vector<long> dist;

// dijkstra shortest path O(E logV)
// from:search starting node, d:array of distance, vec: directed graph
void dijkstra(int from, vector<long> &d, vector<vector<edge> > &vec){
  priority_queue<pair<long, int>, vector<pair<long, int> >, greater<pair<long, int> > > que; //pair<dist, idx>
  fill(all(d), INF);
  d[from] = 0;
  que.push(mp(0,from));
  while(!que.empty()){
    pair<long, int> p = que.top(); que.pop();
    int v = p.second;
    long dis = p.first;
    if(d[v] < dis) continue;
    for(edge e : vec[v]){
      long target = dis + e.cost;
      if(d[e.to] > target){
        d[e.to] = target;
        que.push(mp(d[e.to], e.to));
      }
    }
  }
}
// END dijkstra
