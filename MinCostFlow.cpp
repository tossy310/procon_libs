template<typename T>
class MinCostFlow{
private:
  struct edge{int to; T cap, cost; int rev;};
  using P = pair<T,int>;
  vector<vector<edge> > Graph;
  vector<int> prevv, preve;
  vector<T> h, d; // ポテンシャル，最短距離
public:
  MinCostFlow(int v){
    // 頂点数vで初期化
    Graph.resize(v);
    prevv.resize(v);
    preve.resize(v);
    h.resize(v);
    d.resize(v);
  }
  T min_cost_flow(int s, int t, T f){
    T res = 0;
    fill(all(h), 0);
    // 負の辺除去が必要なとき
    // rep(v,Graph.size()){
    //   rep(j,Graph[v].size()){
    //     edge &e = Graph[v][j];
    //     if(e.cap==0) continue;
    //     int u = e.to;
    //     h[u] = min(h[u],h[v]+e.cost);
    //   }
    // }
    while(f>0){
      priority_queue<P, vector<P>, greater<P>> pq;
      fill(all(d), INF);
      d[s] = 0;
      pq.push(mp(0,s));
      while(!pq.empty()){
        auto p = pq.top(); pq.pop();
        int v = p.se;
        if(d[v] < p.fi) continue;
        rep(i,Graph[v].size()){
          edge &e = Graph[v][i];
          if(e.cap > 0 && d[e.to] > d[v] + e.cost + h[v] - h[e.to]){
            d[e.to] = d[v] + e.cost + h[v] - h[e.to];
            prevv[e.to] = v;
            preve[e.to] = i;
            pq.push(mp(d[e.to], e.to));
          }
        }
      }
      if(d[t] == INF) return -1;
      rep(i,Graph.size()) h[i] += d[i];

      T nf = f;
      for(int v=t; v!=s; v = prevv[v]){
        nf = min(nf, Graph[prevv[v]][preve[v]].cap);
      }
      f -= nf;
      res += nf * h[t];
      for(int v=t; v!=s; v=prevv[v]){
        edge &e = Graph[prevv[v]][preve[v]];
        e.cap -= nf;
        Graph[v][e.rev].cap += nf;
      }
    }
    return res;
  }
  void add_edge(int from ,int to, T cap, T cost){
    Graph[from].pb(((edge){to, cap, cost, (int)Graph[to].size()}));
    Graph[to].pb(((edge){from, 0, -cost, (int)Graph[from].size()-1}));
  }
};
