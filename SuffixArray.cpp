// TODO SA-IS

/*
 sa[] SuffixArrayのi番目は，文字列sの何文字目から始まるか
  ** sa[0] = s.length
 lcp[] SuffixArray配列上のi番目とi+1番目の共通接頭文字列長
  ** lcp[0] = 0
 rank[] 文字列sのi番目の文字(0-indexed)は，sa上で何番目に配置されるか
  ** rank[s.length] = 0
 match(i,j) 文字列sのi文字目からのsuffixとj文字目からのsuffixの共通接頭文字列長
*/

class RMQ {
public:
  int n;
  vector<int> data;
  RMQ(){}
  RMQ(const vector<int> &lcp){
    int m = lcp.size();
    n=1;
    while(n<m) n*=2;
    data.resize(2*n, INF);
    rep(i,m) data[i+n] = lcp[i];
    for(int i=n-1; i>0; i--) data[i] = min(data[2*i], data[2*i+1]);
  }
  inline int query(int l, int r) const {
    int ret = INF;
    for(l+=n, r+=n; l<r; l/=2, r/=2){
      if(l&1) ret = min(ret, data[l++]);
      if(r&1) ret = min(ret, data[--r]);
    }
    return ret;
  }
};

class SuffixArray {
public:
  const string &s;
  vector<int> sa, lcp, rank;
  int n;
  RMQ rmq;
  void construct_sa(){
    sa.resize(n+1);
    rank.resize(n+1);
    rep(i,n+1){
      sa[i]=i;
      rank[i] = (i<n)?s[i]:(-1);
    }
    vector<int> tmp(n+1,0);
    for(int k=1; k<=n; k*=2){
      auto comp = [&](const int i, const int j){
        if(rank[i]!=rank[j]) return rank[i]<rank[j];
        int ri = (i+k<=n)?(rank[i+k]):-1;
        int rj = (j+k<=n)?(rank[j+k]):-1;
        return ri < rj;
      };
      sort(all(sa), comp);
      tmp[sa[0]]=0;
      rep(i,1,n+1){
        tmp[sa[i]] = tmp[sa[i-1]];
        if(comp(sa[i-1], sa[i])) tmp[sa[i]]++;
      }
      swap(tmp, rank);
    }
  }
  void construct_lcp(){
    rep(i,n+1) rank[sa[i]] = i;
    lcp.resize(n+1);
    int h=0;
    lcp[0]=0;
    rep(i,n){
      int j = sa[rank[i]-1];
      if(h>0) h--;
      for(; j+h<n && i+h<n; h++){
        if(s[j+h] != s[i+h]) break;
      }
      lcp[rank[i]-1]=h;
    }
  }
  void construct_rmq(){
    rmq = RMQ(lcp);
  }
  SuffixArray(const string &str) : s(str){
    n = s.size();
    construct_sa();
    construct_lcp();
    construct_rmq();
  }
  bool contain(const string &t) const {
    int m = t.size();
    int l = 0, r = n;
    while(r-l>1){
      int mid = (r+l)/2;
      if(s.compare(sa[mid], m, t) < 0) l = mid;
      else r = mid;
    }
    return s.compare(sa[r], m, t) == 0;
  }
  int lower_bound(const string &t) const {
    int m = t.size();
    int l = 0, r = n+1;
    while(r-l>1){
      int mid = (r+l)/2;
      if(s.compare(sa[mid], m, t) < 0) l = mid;
      else r = mid;
    }
    return r; // return rank
  }
  inline int match(int a, int b) const {
    int x = rank[a], y = rank[b];
    return rmq.query(min(x,y), max(x,y));
  }
  inline char operator[](int idx) const { return s[idx]; }
};
