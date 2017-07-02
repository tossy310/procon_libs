// 蟻本写経しただけでverityしてないので注意！
class SuffixArray {
public:
  string &str;
  vector<int> sa, lcp;
  int n;
  SuffixArray(string &s) : str(s){
    n = str.size();
  }
  void construct_sa(){
    sa.resize(n+1);
    vector<int> rank(n+1);
    rep(i,n+1){
      sa[i]=i;
      rank[i] = (i<n)?str[i]:(-1);
    }
    for(int k=1; k<=n; k*=2){
      function<bool(int,int)> comp = [&](int i, int j){
        if(rank[i]!=rank[j]) return rank[i]<rank[j];
        int ri = (i+k<=n)?(rank[i+k]):-1;
        int rj = (j+k<=n)?(rank[j+k]):-1;
        return ri < rj;
      };
      sort(all(sa), comp);
      vector<int> tmp(n+1,0);
      tmp[sa[0]]=0;
      repl(i,1,n+1){
        tmp[sa[i]] = tmp[sa[i-1]];
        if(comp(sa[i-1], sa[i])) tmp[sa[i]]++;
      }
      swap(tmp, rank);
    }
  }
  void construct_lcp(){
    vector<int> rank(n+1);
    rep(i,n+1) rank[sa[i]] = i;
    lcp.resize(n+1);
    int h=0;
    lcp[0]=0;
    rep(i,n){
      int j = sa[rank[i]-1];
      if(h>0) h--;
      for(; j+h<n && i+h<n; h++){
        if(str[j+h] != str[i+h]) break;
      }
      lcp[rank[i]-1]=h;
    }
  }
  void construct(){
    construct_sa();
    construct_lcp();
  }
};
