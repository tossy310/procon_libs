template<typename T>
class SparseTable {
public:
  vector<T> tbl[20];
  vector<T> tbr[20];
  int n;
  T (*op)(T, T);
  SparseTable(const vector<T> &v, T (*_op)(T, T)) : op(_op){
    // build table O(NlogN)
    n = v.size();
    tbl[0].resize(n+1);
    tbr[0].resize(n+1);
    rep(i,n) tbl[0][i] = tbr[0][i] = v[i];
    for(int k=1; (1<<k)<=n; k++){
      tbl[k].resize(n+1); tbr[k].resize(n+1);
      int mask = (1<<k)-1;
      rep(i,n){
        if((i&mask) == 0) tbr[k][i] = v[i];
        else tbr[k][i] = op(tbr[k][i-1], v[i]);
      }
      for(int i=n-1; i>=0; i--){
        if(((i+1)&mask) == 0) tbl[k][i] = v[i];
        else tbl[k][i] = op(tbl[k][i+1], v[i]);
      }
    }
  }
  T query(int l, int r){ // [l,r)
    if(l+1 == r) return tbl[0][l];
    r--; // [l,r]
    int k = 31 - __builtin_clz(l^r); // 2^k <= l^r < 2^{k+1}
    return op(tbl[k][l], tbr[k][r]);
  }
};