// とりあえずACコードのコピペを貼っておく．
// TODO: classとしてのインタフェース決定

// === global ===
long pow1[100005], pow2[100005];
long hashs1[100005], hashs2[100005], hasht1[100005], hasht2[100005];

// === loal ===
string s,t;
cin>>s>>t;
int n = s.size();
int m = t.size();

pow1[0]=1;
pow2[0]=1;
rep(i,100000) pow1[i+1] = pow1[i]*MOD % BASE;
rep(i,100000) pow2[i+1] = pow2[i]*DOM % BASE;

auto hashinit = [&](long *arr, string &str, const long mod){
  int l = str.size();
  arr[0] = str[0];
  repl(i,1,l) arr[i] = (arr[i-1]*mod + str[i])%BASE;
};

hashinit(hashs1, s, MOD);
hashinit(hashs2, s, DOM);
hashinit(hasht1, t, MOD);
hashinit(hasht2, t, DOM);

auto gethash = [&](long *arr, int from, int l, long *powarr){
  if(from==0) return arr[l-1];
  long ret = arr[from+l-1] - arr[from-1]*powarr[l];
  ret %= BASE;
  if(ret<0) ret += BASE;
  return ret;
};

// s[si], t[ti] からはじまるl文字
auto hasheq = [&](int si, int ti, int l){
  long h1 = gethash(hashs1, si, l, pow1);
  long h2 = gethash(hasht1, ti, l, pow1);
  if(h1!=h2) return false;
  long h3 = gethash(hashs2, si, l, pow2);
  long h4 = gethash(hasht2, ti, l, pow2);
  return h3==h4;
};

auto eqlen = [&](int si, int ti){
  if(s[si]!=t[ti]) return 0;
  int l=1, r=min(n-si, m-ti)+1;
  while(r-l>1){
    int mid = (l+r)/2;
    if(hasheq(si,ti,mid)) l = mid;
    else r = mid;
  }
  return l;
};
