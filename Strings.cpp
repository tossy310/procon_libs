// MP配列 : a[i] = s[0,i)の接頭辞と接尾辞が最大何文字一致しているか O(N)
vector<int> build_kmp(const string &s){
  int n = s.size();
  vector<int> a(n+1, -1);
  int j = -1;
  rep(i,n){
    while(j>=0 && s[i]!=s[j]) j = a[j];
    j++;
    a[i+1]=j;
  }
  return a;
}

// z-algorithm : a[i] = s と s[i,n) の最長共通接頭辞の長さ O(N)
vector<int> zAlgorithm(const string &s){
  int n = s.size();
  vector<int> a(n+1,0);
  int i=1, j=0;
  while(i<n){
    while(i+j<n && s[j] == s[i+j]) j++;
    a[i] = j;
    if(j==0){i++; continue;}
    int k = 1;
    while(i+k<n && k+a[k]<j) a[i+k]=a[k], k++;
    i += k;
    j -= k;
  }
  return a;
}

// manacher : a[i] = s[i] を中心とする最長の回文の半径 O(N)
// [TODO] verify
vector<int> manacher(const string &s){
  int n = s.size();
  int i=0, j=0;
  vector<int> a(n);
  while(i-j>=0 && i+j<n && s[i-j]==s[i+j]) j++;
  a[i]=j;
  int k = 1;
  while(i-k>=0 && i+k<n && k+a[i-k]<j) a[i+k]=a[i-k], k++;
  i += k;
  j -= k;
  return a;
}
// 奇数長の回文しか対応していないが，a$b$c$...みたいにすれば偶数長にも対応可
