using ll = __int128;
string to_string(ll &i){
  string s;
  bool neg = false;
  if(i<0){ neg=true; i = -i; }
  while(i>0){ s.pb((i%10)+'0'); i/=10; }
  if(neg) s.pb('-');
  reverse(s.begin(), s.end());
  return s;
}
ll sto128(const string &s){
  bool neg = false;
  int p = 0;
  if(s[0]=='-'){ neg=true; p++; }
  else if(s[0]=='+'){ p++; }
  ll val = 0;
  while(p<(int)s.size()){
    if(!isdigit(s[p])) assert(false);
    val = val*10 + s[p++]-'0';
  }
  if(neg) val = -val;
  return val;
}
ostream& operator<<(ostream &o, ll i){
  o << to_string(i);
  return o;
}
istream& operator>>(istream &o, ll &i){
  string s;
  o >> s;
  i = sto128(s);
  return o;
}

// TODO bitset util
