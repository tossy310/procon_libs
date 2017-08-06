vector<string> split(const string &str, char sep){
  vector<string> v;
  stringstream ss(str);
  string buffer;
  while( getline(ss, buffer, sep ) ) v.pb(buffer);
  if(str[str.size()-1] == sep) v.pb("");
  return v;
}
