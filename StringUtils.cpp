vector<string> split(const string &str, char sep){
  vector<string> v;
  stringstream ss(str);
  string buffer;
  while( getline(ss, buffer, sep ) ) v.pb(buffer);
  if(str[str.size()-1] == sep) v.pb("");
  return v;
}

template <class T>
string join(const vector<T> &v, const char* delim = ""){
  ostringstream os;
  for(auto &i : v){ os << i << delim; }
  string str = os.str();
  str.erase(str.size() - strlen(delim));
  return str;
}
template <class T>
string join(const vector<T> &v, const char delim){
  return join(v, string(1, delim).c_str());
}
