template<class T>
class RangeMap : public map<T, T> {
// this internally stores range [l,r).
public:
  auto get(T v) const {
    auto it = this->upper_bound(v);
    if(this->empty() || it == this->begin() || (--it)->second <= v) return this->end();
    return it;
  }
  auto add(T l, T r) { // [l,r), returns new range
    assert(l<r);
    auto itl = this->upper_bound(l), itr = this->upper_bound(r);
    if(itl != this->begin() && prev(itl)->second >= l) {
      itl--;
    }
    if(itl != itr) {
      l = min(l, itl->first);
      r = max(r, prev(itr)->second);
      this->erase(itl, itr);
    }
    return this->emplace(l, r).first;
  }
  // void remove(T l, T r) { /* TODO */ }
};
