class UnionFind() :
    def __init__(self, size) :
        self.par = [-1 for _  in range(size)]
        self.rank = [0 for _ in range(size)]
    def find(self, x) :
        if self.par[x] < 0 :
            return x;
        else :
            self.par[x] = self.find(self.par[x])
            return self.par[x]
    def unite(self, x, y) :
        s1 = self.find(x)
        s2 = self.find(y)
        if s1 != s2:
            if self.rank[s1] < self.rank[s2] :
                self.par[s2] += self.par[s1]
                self.par[s1] = s2
            else:
                self.par[s1] += self.par[s2]
                self.par[s2] = s1
                if self.rank[s1] == self.rank[s2] :
                    self.rank[s1] += 1
        return
    def same(self, x, y) :
        return self.find(x) == self.find(y)
    def size(self, x) :
        return -self.par[self.find(x)]
