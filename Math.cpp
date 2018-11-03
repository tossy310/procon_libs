long large_primes[] = {1000000007L, 1000000009L, 1000000021L, 1000000033L, 1000000097L};
long larger_primes[] = {100000000000000003LL, 10000000000000000051LL};

const int MAX_PRIME = 300000;
vector<int> primes;
bool isPrime[MAX_PRIME+1];
// sieve of Erastosthenes
void prime_list(){
  const int sz = MAX_PRIME;
  fill(isPrime, isPrime+sz+1, true);
  isPrime[0]=isPrime[1]=false;
  for(int i=4;i<=sz;i+=2) isPrime[i]=false;
  primes.pb(2);
  for(int i=3; i<=sz; i+=2){
    if(isPrime[i]){
      primes.pb(i);
      for(long j=(long)i*i; j<=sz; j+=i) isPrime[j]=false;
    }
  }
}

// sieve of Erastosthenes
void prime_list(bool *arr, int sz){
  fill(arr, arr+sz+1, true);
  arr[0]=arr[1]=false;
  for(int i=4;i<=sz;i+=2) arr[i]=false;
  for(int i=3; i*i<=sz; i+=2){
    if(arr[i]) for(long j=(long)i*i; j<=sz; j+=i) arr[j]=false;
  }
}

// x^n mod p
long mod_pow(long x, long n, long p=MOD){
  if(x==0) return 0;
  long res=1;
  x %= p;
  while(n>0){
    if(n&1) res=res*x%p;
    x=x*x%p;
    n>>=1;
  }
  return res;
}
inline long mod_inv(long x, long p=MOD){ return mod_pow(x%p, p-2, p); }

long fact(long n, long p=MOD){
  long ret=1;
  while(n>1) ret=ret*n%p, n--;
  return ret;
}
inline long comb(long n, long k, long p=MOD){ return fact(n,p)*mod_inv(fact(k,p),p)%p*mod_inv(fact(n-k,p),p)%p; }

long comb[300][300];
void init_comb(){
  fill(comb[0],comb[300],0);
  comb[0][0] = 1;
  rep(i,1,300){
    comb[i][0] = 1;
    rep(j,i) comb[i][j+1] = (comb[i-1][j]+comb[i-1][j+1])%MOD;
  }
}

// a*x + b*y = gcd(a, b)
long extgcd(long a, long b, long& x, long& y){
  long d=a;
  if(b!=0){
    d=extgcd(b, a%b, y, x);
    y-=(a/b)*x;
  } else {
    x=1; y=0;
  }
  return d;
}

// Non-recursive extgcd
// a*x + b*y = gcd(a, b)
inline long extgcd(long a, long b, long& x, long& y){
  x = 1; y = 0;
  long tx = 0, ty = 1;
  while (b != 0){
    long q = a/b, r = a%b;
    long x2 = x - q*tx;
    long y2 = y - q*ty;
    a = b; b = r;
    x = tx; tx = x2;
    y = ty; ty = y2;
  }
  return a;
}

// Chinese Remainder Theorem
// u,v,mm are precomputable
// require gcd(m1, m2) == 1
// compute x (mod m1*m2) : x = a mod m1, x = b mod m2
inline long CRT(long a, long m1, long b, long m2){
  long u, v;
  extgcd(m1, m2, u, v);
  long mm = m1*m2;
  return ( (b*u %mm *m1 %mm) + (a*v %mm *m2 %mm) + mm) % mm;
}


// n!=a p^e, return a mod p O(log_p n)
long fact[MX + 5];
long mod_fact(long n, long p, long &e){
  e = 0;
  if(n == 0) return 1;
  long res = mod_fact(n/p, p, e);
  e += n/p;
  if(n/p %2 != 0) return res * (p - fact[n%p]) %p;
  return res * fact[n%p] %p;
}
// nCk mod p O(log_p n)
long mod_comb(long n, long k, long p) {
  if (n < 0 || k < 0 || n < k) return 0;
  long e1,e2,e3;
  long a1 = mod_fact(n,p,e1);
  long a2 = mod_fact(k,p,e2);
  long a3 = mod_fact(n-k,p,e3);
  if(e1 > e2+e3) return 0;
  return a1 * mod_inv(a2 * a3 % p, p) % p;
}


// Euler's totient function
// O(sqrt(n))
long eulerPhi(long n) {
  if(n==0) return 0;
  long ans = n;
  for(long x=2; x*x<=n; ++x) {
    if(n%x==0){
      ans -= ans/x;
      while(n%x==0) n/=x;
    }
  }
  if(n>1) ans -= ans/n;
  return ans;
}


// Cipolla's algorithm (modulous p should be a prime)
// calculate a : a^2 = x mod p
long mod_sqrt(long x, long p){
  x%=p;
  if(x==0) return 0;
  if(p==2) return x;
  auto is_square = [&](long a){ return mod_pow(a, (p-1)/2, p) == 1; };
  if(!is_square(x)) return -1;
  long y = 2;
  while(is_square((y*y-x+p)%p)) y++;
  long r = (y*y-x+p)%p;
  auto mul = [&](pair<long,long> u, pair<long,long> v){
    long s = (u.fi*v.fi + u.se*v.se%p*r)%p;
    long t = (u.se*v.fi + u.fi*v.se)%p;
    return mp(s,t);
  };
  long e = (p+1)/2;
  auto ret = mp(1LL,0LL);
  auto v = mp(y, 1LL);
  while(e>0){
    if(e&1) ret = mul(ret, v);
    v = mul(v, v);
    e /=2;
  }
  return ret.fi;
}

// combination nCr with no modulous
long comb(long n, long r){
  if(r>n-r) r=n-r;
  long p=1;
  for (long i=1; i<=r; i++) p=p*(n-i+1)/i;
  return p;
}


// combinarion nCr mod M 素因数分解を用いる 事前計算のprimesが必要
long comb(long n, long r, long mod){
  if(r>n-r) r=n-r;
  vector<int> a(n);
  rep(i,r) a[i]=n-i; //階乗表記での分子における積の表現に対応
  for(int p : primes){
    if(p>r) break;
    for(long q=p; q<=r; q*=p){
      for(int i=n%q, j=0; j<r/q; i+=q,j++){
        // r! のうち，qで割れる値はr/q個ある．また，a[i]はpでわりきれる
        a[i] /= p;
      }
    }
  }
  long mul=1;
  rep(i,r) mul = mul*a[i]%mod;
  return mul;
}


// x*y % mod p (for large x,y)
long mod_mul(long x, long y, long p){
  x %= p;
  y %= p;
  long res=0;
  while(y){
    if(y&1){
      res += x;
      if(res > p) res -= p;
    }
    x <<= 1;
    if(x>p) x -= p;
    y >>= 1;
  }
  return res;
}

// x^n % mod p (for large x)
long mod_pow(long x, long n, long p){
  long res=1;
  x %= p;
  while(n>0){
    if(n&1) res = mod_mul(res, x, p);
    x = mod_mul(x, x, p);
    n>>=1;
  }
  return res;
}


// Miller-Rabin Primality Test
bool millerRabin(long n){
  // using primes up to 37 is enough to test 64 bit integer
  const int primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
  const int nprime = 12;
  if(n<=1) return false;  // 1 or 0 or negative
  for(int i=0; i<nprime; i++){ // embeded primes and its multiple
    if(n%primes[i] == 0) return (n == primes[i]);
  }
  long m = n-1, s=0;
  while(m%2==0){ s++; m>>=1; } // n=m*2^s+1
  // Miller-Rabin primality test
  for(int i=0; i<nprime; i++){
    long a = primes[i];
    long y = mod_pow(a, m, n);
    if(y==1) continue;
    bool flg=true;
    for(int j=0; flg && (j<s); j++){
      if(y==n-1) flg=false;
      y = mod_mul(y, y, n);
    }
    if(flg==false) continue;
    return false;
  }
  return true;
}


//Pollard's Rho Algorithm  (find a divisor of n)
long pollardRho(long n, int c=1){
  long x=2, y=2, d=1;
  while(d==1){
    x = mod_mul(x, x, n)+c;
    y = mod_mul(y, y, n)+c;
    y = mod_mul(y, y, n)+c;
    d = __gcd( ((x>y)?(x-y):(y-x)), n);
  }
  if(d==n) return pollardRho(n, c+1);
  else return d;
}


// =====================素因数分解ここから
const int MAX_PRIME = 300000;
vector<int> primes;
vector<bool> isPrime(MAX_PRIME+1, true);

// sieve of Erastosthenes
inline void prime_list(){
  const int sz = MAX_PRIME;
  isPrime[0]=isPrime[1]=false;
  for(int i=4;i<=sz;i+=2) isPrime[i]=false;
  primes.pb(2);
  for(int i=3; i<=sz; i+=2){
    if(isPrime[i]){
      primes.pb(i);
      for(long j=(long)i*i; j<=sz; j+=i) isPrime[j]=false;
    }
  }
}

bool is_prime(long n){
  if( n <= MAX_PRIME) return isPrime[n];
  else return millerRabin(n);
}

void factorize(long n, map<long, int> & factors){
  if(is_prime(n)){
    factors[n]++;
    return;
  }
  for(int i=0; i<primes.size(); i++){
    int p=primes[i];
    while(n%p==0){
      factors[p]++;
      n /= p;
    }
  }
  if(n != 1){
    if(is_prime(n)){
      factors[n]++;
    } else {
      long d = pollardRho(n);
      factorize(d, factors);
      factorize(n/d, factors);
    }
  }
}
//==================素因数分解ここまで
// mainの最初でprime_list()を呼ぶこと！
// Miller-RabinとPollard-Rhoが必要．





// 埋め込み階乗MOD計算
const vector<long> arr = {1,682498929,491101308,76479948,723816384,67347853,27368307,625544428,199888908,888050723,927880474,281863274,661224977,623534362,970055531,261384175,195888993,66404266,547665832,109838563,933245637,724691727,368925948,268838846,136026497,112390913,135498044,217544623,419363534,500780548,668123525,128487469,30977140,522049725,309058615,386027524,189239124,148528617,940567523,917084264,429277690,996164327,358655417,568392357,780072518,462639908,275105629,909210595,99199382,703397904,733333339,97830135,608823837,256141983,141827977,696628828,637939935,811575797,848924691,131772368,724464507,272814771,326159309,456152084,903466878,92255682,769795511,373745190,606241871,825871994,957939114,435887178,852304035,663307737,375297772,217598709,624148346,671734977,624500515,748510389,203191898,423951674,629786193,672850561,814362881,823845496,116667533,256473217,627655552,245795606,586445753,172114298,193781724,778983779,83868974,315103615,965785236,492741665,377329025,847549272,698611116};
const vector<long> vec = { 0,10000000,20000000,30000000,40000000,50000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000,130000000,140000000,150000000,160000000,170000000,180000000,190000000,200000000,210000000,220000000,230000000,240000000,250000000,260000000,270000000,280000000,290000000,300000000,310000000,320000000,330000000,340000000,350000000,360000000,370000000,380000000,390000000,400000000,410000000,420000000,430000000,440000000,450000000,460000000,470000000,480000000,490000000,500000000,510000000,520000000,530000000,540000000,550000000,560000000,570000000,580000000,590000000,600000000,610000000,620000000,630000000,640000000,650000000,660000000,670000000,680000000,690000000,700000000,710000000,720000000,730000000,740000000,750000000,760000000,770000000,780000000,790000000,800000000,810000000,820000000,830000000,840000000,850000000,860000000,870000000,880000000,890000000,900000000,910000000,920000000,930000000,940000000,950000000,960000000,970000000,980000000,990000000,1000000000};
long embeded_fact(long n){
  if(n>=MOD){ cout<<0<<endl; return 0;}
  int i = upper_bound(all(vec), n) - vec.begin();
  i--;
  long ret = arr[i];
  i = vec[i];
  i++;
  while(i<=n){ret = (ret*i)%MOD; i++;}
  return ret;
}
