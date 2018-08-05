const MOD = 1000000007;

function mul(x: number, y: number){
  if (x == 1) return y;
  if (y == 1) return x;
  if (x == 0 || y == 0) return 0;
  let ans = 0;
  while(y > 0){
    ans += x * (y%1000) % MOD;
    ans %= MOD;
    x = x * 1000 % MOD;
    y = Math.floor(y/1000);
  }
  return ans;
}
