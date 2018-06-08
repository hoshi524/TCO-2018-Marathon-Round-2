#include <bits/stdc++.h>
#include <sys/time.h>
using namespace std;

class Timer {
 public:
  void restart();
  double getElapsed();

  Timer();

 private:
  static double rdtsc_per_sec_inv;

  double getTimeOfDay();
  unsigned long long int getCycle();

  double start_time;
  unsigned long long int start_clock;
};
double Timer::rdtsc_per_sec_inv = -1;

inline double Timer::getElapsed() {
  if (rdtsc_per_sec_inv != -1)
    return (double)(getCycle() - start_clock) * rdtsc_per_sec_inv;

  const double RDTSC_MEASUREMENT_INTERVAL = 0.1;
  double res = getTimeOfDay() - start_time;
  if (res <= RDTSC_MEASUREMENT_INTERVAL) return res;

  rdtsc_per_sec_inv = 1.0 / (getCycle() - start_clock);
  rdtsc_per_sec_inv *= getTimeOfDay() - start_time;
  return getElapsed();
}

inline void Timer::restart() {
  start_time = getTimeOfDay();
  start_clock = getCycle();
}

Timer::Timer() { restart(); }

inline double Timer::getTimeOfDay() {
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec * 0.000001;
}

inline unsigned long long int Timer::getCycle() {
  unsigned int low, high;
  __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
  return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
}

inline unsigned get_random() {
  static unsigned y = 2463534242;
  return y ^= (y ^= (y ^= y << 13) >> 17) << 5;
}
inline double get_random_double() { return get_random() / UINT32_MAX; }

/*
0     = none
1 ~6  = crystal
9 ~14 = lantern
16    = obstacle
17~18 = mirror
*/

constexpr double TIME_LIMIT = 2.5;
constexpr int N = 1 << 7;
constexpr int M = N * N;
constexpr int DIR[] = {1, N, -1, -N};
int H, W;
int CL, CM, CO, MM, MO;
int BOARD[M];

inline int to(int x, int y) { return (x << 7) | y; }
inline void to(int p, int& x, int& y) { x = p >> 7, y = p & (N - 1); }
inline bool in(int p) { return p > -1 && p < N * H && (p & (N - 1)) < W; }
inline int crystalScore(int c) { return bitset<3>(c).count() == 1 ? 20 : 30; }
int cost(int t) {
  if (t == 0) return 0;
  if (t < 16) return CL;
  if (t == 16) return CO;
  if (t > 16) return CM;
  assert(false);
}

struct State {
  int score = 0, mirrors = 0, obstacles = 0;
  int board[M];
  int light[M][4];

  void calcLight(int p) {
    int s[4];
    for (int i = 0; i < 4; ++i) {
      if (board[p] == 0) {
        s[i] = light[p][i];
      } else if (board[p] == 17) {
        s[i] = light[p][3 - i];
      } else if (board[p] == 18) {
        s[i] = light[p][i ^ 1];
      } else {
        s[i] = p;
      }
    }
    for (int d = 0; d < 4; ++d) {
      int a = p, b = d;
      while (true) {
        a += DIR[b];
        if (not in(a)) break;
        light[a][(b + 2) % 4] = s[d];
        if (board[a] == 0) continue;
        if (board[a] == 17) {
          b = 3 - b;
        } else if (board[a] == 18) {
          b = b ^ 1;
        } else {
          break;
        }
      }
    }
  }

  void calcLight() {
    memset(light, -1, sizeof(light));
    for (int i = 0; i < H; ++i) {
      for (int j = 0; j < W; ++j) {
        int p = to(i, j);
        if (board[p] == 0 or board[p] > 16) continue;
        calcLight(p);
      }
    }
  }

  int lightBit(int p) {
    int bit = 0;
    for (int k = 0; k < 4; ++k) {
      int n = light[p][k];
      if (n == -1) continue;
      int t = board[n];
      if (t > 8 && t < 16) bit |= t ^ 8;
    }
    return bit;
  }

  int lightBit(int p, int d, int t) {
    int x = light[p][d];
    light[p][d] = t;
    int b = lightBit(p);
    light[p][d] = x;
    return b;
  }

  void putItem(int p, int t) {
    board[p] = t;
    if (t == 16) obstacles++;
    if (t > 16) mirrors++;
    calcLight(p);
  }

  int diffScore(int p, int t) {
    int x = -cost(t);
    bool ok = false;
    if (t < 16) {
      int b = t ^ 8;
      for (int i = 0; i < 4; ++i) {
        int n = light[p][i];
        if (n == -1) continue;
        int u = board[n];
        if (u > 8 && u < 16) return -99;
        if (u < 8) {
          int bit = lightBit(n);
          if (bit == (bit | b)) continue;
          if (u == (u | b)) ok = true;
          if (u == (bit | b)) {
            int c = bitset<3>(u).count();
            x += c == 1 ? 20 : 40;
          } else if (u == bit) {
            x -= crystalScore(u) + 10;
          } else if (bit == 0) {
            x -= 10;
          }
        }
      }
    } else if (t == 16) {
    } else if (t > 16) {
    }
    return ok ? x : -99;
  }

  void calcScore() {
    score = 0, mirrors = 0, obstacles = 0;
    for (int i = 0; i < H; ++i) {
      for (int j = 0; j < W; ++j) {
        int p = to(i, j);
        if (board[p] > 0 and board[p] < 8) {
          int bit = lightBit(p);
          if (bit == 0) continue;
          if (bit == board[p]) {
            score += crystalScore(bit);
          } else {
            score -= 10;
          }
        }
        if (board[p] == BOARD[p])
          continue;
        else if (board[p] < 16)
          score -= CL;
        else if (board[p] == 16)
          score -= CO, obstacles++;
        else if (board[p] > 16)
          score -= CM, mirrors++;
      }
    }
  }

  void replace(int p1, int p2) {
    int h1, w1, h2, w2;
    to(p1, h1, w1);
    to(p2, h2, w2);
    for (int i = h1; i < h2; ++i) {
      for (int j = w1; j < w2; ++j) {
        int p = to(i, j);
        if (board[p] != BOARD[p]) board[p] = 0;
      }
    }
    calcLight();
    calcScore();
    static int edge[M * 4];
    int es = 0;
    for (int i = h1; i <= h2; ++i) {
      for (int j = w1; j <= w2; ++j) {
        int p = to(i, j);
        for (int k = 1; k < 8; k <<= 1) {
          int t = 8 | k;
          if (diffScore(p, t) > -99) edge[es++] = (p << 8) | t;
        }
      }
    }
    while (es > 0) {
      int i = get_random() % es;
      int p = edge[i] >> 8;
      int t = edge[i] & 0xff;
      edge[i] = edge[--es];
      if (board[p] != 0) continue;
      int v = diffScore(p, t);
      if (v > -99 and v > log(get_random_double())) {
        score += v;
        putItem(p, t);
        if (false) {
          int t = score;
          calcLight();
          calcScore();
          if (t != score) cerr << v << " " << t << " " << score << endl;
          assert(t == score);
        }
      }
    }
  }
};
State tmp, cur;

class CrystalLighting {
 public:
  vector<string> placeItems(vector<string>& targetBoard, int costLantern,
                            int costMirror, int costObstacle, int maxMirrors,
                            int maxObstacles) {
    Timer timer;
    {
      CL = costLantern;
      CM = costMirror;
      CO = costObstacle;
      MM = maxMirrors;
      MO = maxObstacles;
      H = targetBoard.size();
      W = targetBoard[0].size();
      for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
          char c = targetBoard[i][j];
          int k = to(i, j);
          if (c == '.') {
            BOARD[k] = 0;
          } else if (c == 'X') {
            BOARD[k] = 16;
          } else {
            BOARD[k] = c - '0';
          }
        }
      }
      memcpy(cur.board, BOARD, sizeof(BOARD));
    }
    {
      while (timer.getElapsed() < TIME_LIMIT) {
        memcpy(&tmp, &cur, sizeof(cur));
        constexpr int MASK = 5;
        int h = get_random() % (H - MASK);
        int w = get_random() % (W - MASK);
        tmp.replace(to(h, w), to(h + MASK, w + MASK));
        if (tmp.score - cur.score >= 0) {
          memcpy(&cur, &tmp, sizeof(tmp));
        }
      }
    }
    {
      vector<string> result;
      for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
          int k = to(i, j);
          if (cur.board[k] == BOARD[k]) continue;
          int x = cur.board[k];
          char c = '-';
          if (x < 16) {
            c = x - 8 + '0';
          } else if (x == 16) {
            c = 'X';
          } else if (x == 17) {
            c = '/';
          } else if (x == 18) {
            c = '\\';
          }
          ostringstream stream;
          stream << i << " " << j << " " << c;
          result.push_back(stream.str());
        }
      }
      return result;
    }
  }
};

// -------8<------- end of solution submitted to the website -------8<-------
int main() {
  CrystalLighting cl;
  int H;
  cin >> H;
  vector<string> t(H);
  for (int i = 0, s = t.size(); i < s; ++i) cin >> t[i];
  int a, b, c, d, e;
  cin >> a >> b >> c >> d >> e;
  vector<string> ret = cl.placeItems(t, a, b, c, d, e);
  cout << ret.size() << endl;
  for (int i = 0; i < (int)ret.size(); ++i) cout << ret[i] << endl;
  cout.flush();
}
