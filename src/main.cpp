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
inline double get_random_double() { return (double)get_random() / UINT32_MAX; }

/*
0     = none
1 ~6  = crystal
9 ~14 = lantern
16    = obstacle
17~18 = mirror
*/

constexpr double TIME_LIMIT = 2.8;
constexpr int N = 1 << 7;
constexpr int M = N * N;
constexpr int DIR[] = {1, N, -1, -N};
constexpr int MINI = -0xffff;
int H, W;
int CL, CM, CO, MM, MO;
int BOARD[M];
double remain;

inline int to(int x, int y) { return (x << 7) | y; }
inline void to(int p, int& x, int& y) { x = p >> 7, y = p & (N - 1); }
inline bool in(int p) { return p > -1 && p < N * H && (p & (N - 1)) < W; }
inline int crystalScore(int c) { return bitset<3>(c).count() == 1 ? 20 : 30; }
inline bool isC(int t) { return 0 < t && t < 8; }
inline bool isL(int t) { return 8 < t && t < 16; }
inline bool isO(int t) { return t == 16; }
inline bool isM(int t) { return 16 < t && t < 19; }
int cost(int t) {
  if (t == 0) return 0;
  if (t < 16) return CL;
  if (t == 16) return CO;
  if (t > 16) return CM;
  assert(false);
}

struct State {
  int score1 = 0, score2 = 0, mirrors = 0, obstacles = 0;
  double score = 0;
  int board[M];
  int light[M][4];

  void calcLight(int p, bool on = true) {
    for (int d = 0; d < 4; ++d) {
      int a = p, b = d;
      int l = [&]() {
        if (not on) return -1;
        if (board[p] == 0) {
          return light[p][d];
        } else if (board[p] == 17) {
          return light[p][3 - d];
        } else if (board[p] == 18) {
          return light[p][d ^ 1];
        } else {
          return p;
        }
      }();
      while (true) {
        a += DIR[b];
        if (not in(a)) break;
        light[a][(b + 2) % 4] = l;
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
        if (board[p] > 0 && board[p] <= 16) calcLight(p);
      }
    }
  }

  int lightBit(int p) {
    int bit = 0;
    for (int k = 0; k < 4; ++k) {
      int n = light[p][k];
      if (n == -1) continue;
      int t = board[n];
      if (isL(t)) bit |= t ^ 8;
    }
    return bit;
  }

  int lightBit(int p, int i, int j) {
    static int tmp[4];
    memcpy(tmp, light[p], sizeof(tmp));
    for (int k = 0; k < 4; ++k) {
      if (light[p][k] == i) {
        light[p][k] = j;
        break;
      }
    }
    int b = lightBit(p);
    memcpy(light[p], tmp, sizeof(tmp));
    return b;
  }

  int calcScore1(int t, int x) {
    if (x == 0) return 0;
    if (x == t) return crystalScore(x);
    return -10;
  }

  int calcScore2(int t, int x) {
    if (x == 0) return 0;
    if (x == t) return crystalScore(x);
    if (t == (t | x)) return 10;
    return -10;
  }

  double calcScore(int v1, int v2) {
    return (score1 + v1) * (1.0 - remain) + (score2 + v2) * remain;
  }

  void putItem(int p, int t) {
    assert(board[p] == 0);
    assert(t > 0);
    static int _light[4];
    if (isM(t)) {
      mirrors++;
      memcpy(_light, light[p], sizeof(_light));
      for (int i = 0; i < 4; ++i) {
        int x = _light[i];
        if (x != -1) {
          calcLight(x, false);
        }
      }
    }
    board[p] = t;
    if (t == 16) obstacles++;
    if (t <= 16) {
      calcLight(p);
    } else {
      for (int i = 0; i < 4; ++i) {
        int x = _light[i];
        if (x != -1) calcLight(x);
      }
    }
  }

  tuple<int, int> diffScore(int p, int t) {
    tuple<int, int> invalid = forward_as_tuple(MINI, MINI);
    if (board[p] != 0) return invalid;
    int v1 = -cost(t), v2 = -cost(t);
    auto add = [&](int t, int cur, int nxt) {
      v1 += calcScore1(t, nxt) - calcScore1(t, cur);
      v2 += calcScore2(t, nxt) - calcScore2(t, cur);
    };
    bool ok = false;
    auto uniqueC = [&](int n, int i) {
      for (int j = 0; j < i; ++j) {
        if (n == light[p][j]) return false;
      }
      return true;
    };
    if (t < 16) {
      int b = t ^ 8;
      for (int i = 0; i < 4; ++i) {
        int n = light[p][i];
        if (n == -1) continue;
        int u = board[n];
        if (isL(u)) return invalid;
        if (isC(u) && uniqueC(n, i)) {
          int bit = lightBit(n);
          if (bit == (bit | b)) continue;
          if (u == (u | b)) ok = true;
          add(u, bit, bit | b);
        }
      }
    } else if (t == 16) {
      auto diff = [&](int d1, int d2) {
        int p1 = light[p][d1];
        int p2 = light[p][d2];
        if (p1 == -1) return;
        if (p2 == -1) return;
        int t1 = board[p1];
        int t2 = board[p2];
        if (isC(t1) && isL(t2)) {
          if ((t1 & t2) == 0) ok = true;
          add(t1, lightBit(p1), lightBit(p1, p2, -1));
        }
        if (isC(t2) && isL(t1)) {
          if ((t1 & t2) == 0) ok = true;
          add(t2, lightBit(p2), lightBit(p2, p1, -1));
        }
      };
      diff(0, 2), diff(1, 3);
    } else if (t > 16) {
      for (int i = 0; i < 4; i += 2) {
        int d1, d2;
        if (t == 17) {
          d1 = i, d2 = i ^ 1;
        } else {
          d1 = i, d2 = 3 - i;
        }
        int p1 = light[p][d1];
        int p2 = light[p][d2];
        if (p1 == -1 && p2 == -1) continue;
        int t1 = p1 == -1 ? 16 : board[p1];
        int t2 = p2 == -1 ? 16 : board[p2];
        if (isL(t1) && isL(t2)) return invalid;
        if (isC(t1)) {
          if (isL(t2) && (t1 & t2)) ok = true;
          add(t1, lightBit(p1), lightBit(p1, light[p][(d1 + 2) % 4], p2));
        }
        if (isC(t2)) {
          if (isL(t1) && (t1 & t2)) ok = true;
          add(t2, lightBit(p2), lightBit(p2, light[p][(d2 + 2) % 4], p1));
        }
      }
      if (ok) {
        for (int i = 0; i < 4; ++i) {
          int t1 = board[light[p][i]];
          if (not isC(t1)) continue;
          for (int j = 0; j < i; ++j) {
            int t2 = board[light[p][j]];
            if (t1 == t2) return invalid;
          }
        }
      }
    }
    return ok ? forward_as_tuple(v1, v2) : invalid;
  }

  void calcScore() {
    score1 = 0, score2 = 0, mirrors = 0, obstacles = 0;
    for (int i = 0; i < H; ++i) {
      for (int j = 0; j < W; ++j) {
        int p = to(i, j);
        int t = board[p];
        if (isC(t)) {
          int bit = lightBit(p);
          if (bit == 0) continue;
          score1 += calcScore1(t, bit);
          score2 += calcScore2(t, bit);
        }
        if (t == BOARD[p])
          continue;
        else if (isL(t))
          score1 -= CL, score2 -= CL;
        else if (isO(t))
          score1 -= CO, score2 -= CO, obstacles++;
        else if (isM(t))
          score1 -= CM, score2 -= CM, mirrors++;
      }
    }
    score = calcScore(0, 0);
  }

  void replace(int p1, int p2) {
    int h1, w1, h2, w2;
    to(p1, h1, w1);
    to(p2, h2, w2);
    {  // remove
      for (int i = h1; i < h2; ++i) {
        for (int j = w1; j < w2; ++j) {
          int p = to(i, j);
          int t = board[p];
          if (t != BOARD[p]) board[p] = 0;
        }
      }
      calcLight();
      for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
          int p = to(i, j);
          int t = board[p];
          if (isL(t)) {
            for (int k = 0; k < 4; ++k) {
              int n = light[p][k];
              if (n != -1 && isL(board[n])) board[n] = 0;
            }
          }
        }
      }
      calcLight();
      calcScore();
    }
    static int edge[M * 4];
    int es = 0;
    for (int i = h1; i <= h2; ++i) {
      for (int j = w1; j <= w2; ++j) {
        int p = to(i, j);
        if (board[p] != 0) continue;
        for (int k = 1; k < 8; k <<= 1) {
          int t = 8 | k;
          edge[es++] = (p << 8) | t;
        }
      }
    }
    while (es > 0) {
      double _score = -1e10;
      int p = 0, t = 0, v1 = 0, v2 = 0;
      if ((obstacles < MO || mirrors < MM) && get_random() % 100 == 0) {
        for (int i = h1; i <= h2; ++i) {
          for (int j = w1; j <= w2; ++j) {
            int tp = to(i, j);
            if (board[tp] != 0) continue;
            for (int tt = 16; tt < 19; ++tt) {
              if (tt == 16 && obstacles == MO) continue;
              if (tt > 16 && mirrors == MM) continue;
              int _v1, _v2;
              tie(_v1, _v2) = diffScore(tp, tt);
              if (_v1 == MINI) continue;
              double x = calcScore(_v1, _v2);
              if (_score < x) {
                _score = x;
                p = tp;
                t = tt;
                v1 = _v1;
                v2 = _v2;
              }
            }
          }
        }
      } else {
        int i = get_random() % es;
        p = edge[i] >> 8;
        t = edge[i] & 0xff;
        edge[i] = edge[--es];
        if (board[p] != 0) continue;
        tie(v1, v2) = diffScore(p, t);
        if (v1 == MINI) continue;
        _score = calcScore(v1, v2);
      }
      if (_score - score > remain * log(get_random_double())) {
        score = _score;
        score1 += v1;
        score2 += v2;
        putItem(p, t);
        if (false) {
          int p1 = score1;
          int p2 = score2;
          static int tmp[M][4];
          memcpy(tmp, light, sizeof(tmp));
          calcLight();
          calcScore();
          assert(p1 == score1);
          assert(p2 == score2);
          for (int i = 0; i < M; ++i)
            for (int j = 0; j < 4; ++j) assert(tmp[i][j] == light[i][j]);
        }
      }
    }
  }
};
State cur;

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
    State tmp, bst;
    {
      cur.replace(to(0, 0), to(H, W));
      while (true) {
        remain = 1.0 - timer.getElapsed() / TIME_LIMIT;
        if (remain < 0) break;
        memcpy(&tmp, &cur, sizeof(cur));
        constexpr int MASK = 6;
        int h = get_random() % (H - MASK);
        int w = get_random() % (W - MASK);
        tmp.replace(to(h, w), to(h + MASK, w + MASK));
        cur.score = cur.calcScore(0, 0);
        if (tmp.score - cur.score > 5 * remain * log(get_random_double())) {
          memcpy(&cur, &tmp, sizeof(tmp));
        }
        if (bst.score1 < tmp.score1) {
          memcpy(&bst, &tmp, sizeof(tmp));
        }
      }
    }
    {
      vector<string> result;
      for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
          int k = to(i, j);
          int x = bst.board[k];
          if (x == BOARD[k]) continue;
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
