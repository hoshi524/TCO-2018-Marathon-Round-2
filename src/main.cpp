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

constexpr double TIME_LIMIT = 1.8;
constexpr int N = 1 << 7;
constexpr int M = N * N;
constexpr int DIR[] = {1, N, -1, -N};
constexpr double invalid = -1e10;
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
inline int rev(int d) { return (d + 2) % 4; }
int cost(int t) {
  if (isL(t)) return CL;
  if (isO(t)) return CO;
  if (isM(t)) return CM;
  return 0;
}

struct State {
  int score1 = 0, score2 = 0, mirrors = 0, obstacles = 0;
  int board[M];
  int light[M][4];

  void init() {
    calcLight();
    calcScore();
    replace(to(0, 0), to(H, W));
  }

  int bend(int p, int d) { return bend(p, d, board[p]); }

  int bend(int p, int d, int t) {
    if (t == 0) {
      return light[p][rev(d)];
    } else if (t == 17) {
      return light[p][d ^ 1];
    } else if (t == 18) {
      return light[p][3 - d];
    } else {
      return p;
    }
  }

  void calcLight(int p) {
    for (int d = 0; d < 4; ++d) {
      int a = p, b = d, l = bend(p, d);
      while (true) {
        a += DIR[b];
        if (not in(a)) break;
        light[a][(b + 2) % 4] = l;
        if (p == a) break;
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

  double score(int v1 = 0, int v2 = 0) {
    return (score1 + v1) * (1.0 - remain) + (score2 + v2) * remain;
  }

  void putItem(int p, int t) {
    static bool ok[4];
    for (int i = 0; i < 4; ++i) {
      ok[i] = [&]() {
        int a = light[p][i];
        if (a == -1) return false;
        int t = board[a];
        if (!isC(t)) return false;
        for (int j = 0; j < i; ++j)
          if (light[p][i] == light[p][j]) return false;
        return true;
      }();
    }
    auto calc = [&](bool add) {
      for (int i = 0; i < 4; ++i) {
        if (ok[i]) {
          int a = light[p][i];
          int b = lightBit(a);
          if (add) {
            score1 += calcScore1(board[a], b);
            score2 += calcScore2(board[a], b);
          } else {
            score1 -= calcScore1(board[a], b);
            score2 -= calcScore2(board[a], b);
          }
        }
      }
      int t = board[p];
      int c = cost(t);
      if (add) {
        score1 -= c, score2 -= c;
        if (isO(t)) obstacles++;
        if (isM(t)) mirrors++;
      } else {
        score1 += c, score2 += c;
        if (isO(t)) obstacles--;
        if (isM(t)) mirrors--;
      }
    };
    calc(false);
    board[p] = t;
    calcLight(p);
    calc(true);
  }

  double diffScore(int p, int t) {
    int pt = board[p];
    if (t == pt) return invalid;
    int v1 = cost(pt) - cost(t);
    int v2 = cost(pt) - cost(t);
    static bool used[4];
    memset(used, false, sizeof(used));
    for (int d = 0; d < 4; ++d) {
      if (used[d]) continue;
      int cp = light[p][d];
      if (cp == -1) continue;
      int ct = board[cp];
      if (!isC(ct)) continue;
      {
        int cb = lightBit(cp);
        v1 -= calcScore1(ct, cb);
        v2 -= calcScore2(ct, cb);
      }
      static int dl[4];
      int ds = 0;
      for (int i = d; i < 4; ++i) {
        if (cp == light[p][i]) {
          dl[ds++] = i;
          used[i] = true;
        }
      }
      {
        static int tmp[4];
        memcpy(tmp, light[cp], sizeof(tmp));
        for (int i = 0; i < ds; ++i) {
          int pp = bend(p, dl[i], pt);
          int np = bend(p, dl[i], t);
          for (int j = 0; j < 4; ++j) {
            if (light[cp][j] == pp) {
              light[cp][j] = np;
              break;
            }
          }
        }
        board[p] = t;
        int nb = lightBit(cp);
        board[p] = pt;
        memcpy(light[cp], tmp, sizeof(tmp));
        v1 += calcScore1(ct, nb);
        v2 += calcScore2(ct, nb);
      }
    }
    return score(v1, v2);
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
  }

  void replace(int p1, int p2) {
    int h1, w1, h2, w2;
    to(p1, h1, w1);
    to(p2, h2, w2);
    {  // remove
      double rexp = 0.5 + remain * 0.5;
      for (int i = h1; i < h2; ++i) {
        for (int j = w1; j < w2; ++j) {
          int p = to(i, j);
          int t = board[p];
          if (t != BOARD[p] && get_random_double() < rexp) putItem(p, 0);
        }
      }
      for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
          int p = to(i, j);
          int t = board[p];
          if (isL(t)) {
            for (int k = 0; k < 4; ++k) {
              int n = light[p][k];
              if (n != -1 && isL(board[n])) putItem(n, 0);
            }
          }
        }
      }
    }
    static int edge[M * 4];
    static int pl[M];
    int es = 0, ps = 0;
    for (int i = h1; i <= h2; ++i) {
      for (int j = w1; j <= w2; ++j) {
        int p = to(i, j);
        if (BOARD[p] != 0) continue;
        pl[ps++] = p;
        constexpr int X[] = {0, 9, 10, 12, 16, 17, 18};
        for (int t : X) {
          edge[es++] = (p << 8) | t;
        }
      }
    }
    assert(es < M * 4);
    while (es > 0) {
      double _score = invalid;
      int p = 0, t = 0;
      if (get_random() % 20 == 0) {
        int z = ps;
        while (z > 0) {
          int i = get_random() % z;
          p = pl[i];
          t = board[p];
          if (isL(t)) break;
          swap(pl[i], pl[--z]);
        }
        if (!isL(t)) continue;
        _score = score();
        putItem(p, 0);
        z = p;
        for (int d = 0; d < 4; ++d) {
          int pp = light[z][d];
          if (pp == -1) continue;
          int tt = board[pp];
          if (!isC(tt)) continue;
          if (tt != (tt | (t ^ 8))) continue;
          auto calc = [&](int a, int b) {
            while (true) {
              a += DIR[b];
              if (not in(a)) break;
              if (a == z) break;
              if (BOARD[a] == 0 && isValid(a, t)) {
                double x = diffScore(a, t) + get_random_double();
                if (_score < x) {
                  _score = x;
                  p = a;
                }
              }
              if (board[a] == 0) continue;
              if (board[a] == 17) {
                b = 3 - b;
              } else if (board[a] == 18) {
                b = b ^ 1;
              } else {
                break;
              }
            }
          };
          calc(z, d);
          calc(z, rev(d));
        }
      } else {
        int i = get_random() % es;
        p = edge[i] >> 8;
        t = edge[i] & 0xff;
        edge[i] = edge[--es];
        if (!isValid(p, t)) continue;
        _score = diffScore(p, t);
      }
      if (_score - score() > remain * log(get_random_double())) {
        putItem(p, t);
        // assert(_score == score());
      }
    }
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

  bool isValid(int p, int t) {
    // if (BOARD[p] != 0) return false;
    int* l = light[p];
    if (t == 0) {
      for (int i = 0; i < 2; ++i) {
        if (isL(board[l[i]]) && isL(board[l[i + 2]])) return false;
      }
    } else if (isL(t)) {
      for (int i = 0; i < 4; ++i) {
        if (isL(board[l[i]])) return false;
      }
    } else if (isO(t)) {
      if (obstacles == MO) return false;
    } else if (isM(t)) {
      if (mirrors == MM) return false;
      for (int i = 0; i < 4; i += 2) {
        if (t == 17) {
          if (isL(board[l[i]]) && isL(board[l[i ^ 1]])) return false;
        } else {
          if (isL(board[l[i]]) && isL(board[l[3 - i]])) return false;
        }
      }
    }
    return true;
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
      cur.init();
      while (true) {
        remain = 1.0 - timer.getElapsed() / TIME_LIMIT;
        if (remain < 0) break;
        memcpy(&tmp, &cur, sizeof(cur));
        constexpr int MASK = 6;
        int h = get_random() % (H - MASK);
        int w = get_random() % (W - MASK);
        tmp.replace(to(h, w), to(h + MASK, w + MASK));
        if (tmp.score() - cur.score() > 5 * remain * log(get_random_double())) {
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
