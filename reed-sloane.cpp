#include <vector>
#include <iostream>
#include <utility>
#include <functional>
#include <cmath>
#include <cassert>

using std::vector;
using std::istream, std::ostream;
using std::tuple, std::pair, std::tie;
using std::function;
using std::max, std::swap;

template <typename T>
class Modular {
public:
    constexpr Modular() : v() {}
    template <typename U> Modular(const U &u) { v = (0 <= u && u < MOD ? u : (u%MOD+MOD)%MOD); }
    template <typename U> explicit operator U() const { return U(v); }
    Modular &operator +=(const Modular &rhs) {
        return v += rhs.v - MOD, v += (v >> width) & MOD, *this;
    }
    Modular &operator -=(const Modular &rhs) {
        return v -= rhs.v, v += (v >> width) & MOD, *this;
    }
    // fits for MOD^2 <= ULLONG_MAX
    Modular &operator *=(const Modular &rhs) {
        return v = fm.reduce(1ULL * v * rhs.v), *this;
    }
    Modular &operator /=(const Modular &rhs) {
        return *this *= inverse(rhs.v);
    }
    friend Modular operator +(Modular a, const Modular &b) { return a += b; }
    friend Modular operator -(Modular a, const Modular &b) { return a -= b; }
    friend Modular operator *(Modular a, const Modular &b) { return a *= b; }
    friend Modular operator /(Modular a, const Modular &b) { return a /= b; }
    Modular operator-() const { return 0 - *this; }
    friend bool operator == (const Modular &lhs, const Modular &rhs) { return lhs.v == rhs.v; }
    friend bool operator != (const Modular &lhs, const Modular &rhs) { return lhs.v != rhs.v; }
    friend istream & operator>>(istream &I, Modular &m) { T x; I >> x, m = x; return I; }
    friend ostream & operator<<(ostream &O, const Modular &m) { return O << m.v; }
    static void setMod(T x) {
        MOD = x;
        fm.init(x);
    }
private:
    constexpr static int width = sizeof(T) * 8 - 1;
    T v;
    static T inverse(T a) {
        // copy from tourist's template
        T u = 0, v = 1, m = MOD;
        while (a != 0) {
            T t = m / a;
            m -= t * a; swap(a, m);
            u -= t * v; swap(u, v);
        }
        assert(m == 1);
        return u;
    }
    static T MOD;
    // Barret reduction
    // https://gist.github.com/simonlindholm/51f88e9626408723cf906c6debd3814b
    struct FastMod {
        typedef unsigned long long ull;
        typedef __uint128_t L;
	ull b, m;
        void init(ull _b) {
            b = _b;
            m = ull((L(1) << 64) / b);
        }
	ull reduce(ull a) {
		ull q = (ull)((L(m) * a) >> 64), r = a - q * b;
		return r >= b ? r - b : r;
	}
    };
    static FastMod fm;
};
template <typename T> T Modular<T>::MOD;
template <typename T> typename Modular<T>::FastMod Modular<T>::fm;

using Mint = Modular<int>;

template <typename T>
struct Poly : vector<T> {
    Poly (size_t sz = 0, T val = 0) : vector<T>(sz, val) {}
    Poly & operator-=(const Poly &rhs) {
        if (this -> size() < rhs.size()) this -> resize(rhs.size());
        for (size_t i = 0; i < rhs.size(); i++)
            (*this)[i] -= rhs[i];
        return *this;
    }
    friend Poly operator*(T x, Poly p) {
        for (T &v: p) v *= x;
        return p;
    }
    Poly operator<<(int x) {
        Poly res = *this;
        res.insert(res.begin(), x, 0);
        return res;
    }
    T at(size_t idx) const {
        if (idx >= this -> size()) return 0;
        return (*this)[idx];
    }
};

struct Factorization {
    constexpr static int maxn = 100000;
    vector<int> primes;
    bool sieve[maxn];
    Factorization() {
        for (int i = 2; i < maxn; i++) {
            if (!sieve[i]) {
                primes.push_back(i);
                for (int j = i * 2; j < maxn; j += i) {
                    sieve[j] = true;
                }
            }
        }
    }
    vector<tuple<int,int,int>> operator()(int n) const {
        vector<tuple<int,int,int>> result;
        int s = sqrt(n) + 1;
        for (int i: primes) {
            if (i > s) break;
            if (n % i == 0) {
                int c = 0, p = 1;
                while (n % i == 0) n /= i, ++c, p *= i;
                result.emplace_back(i, c, p);
            }
        }
        if (n > 1)
            result.emplace_back(n, 1, n);
        return result;
    }
} factorize;

template <typename T>
pair<T,T> extgcd(T a, T b) {
    if (!b) return { 1, 0 };
    auto [x, y] = extgcd(b, a % b);
    return {y, x-a/b*y};
}
function<int(int,int)> CRT(int m1, int m2) { // assume gcd(m1, m2) == 1
    auto [x, y] = extgcd<int>(m1, m2);
    int M = m1 * m2;
    int64_t c2 = static_cast<int64_t>(x) * m1 % M;
    int64_t c1 = static_cast<int64_t>(y) * m2 % M;
    function<int(int,int)> crt = [M, c1, c2](int r1, int r2) -> int {
        int res = (r1 * c1 + r2 * c2) % M;
        return res < 0 ? res + M : res;
    };
    return crt;
}

void modadd(int64_t &x, int64_t v, int64_t m) {
    // 0 <= x, v < m
    // x := (x + v) % m
    x += v - m;
    x += (x >> 63) & m;
}

vector<int> ReedSloane(const vector<int> &S, int mod) {
    assert (!S.empty());

    const int64_t M2 = static_cast<int64_t>(mod) * mod;
    auto solvePrimePower = [&S, &M2](int p, int e) -> vector<int> {
        vector<Poly<Mint>> a(e), a_new(e), a_old(e);
        vector<int> L(e), L_new(e), L_old(e, -1);
        vector<int> r(e);
        vector<Mint> theta(e), theta_old(e);
        vector<int> u(e), u_old(e);
        vector<int> power(e); // power[i] = p^i
        power[0] = 1;
        for (int i = 1; i < e; i++) power[i] = power[i-1] * p;

        auto invertible = [p, e](Mint v) -> pair<Mint,int> {
            // returns the invertible part and the prime power part.
            int x = static_cast<int>(v);
            if (!x) return { 1, e };
            int c = 0;
            while (x % p == 0) x /= p, ++c;
            return { x, c };
        };

        auto getError = [&S, &a, &M2](int r, int k) -> Mint {
            // Instead of "err = (err + S[i] * a[r].at(k-i)()) % mod",
            // calculate the sum under modulo M2 = mod^2 would be
            // faster using modadd function.
            int64_t err = 0;
            for (int i = 0; i <= k; i++) {
                modadd(err, S[i] * static_cast<int64_t>(a[r].at(k-i)), M2);
            }
            return err;
        };

        for (int i = 0; i < e; i++) {
            a[i] = Poly<Mint>(1, power[i]);
            L[i] = 0;
            a_new[i] = Poly<Mint>(1, power[i]);
            tie(theta[i], u[i]) = invertible(getError(i, 0));
            if (u[i] == e) {
                L_new[i] = 0;
            } else {
                L_new[i] = 1;
            }
        }
        for (int k = 1; k < S.size(); k++) {
            for (int g = 0; g < e; g++) {
                if (L_new[g] > L[g]) {
                    int h = e - 1 - u[g];
                    a_old[g] = a[h];
                    L_old[g] = L[h];
                    theta_old[g] = theta[h];
                    u_old[g] = u[h];
                    r[g] = k-1;
                }
            }

            a = a_new;
            L = L_new;

            for (int i = 0; i < e; i++) {
                tie(theta[i], u[i]) = invertible(getError(i, k));
                if (u[i] == e) continue;
                int g = e - 1 - u[i];
                if (L[g] == 0) {
                    L_new[i] = max(L_new[i], k + 1);
                } else {
                    Mint coef = (theta[i] / theta_old[g]) * power[u[i] - u_old[g]];
                    int shift = k - r[g];
                    a_new[i] -= coef * (a_old[g] << shift);
                    L_new[i] = max(L_new[i], L_old[g] + shift);
                }
            }
        }
        vector<int> A(a_new[0].begin(), a_new[0].end());
        A.resize(L_new[0] + 1);
        return A;
    };

    vector<int> A;
    int M = 1;
    for (auto [p, e, m]: factorize(mod)) {
        Mint::setMod(m);
        auto a = solvePrimePower(p, e);

        if (A.size() < a.size()) A.resize(a.size());
        auto crt = CRT(M, m);
        M *= m;
        for (size_t i = 0; i < A.size(); i++)
            A[i] = crt(A[i], (i < a.size() ? a[i] : 0));
    }
    A.erase(A.begin());
    for (int &x: A)
        x = (x ? M - x : 0);
    return A;
}

template <typename T> T linearReccurenceKthTerm(vector<T> rec, vector<T> init, uint64_t k) {
    // NOTE: not tested
    if (rec.empty())
        return T(0);
    if (k < init.size())
        return init[k];

    assert( init.size() >= rec.size() );

    vector<T> r{1}, e{0, 1};
    auto mul = [&rec](const vector<T> &a, const vector<T> &b) -> vector<T> {
        vector<T> c(a.size() + b.size() - 1);
        for (size_t i = 0; i < a.size(); i++)
            for (size_t j = 0; j < b.size(); j++)
                c[i+j] += a[i] * b[j];
        for (size_t i = c.size()-1; i >= rec.size(); i--)
            for (size_t j = 0; j < rec.size(); j++)
                c[i-j-1] += c[i] * rec[j];
        c.resize(rec.size());
        return c;
    };
    while (k) {
        if (k & 1)
            r = mul(r, e);
        e = mul(e, e);
        k >>= 1;
    }
    T res = 0;
    for (size_t i = 0; i < r.size() && i < init.size(); i++)
        res += r[i] * init[i];
    return res;
}
