#include "reed-sloane.cpp"
#include <iostream>
#include <vector>

using namespace std;

// https://www.spoj.com/status/FINDLR
void solve() {
    ios_base::sync_with_stdio(0), cin.tie(0);
    int T;
    cin >> T;
    while (T--) {
        int n, mod;
        cin >> n >> mod;
        vector<int> S(n * 2);
        for (int i = 0; i < n * 2; i++) cin >> S[i];
        if (mod == 1) {
            cout << 0 << '\n';
            continue;
        }

        auto A = ReedSloane(S, mod);

        const int64_t M2 = static_cast<int64_t>(mod) * mod;
        int64_t ans = 0;
        for (size_t i = 0; i < A.size(); i++)
            modadd(ans, static_cast<int64_t>(S[n * 2 - i - 1]) * A[i], M2);
        ans %= mod;
        cout << ans << '\n';
    }
}

void show(vector<int> seq, int mod, size_t n) {
    assert(!seq.empty());
    assert(n >= seq.size());
    // show the first n terms
    auto rel = ReedSloane(seq, mod);
    Mint::setMod(mod);
    vector<Mint> result(n);
    for (size_t i = 0; i < seq.size(); i++) result[i] = seq[i];
    for (size_t i = seq.size(); i < n; i++) {
        for (size_t j = 0; j < rel.size(); j++) {
            result[i] += result[i-j-1] * rel[j];
        }
    }

    cout << "\033[1;32m";
    cout << "mod = " << mod << '\n';
    cout << "S = [";
    for (size_t i = 0; i < seq.size(); i++)
        cout << result[i] << (i+1==n ? "" : ", ");
    cout << "\033[91m";
    for (size_t i = seq.size(); i < n; i++)
        cout << result[i] << (i+1==n ? "" : ", ");
    cout << "]\n";
    cout << "\033[1;32m";
    if (rel.empty()) {
        cout << "a[n] = 0 for big enough n\n";
    } else {
        cout << "deduced linear relation: ";
        cout << "a[n] = ";
        bool first = true;
        for (size_t i = 0; i < rel.size(); i++) {
            int r = rel[i];
            if (r == 0) continue;
            if (r <= mod / 2) {
                cout << (first ? "" : "+ ") << r << ' ';
            } else {
                cout << "- " << mod-r << ' ';
            }
            first = false;
            cout << "a[n-" << i+1 << "] ";
        }
        cout << '\n';
    }
    cout << "\033[0m";
    cout << '\n';
}

int deduce(vector<int> seq, int mod, int64_t n) {
    // deduce only the n-th term
    auto A = ReedSloane(seq, mod);
    Mint::setMod(mod);
    vector<Mint> rec(A.begin(), A.end());
    vector<Mint> init(seq.begin(), seq.end());
    return static_cast<int>(linearReccurenceKthTerm<Mint>(rec, init, n));
}

signed main() {
    // solve();

    vector<int> mods = {2, 100, 720720, 998244353};
    for (int mod: mods) {
        show({1, 1, 1, 1}, mod, 10);

        show({1, 2, 4, 8}, mod, 10);

        show({0, 1, 0, 1, 0, 1}, mod, 10);

        show({1, 3, 6, 10, 15, 21}, mod, 10);

        show({10000, 100000, 1000000, 10000000}, mod, 10);

        show({0, 2, 16, 72, 256, 800}, mod, 10); // a[n] = n^2 * 2^n

        show({41, 41, 43, 47, 53, 61}, mod, 10); // a[n] = n^2 - n + 41
    }
}
