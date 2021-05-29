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

void show(vector<int> seq, int mod, int n) {
    // show the first n terms
    ;
}

int deduce(vector<int> seq, int mod, int n) {
    // deduce only the n-th term
    ;
}

signed main() {
    // solve();

    /* WIP
    vector<int> mods = {998244353, 2, 77777, 14, 49, 48, 720720};
    for (int mod: mods) {
        show({1, 1, 1, 1}, mod, 10);

        show({1, 2, 4, 8}, mod, 10);

        show({0, 1, 0, 1, 0, 1}, mod, 10);

        show({1, 3, 6 }, mod, 10);
    }
    */
}
