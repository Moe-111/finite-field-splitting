#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace std;

// ------------------------------ 基础数论工具 ------------------------------
// 素性检测
bool is_prime(int n) {
    if (n < 2) return false;
    for (int i = 2; i * i <= n; ++i)
        if (n % i == 0) return false;
    return true;
}

// 分解 q = p^k，返回素数p，k存入exp
int prime_power(int q, int& exp) {
    exp = 1;
    if (q < 2) return -1;
    for (int p = 2; p <= q; ++p) {
        if (!is_prime(p)) continue;
        int tmp = q, cnt = 0;
        while (tmp % p == 0) tmp /= p, cnt++;
        if (tmp == 1) {
            exp = cnt;
            return p;
        }
    }
    return -1;
}

int gcd(int a, int b) { return b == 0 ? a : gcd(b, a % b); }
int lcm(int a, int b) { return a / gcd(a, b) * b; }

// ------------------------------ 有限域F_p上多项式运算 ------------------------------
using Poly = vector<int>; // 系数[低次→高次], 模p

// 多项式次数
int deg(const Poly& f) {
    int d = f.size() - 1;
    while (d > 0 && f[d] == 0) d--;
    return d;
}

// 模p化简
void mod_p(Poly& f, int p) {
    for (int& x : f) x = (x % p + p) % p;
}

// 多项式乘法(模p)
Poly mul(const Poly& f, const Poly& g, int p) {
    Poly res(deg(f) + deg(g) + 1, 0);
    for (int i = 0; i < f.size(); ++i)
        for (int j = 0; j < g.size(); ++j)
            res[i+j] = (res[i+j] + f[i] * g[j]) % p;
    mod_p(res, p);
    return res;
}

// ------------------------------ 有限域多项式不可约判定 ------------------------------
bool is_irreducible(const Poly& f, int p) {
    int n = deg(f);
    if (n == 0) return false;
    if (n == 1) return true;
    // 试除所有次数≤n/2的首一多项式(简化教学版)
    for (int d = 1; d <= n/2; ++d) {
        // 遍历d次首一多项式
        int max_c = 1;
        for (int i = 0; i < d; ++i) max_c *= p;
        for (int c = 0; c < max_c; ++c) {
            Poly g(d+1, 0); g[d] = 1;
            int tmp = c;
            for (int i = 0; i < d; ++i) g[i] = tmp % p, tmp /= p;
            mod_p(g, p);
            if (deg(g) != d) continue;
            // 除法试除(简化：乘积匹配)
            for (int k = 1; k <= n - d + 1; ++k) {
                Poly h(k, 0); h[k-1] = 1;
                if (mul(g, h, p) == f) return false;
            }
        }
    }
    return true;
}

// ------------------------------ 分解多项式，收集不可约因子次数 ------------------------------
vector<int> get_irreducible_degrees(Poly f, int p) {
    vector<int> res;
    mod_p(f, p);
    int n = deg(f);
    if (n == 0) return res;

    for (int d = 1; d <= n; ++d) {
        for (int c = 1; c < p; ++c) {
            Poly fac(d+1, 0); fac[d] = c;
            for (int i = 0; i < d; ++i) fac[i] = rand() % p;
            mod_p(fac, p);
            if (deg(fac) != d || !is_irreducible(fac, p)) continue;
            // 简单因子匹配
            while (true) {
                bool found = false;
                for (int k = 1; k <= deg(f) - d + 1; ++k) {
                    Poly h(k, 0); h[k-1] = 1;
                    if (mul(fac, h, p) == f) {
                        res.push_back(d);
                        f = h;
                        found = true;
                        break;
                    }
                }
                if (!found) break;
            }
        }
    }
    if (deg(f) > 0) res.push_back(deg(f));
    return res;
}

// ------------------------------ 主函数：输入→计算→输出分裂域与维数 ------------------------------
int main() {
    cout << "===== 有限域多项式分裂域计算器 (抽象代数专用) =====" << endl;
    int q;
    cout << "请输入有限域阶 q (必须是素数的幂 p^k): ";
    cin >> q;

    int k;
    int p = prime_power(q, k);
    if (p == -1) {
        cout << "错误：" << q << " 不是素数的幂！" << endl;
        return 1;
    }
    cout << "→ 基域 F_" << q << " = F_" << p << "^" << k << endl;

    // 输入多项式：系数从常数项到最高次
    cout << "\n请输入多项式系数(常数项开始，以-1结束): ";
    Poly f;
    int x;
    while (cin >> x && x != -1) f.push_back(x);
    mod_p(f, p);
    cout << "→ 多项式(模" << p << ")：次数 = " << deg(f) << endl;

    // 不可约因子次数
    vector<int> degs = get_irreducible_degrees(f, p);
    cout << "\n不可约因子次数: ";
    for (int d : degs) cout << d << " ";
    cout << endl;

    // 计算lcm = 分裂维数
    int split_dim = 1;
    for (int d : degs) split_dim = lcm(split_dim, d);
    int split_q = pow(q, split_dim);

    // 最终输出
    cout << "\n===== 结果 =====" << endl;
    cout << "分裂域 L = F_" << q << "^" << split_dim << " = F_" << split_q << endl;
    cout << "分裂维数 [L : F_" << q << "] = " << split_dim << endl;
    cout << "===========================================" << endl;
    return 0;
}
