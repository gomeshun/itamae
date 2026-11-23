# Physlib → Certified H(z): 最小試作

## 状態（2026-11-23）

- **実行確認済み**: Python標準ライブラリによる証明書生成、独立した厳密有理数による証明書検査、8件の単体テスト。
- **ソース確認済み**: Physlibの`reducedHubble`定義と`FriedmannSolution.hubbleConstant_eq_reducedHubble`定理（以下のコミット）。
- **Leanソース生成済み／未コンパイル**: `src/HzCertificate.lean`と`generated/HzExample.lean`。この制作環境にはLean/Lakeがなく、外部ダウンロードも制限されているため、**Leanのカーネルによる証明検査は未実施**。コンパイルの結果を偽って報告していない。実際のLean検査を完了するには、下記手順に従って`lake build`を実行する。

## 想定する物理モデル

平坦な物質＋宇宙定数（放射なし）のΛCDM：

\[
H(z) = H_0 E\left(\frac{1}{1+z}\right),\qquad
E(x) = \sqrt{(1-\Omega_m) + \Omega_m x^{-3}},
\]

すなわち

\[
H(z)=H_0\sqrt{(1-\Omega_m)+\Omega_m(1+z)^3}.
\]

制約：`H0 > 0`, `0 ≤ Omega_m ≤ 1`, `z ≥ 0`。`H0`と`H(z)`の値の単位は`km s^-1 Mpc^-1`、赤方偏移と密度パラメータは無次元。観測誤差とパラメータ事後分布は今回の範囲外。

Physlib（固定コミット`db96e88549c3ed2b1211fe34cd751f30847e5386`）の
`Cosmology.FLRW.FriedmannEquation.reducedHubble`を`modelH`で**直接呼び出す**。
Physlibの理論定理`FriedmannSolution.hubbleConstant_eq_reducedHubble`を元にしたモデルだが、今回の数値証明は`reducedHubble`の具体的なパラメータ指定に対するものであり、特定のFLRW時空を構成してその仮定を証明したものではない。

## 例（H0=70, Omega_m=3/10, z=1）

厳密な計算：

\[
H(1)^2 = 70^2\left(\frac{3}{10}\,2^3+\frac{7}{10}\right)
=15190.
\]

`--digits 8` で得られる上下限（いずれも**正確な有理数**）：

\[
\boxed{123.24771803 \leq H(1) \leq 123.24771804}
\]

検証器が正確な有理数演算だけで判定する条件：

\[
(123.24771803)^2 \leq 15190 \leq (123.24771804)^2.
\]

誤差区間幅は\(10^{-8}\)以下。参考表示値は約 `123.24771803161306`。

## 再現方法

### (A) Pythonだけで証明書を生成し、厳密な有理数計算で検査

Python 3.10+、追加パッケージ不要：

```bash
python3 scripts/hz_cert.py \
  --h0 70 --omega-m 3/10 --z 1 --digits 8 \
  --certificate generated/hz_z1.json --lean generated/HzExample.lean

python3 scripts/hz_cert.py --verify generated/hz_z1.json
python3 -m unittest discover -s tests -v
```

別の入力も指定可能：

```bash
python3 scripts/hz_cert.py --h0 67.4 --omega-m 0.315 --z 2 --digits 10 \
  --certificate generated/other.json --lean generated/HzExample.lean
```

この操作で生成されたLeanファイルは、**入力値の変更後に別途Leanで検査**する。

### (B) Physlibに依存するLeanの証明を検査

Lean 4/Lake/Elanをインストールし、ネットワーク接続がある環境で、このディレクトリから：

```bash
lake update
lake exe cache get    # Mathlibなどのキャッシュが利用できる場合
lake build
```

`lakefile.toml`でPhyslibのコミットを固定し、`lean-toolchain`でLean `v4.36.0`を使用する。Physlibの依存先Mathlibも同じ系列を使用する。`src/HzCertificate.lean`は一般定理`modelH_bounds`を、`generated/HzExample.lean`はデータごとのインスタンス`certified_interval`と`certified_width`を宣言する。最後に`#print axioms`で定理の公理依存を出力する。

**重要**：この制作環境では`lake build`を実行できていないため、上記ファイルはLeanコードの*候補実装*である。型エラーなどが発生する場合はコンパイル出力を用いて修正する。Pythonによる有理数検査の成功を、Leanの形式検証完了と同一視しない。

## 仕組み

1. 物理理論側：Physlibの`reducedHubble`を`ℝ`上の正確な関数として利用。
2. 変形の数学的証明：`x=1/(1+z)`から`x^(-3)=(1+z)^3`を導出。
3. データ専用の計算：Pythonの`fractions.Fraction`と`math.isqrt`で整数演算のみを使い、正確な上下限を生成。
4. 証明書：`H0, Omega_m, z, L, U, H(z)^2`をJSON保存。
5. Python検証器：入力から`H(z)^2`を再計算し`L^2 ≤ H(z)^2 ≤ U^2`と区間幅を検査。改変した証明書は拒否する。
6. Lean定理案：条件が成立すれば`L ≤ modelH H0 Omega_m z ≤ U`とする一般定理を証明し、各インスタンスの前提を`norm_num`で処理。

## 信頼境界と拡張

- Pythonの整数演算・分数演算に基づく検算は実行確認済み。`math.isqrt`で候補境界を作るが、検証器は境界から**独立に平方不等式を再計算**する。
- Leanカーネルで証明を確認したとき、Pythonの生成器を信頼する必要はない。生成されたLeanの全前提を確認する。
- `Float`からの計算誤差を評価する実験ではなく、正確な入力有理数で定義された物理予測を有理数区間に包囲する試作である。
- 実験装置の較正・測定誤差やモデルの物理的適用性は今回の保証対象外。
- 次は`H(z)`の保証付き評価器を再利用し、距離積分`∫ dz/H(z)`の有理数区間証明書と、ガウス観測モデルの対数尤度へ接続する。

## 参照コード

- Physlib: https://github.com/leanprover-community/physlib
- 対応する定義: https://github.com/leanprover-community/physlib/blob/db96e88549c3ed2b1211fe34cd751f30847e5386/Physlib/Cosmology/FLRW/ReducedHubble.lean
- Physlibのツールチェーン: `leanprover/lean4:v4.36.0`
