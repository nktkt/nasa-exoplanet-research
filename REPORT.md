# NASA Exoplanet Archive Research Report

## Snapshot

- 実行日時: 2026-04-10 23:48 JST
- 確認済みサンプル: 2,896 planets
- TOI 候補サンプル: 854 candidates
- データ: NASA Exoplanet Archive `PSCompPars`, `TOI`

## Question

1. 近接小型トランジット惑星の半径谷は、最新サンプルで照射量に応じてどの程度移動するか。
2. 半径谷にホスト星金属量依存の兆候はあるか。
3. その結果を使うと、どの未確認 TOI が最も高い統計学的リターンを持つか。

## Method

- `PSCompPars` から `0.8 < R_p < 4.5 R_earth`, `0.3 < P < 100 d`, `0.5 < S < 3000 S_earth` のトランジット惑星を抽出。
- `st_rad < 1.8 R_sun`, `3200 < T_eff < 7200 K` で主系列寄りのサンプルに制限。
- 半径分布に KDE を当て、super-Earth peak と sub-Neptune peak の間の密度最小点を valley radius として推定。
- TOI については、確認済み惑星分布からの kNN novelty、明るさ、transit depth、半径谷への近さを幾何平均で統合して priority score を作成。

## Main Findings

### 1. 半径谷は照射量とともに有意に上昇する

推定式:

`log10(R_valley / R_earth) = 0.030 + 0.108 * log10(S / S_earth)`

- 傾きの bootstrap 68% 区間: `0.066` to `0.113`
- 最低照射 bin の valley radius: `1.37 R_earth`
- 最高照射 bin の valley radius: `2.16 R_earth`

この信号はかなり強く、現行サンプルでも半径谷が irradiation-dependent であることを再確認できた。

### 2. 金属量依存は示唆的だが、照射量ほどは堅くない

- 低金属量 tertile (`[Fe/H] median = -0.13`): `1.71 R_earth`
- 中央 tertile (`[Fe/H] median = 0.00`): `1.81 R_earth`
- 高金属量 tertile (`[Fe/H] median = +0.1245`): `2.01 R_earth`
- 高金属量 minus 低金属量: base estimate `0.307 R_earth`
- bootstrap median: `0.262 R_earth`
- bootstrap 68% 区間: `0.082` to `0.407 R_earth`

金属量が高いホストで valley radius が上方にずれる傾向は見える。ただし照射量依存より不確実性が大きく、ここは「有望なシグナル」として扱うべきで、断定は避ける。

## Highest-value TOI Targets

### Radius-valley probes

1. TOI 7466.01: score `80.7`, `R=2.08 R_earth`, `P=13.62 d`, `T=9.81`
2. TOI 1783.01: score `76.9`, `R=2.28 R_earth`, `P=1.42 d`, `T=9.27`
3. TOI 426.01: score `72.4`, `R=2.16 R_earth`, `P=1.32 d`, `T=9.60`
4. TOI 494.01: score `71.7`, `R=1.96 R_earth`, `P=1.70 d`, `T=10.09`
5. TOI 119.02: score `71.2`, `R=1.81 R_earth`, `P=10.69 d`, `T=9.28`

### Temperate / cooler sub-Neptune priorities

1. TOI 6678.01: score `72.5`, `R=2.13 R_earth`, `P=3.43 d`, `S=94`
2. TOI 5112.01: score `72.4`, `R=2.18 R_earth`, `P=15.53 d`, `S=31`
3. TOI 1027.01: score `69.2`, `R=2.88 R_earth`, `P=3.28 d`, `S=161`
4. TOI 5807.01: score `67.7`, `R=2.12 R_earth`, `P=14.24 d`, `S=185`
5. TOI 3353.02: score `67.5`, `R=2.36 R_earth`, `P=8.82 d`, `S=86`

## Interpretation

- いまの NASA サンプルでは、半径谷の主軸はまず irradiation で説明できる。
- 一方で host metallicity も valley を押し上げる方向に効いている可能性があり、これは formation history と envelope retention の両方を含む混合シグナルかもしれない。
- したがって、追観測の優先順位は「谷の近くにある」「明るい」「既存確認サンプルから十分に新規性がある」の3条件で決めるのが合理的。

## Files

- 図: `results/radius_valley_overview.png`
- 図: `results/toi_priority_map.png`
- 候補一覧: `results/top_candidates_overall.csv`
- 半径谷プローブ: `results/top_radius_valley_probes.csv`
- 温和な sub-Neptune: `results/top_temperate_sub_neptunes.csv`
- 数値サマリ: `results/research_summary.json`

## Limits

- selection effects の完全補正はしていない
- metallicity dependence は現段階では tentative
- TOI ranking は follow-up triage 用の empirical score であり、validation probability を直接推定したものではない
