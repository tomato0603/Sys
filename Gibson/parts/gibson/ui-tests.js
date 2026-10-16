/* Tests for the text that the request form hands to Claude. The expected
 * strings are written out in full, so any change to the format shows up. */
(function () {
  'use strict';
  var test = PBTESTS.test, eq = PBTESTS.eq, ok = PBTESTS.ok;

  test('リクエスト: 日時は YYYY-MM-DD HH:MM（1 桁は 0 で埋める）', function () {
    eq(GIBSONUI.stamp(new Date(2026, 0, 5, 9, 3)), '2026-01-05 09:03');
    eq(GIBSONUI.stamp(new Date(2026, 11, 31, 23, 59)), '2026-12-31 23:59');
  });

  test('リクエスト: 前後の空白・空行は取り、中の改行は残す', function () {
    eq(GIBSONUI.tidy('\n  1 行目\r\n2 行目  \n\n'), '1 行目\n2 行目');
    eq(GIBSONUI.tidy('   '), '');
  });

  test('リクエスト: 1 件の書式', function () {
    eq(GIBSONUI.formatRequests(
         [{ at: '2026-09-28 14:03', text: 'Tm の列を右寄せにしてほしい' }],
         '2026-09-28 12:07'),
       '【Gibson 修正リクエスト】全 1 件\n'
       + 'ツール: Gibson.html（作製 2026-09-28 12:07）\n'
       + '\n'
       + '[1] 2026-09-28 14:03\n'
       + 'Tm の列を右寄せにしてほしい');
  });

  test('リクエスト: 何件でも番号順につなげる', function () {
    eq(GIBSONUI.formatRequests([
         { at: '2026-09-28 14:03', text: 'A' },
         { at: '2026-09-28 14:10', text: 'B\nC' }
       ], 'X'),
       '【Gibson 修正リクエスト】全 2 件\nツール: Gibson.html（作製 X）\n'
       + '\n[1] 2026-09-28 14:03\nA\n'
       + '\n[2] 2026-09-28 14:10\nB\nC');
  });

  // ---- ナビの自動切り替え --------------------------------------------

  test('ナビ: ヘッダーの下を通り過ぎた最後の区画を選ぶ', function () {
    // 区画の上端（画面の上からの位置）。ヘッダーの下端は 100
    eq(GIBSONUI.sectionInView([120, 900, 1500], 100, false), 0);   // まだ最初の区画
    eq(GIBSONUI.sectionInView([-800, 60, 700], 100, false), 1);
    eq(GIBSONUI.sectionInView([-2000, -900, -10], 100, false), 2);
  });

  test('ナビ: 一番下までスクロールしたら最後の区画を選ぶ', function () {
    // 最後の区画が短く、上端がヘッダーまで届かなくても選ぶ
    eq(GIBSONUI.sectionInView([-2000, -300, 400], 100, true), 2);
    eq(GIBSONUI.sectionInView([], 100, false), -1);
  });

  // ---- 条件の保存 ----------------------------------------------------

  test('保存: 一覧の説明は方式・本数・入力値・Tm の配列数を書く', function () {
    eq(GIBSONUI.conditionSummary({ mode: 'manual', tmText: 'ATGC\n\nGGCC\n',
      rows: [{ concentration: '166', bp: '5238', ratio: '', ul: '0.35' },
             { concentration: '100', bp: '308', ratio: '', ul: '0.8' }] }),
      'Volume を入力 / Insert 1 本 / Vector 166 ng/µL・5238 bp・0.35 µL / '
      + 'Insert 100 ng/µL・308 bp・0.8 µL / Tm の配列 2 本');
    eq(GIBSONUI.conditionSummary({ mode: 'ratio', tmText: '',
      rows: [{ concentration: '166', bp: '5238', ratio: '', ul: '0.35' },
             { concentration: '', bp: '', ratio: '', ul: '' },
             { name: 'GFP', concentration: '100', bp: '1500', ratio: '2.5', ul: '' }] }),
      '倍数から計算 / Insert 2 本 / '
      + 'Vector 166 ng/µL・5238 bp・0.35 µL / Insert A （空欄） / GFP 100 ng/µL・1500 bp・×2.5');
  });

  // ---- 結果の出力 ----------------------------------------------------

  // The bench sheet of 2026-10-10, as GIBSON.reaction returns it.
  var REACTION_REPORT = { tm: null, reaction: {
    rows: [{ name: 'Vector', concentration: 166, bp: 5238, ul: 0.35,
             pmolPerUl: 0.0488, pmol: 0.0171, ng: 58.1 },
           { name: 'Insert', concentration: 100, bp: 308, ul: 0.8,
             pmolPerUl: 0.4995, pmol: 0.3996, ng: 80, ratio: 1.38 }],
    mixUl: 10, waterUl: 8.85, dnaUl: 1.15, mixFold: 2,
    totals: { pmol: 0.4167, ul: 20 },
    checks: [{ label: 'Vector の量', target: '50–100 ng', value: '58.1 ng', ok: true },
             { label: '合計 pmol（3 fragment）', target: '目安なし',
               value: '0.4167 pmol', ok: null },
             { label: 'ピペットで取れる量（最小の Volume）', target: '0.22 µL 以上',
               value: '0.2 µL（Vector）', ok: false }],
    notes: ['Vector は 0.35 µL です。'] } };

  test('出力: 反応の組み立ては実験ノートの表と同じ並び', function () {
    var md = GIBSONUI.formatReport(REACTION_REPORT, 'B', 'T');
    eq(md.split('\n').slice(0, 15), [
      '# Gibson 計算結果', '',
      '- 出力: T',
      '- ツール: Gibson.html（作製 B）', '',
      '## 反応の組み立て', '',
      '|   | Conc ng/µL | Length bp | pmol/µL | target pmols | Volume µL | 確認用 ng |',
      '| --- | --- | --- | --- | --- | --- | --- |',
      '| Vector | 166 | 5238 | 0.0488 | 0.0171 | 0.35 | 58.1 |',
      '| Insert | 100 | 308 | 0.4995 | 0.3996 | 0.8 | 80 |',
      '| G.A. M ×2 | — | — | — | — | 10 | — |',
      '| miliQ | — | — | — | — | 8.85 | — |',
      '| Total | — | — | — | 0.4167 | 20 | — |', '']);
  });

  test('出力: 条件設定は OK / NG / 判定なし で書き、注意書きは出さない', function () {
    var md = GIBSONUI.formatReport(REACTION_REPORT, 'B', 'T');
    ok(md.indexOf('| Vector の量 | 50–100 ng | 58.1 ng | OK |') !== -1, md);
    ok(md.indexOf('| 合計 pmol（3 fragment） | 目安なし | 0.4167 pmol | 判定なし |') !== -1, md);
    ok(md.indexOf('| ピペットで取れる量（最小の Volume） | 0.22 µL 以上 | 0.2 µL（Vector） | NG |') !== -1, md);
    ok(md.indexOf('0.35 µL です') === -1, md);
  });

  test('出力: 倍数方式では計算確認の表も出す', function () {
    var report = JSON.parse(JSON.stringify(REACTION_REPORT));
    report.reaction.vectorNg = 58.1;
    report.reaction.targets = [
      { name: 'Vector', targetNg: 58.1, ul: 0.35, actualNg: 58.1, targetRatio: null, actualRatio: null },
      { name: 'Insert', targetNg: 145.3, ul: 1.45, actualNg: 145, targetRatio: 2.5, actualRatio: 2.5 }];
    var md = GIBSONUI.formatReport(report, 'B', 'T');
    ok(md.indexOf('### 計算確認（Vector 58.1 ng と倍数から計算）') !== -1, md);
    ok(md.indexOf('| Insert | 145.3 | 1.45 | 145 | ×2.5 | ×2.5 |') !== -1, md);
    ok(GIBSONUI.tablesOnly(md).indexOf('| Insert | 145.3 |') !== -1, '印刷にも出る');
  });

  test('出力: 計算していない区画と、その計算式は出さない', function () {
    var md = GIBSONUI.formatReport(REACTION_REPORT, 'B', 'T');
    ok(md.indexOf('## Tm値の計算') === -1 && md.indexOf('簡易式') === -1, md);
    ok(md.indexOf('全量 ÷ 2') !== -1, md);
  });

  test('出力: Tm は表と Tm 差、適用範囲を書く', function () {
    var md = GIBSONUI.formatReport({ reaction: null, tm: {
      rows: [{ sequence: 'ATGC', length: 4, tm: 7, gcPercent: 50 },
             { sequence: 'ATGCATGC', length: 8, tm: 19, gcPercent: 50 }],
      diff: { gap: 12, a: 1, b: 2, spread: false } } }, 'B', 'T');
    ok(md.indexOf('## Tm値の計算') !== -1, md);
    ok(md.indexOf('| 2 | ATGCATGC | 8 | 19 | 50 |') !== -1, md);
    ok(md.indexOf('PCR') === -1 && md.indexOf('アニーリング') === -1, 'PCR 条件は出さない');
    ok(md.indexOf('Tm 差: 12 ℃（#1 と #2）') !== -1, md);
    ok(md.indexOf('25 mer 以下') !== -1, md);
    ok(md.indexOf('650') === -1, '量の式は出さない');
  });

  test('出力: 名前の | は表を壊さないようにエスケープし、印刷では元に戻す', function () {
    var report = JSON.parse(JSON.stringify(REACTION_REPORT));
    report.reaction.rows[0].name = 'A|B';
    var md = GIBSONUI.formatReport(report, 'B', 'T');
    ok(md.indexOf('| A\\|B | 166 |') !== -1, md);
    ok(GIBSONUI.reportToHtml(md).indexOf('<td>A|B</td><td>166</td>') !== -1);
  });

  test('出力: 印刷は見出しと表だけ（計算式・文・箇条書きは出さない）', function () {
    var printed = GIBSONUI.tablesOnly(GIBSONUI.formatReport(REACTION_REPORT, 'B', 'T'));
    ok(printed.indexOf('## 計算式と前提') === -1 && printed.indexOf('650') === -1, printed);
    ok(printed.indexOf('### 条件設定') !== -1, printed);
    ok(printed.indexOf('| Vector | 166 |') !== -1, printed);
    eq(printed.split('\n').slice(0, 4),
       ['# Gibson 計算結果', '', '- 出力: T', '- ツール: Gibson.html（作製 B）']);
    eq(GIBSONUI.tablesOnly('# A\n\n## S\n\n文\n\n- x\n\n| a |\n| --- |\n| 1 |'),
       '# A\n\n## S\n\n| a |\n| --- |\n| 1 |');
  });

  test('出力: 印刷用 HTML は見出し・箇条書き・表を作り、文字はエスケープする', function () {
    var html = GIBSONUI.reportToHtml(
      '# T\n\n- a<b\n- c\n\n| x | y |\n| --- | --- |\n| 1 | 2 |\n\n文');
    eq(html, '<h1>T</h1>\n<ul><li>a&lt;b</li><li>c</li></ul>\n'
       + '<table><thead><tr><th>x</th><th>y</th></tr></thead>'
       + '<tbody><tr><td>1</td><td>2</td></tr></tbody></table>\n<p>文</p>');
  });
})();
