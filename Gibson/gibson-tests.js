
/* Tests for the Gibson calculations. Every expected value is worked out by
 * hand, and the inputs are chosen so the arithmetic comes out even:
 * 0.05 pmol x 5000 bp x 650 / 1000 = 162.5 ng, and 162.5 / 50 = 3.25 uL. */
(function () {
  'use strict';
  var test = PBTESTS.test, eq = PBTESTS.eq, ok = PBTESTS.ok;

  function message(fn) {
    try { fn(); } catch (error) { return error.message; }
    return '';
  }

  // ---- ng / pmol / microlitres -------------------------------------

  test('Gibson: ng から pmol（100 ng の 5000 bp = 0.0308 pmol）', function () {
    eq(GIBSON.significant(GIBSON.pmolFromNg(100, 5000), 3), 0.0308);
  });

  test('Gibson: pmol から ng（0.05 pmol の 5000 bp = 162.5 ng）', function () {
    eq(GIBSON.ngFromPmol(0.05, 5000), 162.5);
  });

  test('Gibson: ng と pmol は往復しても戻る', function () {
    eq(GIBSON.pmolFromNg(GIBSON.ngFromPmol(0.2, 3400), 3400), 0.2);
  });

  test('Gibson: 162.5 ng を 50 ng/µL から取ると 3.25 µL', function () {
    eq(GIBSON.volumeUl(162.5, 50), 3.25);
  });

  test('Gibson: 長さが 0 なら理由を言って断る', function () {
    ok(message(function () { GIBSON.pmolFromNg(100, 0); })
       .indexOf('0 より大きい') !== -1);
  });

  test('Gibson: 長さが整数でなければ断る', function () {
    ok(message(function () { GIBSON.pmolFromNg(100, 1000.5); })
       .indexOf('整数') !== -1);
  });

  test('Gibson: 空欄は「入力してください」と言う', function () {
    ok(message(function () { GIBSON.volumeUl(100, ''); })
       .indexOf('入力してください') !== -1);
  });

  test('Gibson: .8 は 0.8、日本語入力の全角数字・「。」「．」も数として読む', function () {
    eq(GIBSON.normalizeNumberText('.8'), '0.8');
    eq(GIBSON.normalizeNumberText('．３５'), '0.35');
    eq(GIBSON.normalizeNumberText('。8'), '0.8');
    eq(GIBSON.normalizeNumberText('１６６'), '166');
    eq(GIBSON.normalizeNumberText(' 12.5 '), '12.5');
    eq(GIBSON.volumeUl(80, '．8'), 100);
    var rows = GIBSON.reaction({ fragments: [
      { concentration: '１６６', bp: '５２３８', ul: '。35' },
      { concentration: '100', bp: '308', ul: '.8' }] }).rows;
    eq(rows.map(function (r) { return r.ul; }), [0.35, 0.8]);
  });

  // ---- 有効数字 ------------------------------------------------------

  test('Gibson: 有効数字は桁の大小で変わる', function () {
    eq(GIBSON.significant(0.0304999, 3), 0.0305);
    eq(GIBSON.significant(12345, 3), 12300);
    eq(GIBSON.significant(0, 3), 0);
  });

  // ---- 断片の量 ------------------------------------------------------

  var SPEC = { basePmol: 0.05, fragments: [
    { name: 'ベクター', bp: 5000, ratio: 1, concentration: 50 },
    { name: '挿入', bp: 1000, ratio: 2, concentration: 20 }
  ] };

  test('Gibson: モル比のぶんだけ pmol が増える', function () {
    eq(GIBSON.plan(SPEC).rows.map(function (r) { return r.pmol; }),
       [0.05, 0.1]);
  });

  test('Gibson: 各断片の ng と µL を出す', function () {
    var rows = GIBSON.plan(SPEC).rows;
    eq(rows.map(function (r) { return r.ng; }), [162.5, 65]);
    eq(rows.map(function (r) { return r.ul; }), [3.25, 3.25]);
  });

  test('Gibson: 合計は pmol・ng・µL の 3 つとも出す', function () {
    eq(GIBSON.plan(SPEC).totals, { pmol: 0.15, ng: 227.5, ul: 6.5 });
  });

  test('Gibson: 濃度が 1 つでも無ければ µL は出さない', function () {
    var result = GIBSON.plan({ basePmol: 0.05, fragments: [
      { name: 'ベクター', bp: 5000, ratio: 1, concentration: 50 },
      { name: '挿入', bp: 1000, ratio: 2 }
    ] });
    eq(result.rows[1].ul, null);
    eq(result.totals.ul, null);
    // ng までは出せるので、そこまでは答える
    eq(result.totals.ng, 227.5);
  });

  test('Gibson: 名前を書かなければ「断片 n」と呼ぶ', function () {
    eq(GIBSON.plan({ basePmol: 0.05,
                     fragments: [{ bp: 1000, ratio: 1 }] }).rows[0].name,
       '断片 1');
  });

  test('Gibson: 0.5 µL 未満はピペットの誤差を注意する', function () {
    var result = GIBSON.plan({ basePmol: 0.05, fragments: [
      { name: '濃い挿入', bp: 1000, ratio: 1, concentration: 500 }
    ] });
    eq(result.rows[0].ul, 0.07);   // 32.5 ng / 500 = 0.065 -> 0.07
    ok(result.notes.join(' ').indexOf('濃い挿入') !== -1,
       result.notes.join(' '));
  });

  test('Gibson: 体積の上限を書いたときだけ、超過を知らせる', function () {
    eq(GIBSON.plan(SPEC).notes, []);
    var over = GIBSON.plan({ basePmol: 0.05, dnaVolumeLimitUl: 5,
                             fragments: SPEC.fragments });
    ok(over.notes.join(' ').indexOf('6.5 µL') !== -1, over.notes.join(' '));
  });

  test('Gibson: 断片が無ければ計算しない', function () {
    ok(message(function () { GIBSON.plan({ basePmol: 0.05, fragments: [] }); })
       .indexOf('断片') !== -1);
  });

  test('Gibson: モル比 1 にあたる pmol が無ければ計算しない', function () {
    ok(message(function () {
      GIBSON.plan({ fragments: [{ bp: 1000, ratio: 1 }] });
    }).indexOf('pmol') !== -1);
  });

  // ---- 反応の表（Amount を手で入れる） --------------------------------

  // The bench sheet of 2026-10-10: every number below is on it.
  var BENCH = { fragments: [
    { concentration: 166, bp: 5238, ul: 0.35 },
    { concentration: 100, bp: 308, ul: 0.8 }
  ] };

  test('反応: pmol/µL・pmol・ng は実験ノートの表と一致する', function () {
    var rows = GIBSON.reaction(BENCH).rows;
    eq(rows.map(function (r) { return r.name; }), ['Vector', 'Insert']);
    eq(rows.map(function (r) { return r.pmolPerUl; }), [0.0488, 0.4995]);
    eq(rows.map(function (r) { return r.pmol; }), [0.0171, 0.3996]);
    eq(rows[0].ng, 58.1);
    eq(rows[1].ng, 80);          // 0.80 µL x 100 ng/µL
  });

  test('反応: Master Mix は 2× なので全量の半分、miliQ は残り', function () {
    var result = GIBSON.reaction(BENCH);
    eq([result.mixUl, result.dnaUl, result.waterUl], [10, 1.15, 8.85]);
    eq(result.totals, { pmol: 0.4167, ul: 20, pmolOk: true });
  });

  test('反応: 束縛条件は判定するだけで、量を変えない', function () {
    var checks = GIBSON.reaction(BENCH).checks;
    eq(checks.map(function (c) { return c.ok; }), [true, false, true, true]);
    eq(checks[1].value, '×1.38');   // 80 / 58.1
    eq(GIBSON.reaction(BENCH).rows[1].ul, 0.8);
  });

  test('反応: 数値ごとの判定（色分け用）は束縛条件の判定と一致する', function () {
    var result = GIBSON.reaction(BENCH);
    eq(result.rows[0].ngOk, true);          // Vector 58.1 ng
    eq(result.rows[1].ratioOk, false);      // x1.38
    eq(result.totals.pmolOk, true);         // 0.4167 pmol
    eq(result.rows.map(function (r) { return r.ulOk; }), [true, true]);
    eq(GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 1 }, { concentration: 100, bp: 1000, ul: 1 },
      { concentration: 100, bp: 1000, ul: 1 }] }).totals.pmolOk, true);   // 3 fragment
  });

  test('条件設定の一覧は判定に使う値から作る', function () {
    var text = GIBSON.guide().map(function (g) { return g.label + ' ' + g.target; }).join('\n');
    ['50–100 ng', '×2–3', '×5 以上', '2–3 fragment', '0.02–0.5 pmol', '4–6 fragment',
     '0.2–1 pmol'].forEach(function (piece) {
      ok(text.indexOf(piece) !== -1, piece + ' が一覧に無い: ' + text);
    });
    // 載せないよう頼まれた 3 項目（判定はする）
    ['0.22', '10 µL', '20 µL', '1–2 fragment'].forEach(function (piece) {
      ok(text.indexOf(piece) === -1, piece + ' は一覧に載せない: ' + text);
    });
  });

  test('反応: 200 bp 未満の Insert は ×5 以上で判定する', function () {
    var result = GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 0.5 },     // 50 ng
      { concentration: 100, bp: 199, ul: 2.5 }       // 250 ng = x5
    ] });
    eq(result.checks[1].target, '×5 以上');
    eq(result.checks[1].ok, true);
  });

  test('反応: 3 fragment は 2–3 fragment の範囲（0.02–0.5 pmol）で判定する', function () {
    var result = GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 1 },
      { concentration: 100, bp: 1000, ul: 1 },
      { concentration: 100, bp: 1000, ul: 1 }
    ] });
    eq(result.rows.map(function (r) { return r.name; }),
       ['Vector', 'Insert A', 'Insert B']);
    eq(result.checks[3].target, '0.02–0.5 pmol');
    eq(result.checks[3].ok, true);       // 0.0308 + 0.1538 x 2 = 0.3385 pmol
  });

  test('反応: Insert に名前を付けられる（無ければ Insert / Insert A, B…）', function () {
    eq(GIBSON.fragmentLabel(0, 4), 'Vector');
    eq(GIBSON.fragmentLabel(1, 2), 'Insert');
    eq(GIBSON.fragmentLabel(3, 4), 'Insert C');
    var result = GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 1 },
      { name: 'GFP', concentration: 100, bp: 1000, ul: 1 },
      { name: '', concentration: 100, bp: 1000, ul: 1 }] });
    eq(result.rows.map(function (r) { return r.name; }), ['Vector', 'GFP', 'Insert B']);
    ok(result.checks[1].label.indexOf('GFP / Vector') === 0, result.checks[1].label);
  });

  test('反応: 4–6 fragment は 0.2–1.0 pmol で判定する', function () {
    var result = GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 1 }, { concentration: 100, bp: 1000, ul: 0.3 },
      { concentration: 100, bp: 1000, ul: 0.3 }, { concentration: 100, bp: 1000, ul: 0.3 }] });
    eq(result.checks[4].target, '0.2–1 pmol');
    eq(result.checks[4].ok, false);      // 0.0308 + 0.0462 x 3 = 0.1692 < 0.2
    eq(GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 1 }, { concentration: 100, bp: 1000, ul: 0.5 },
      { concentration: 100, bp: 1000, ul: 0.5 }, { concentration: 100, bp: 1000, ul: 0.5 }]
    }).checks[4].ok, true);              // 0.0308 + 0.0769 x 3 = 0.2615
  });

  test('反応: ピペットで取れるのは 0.22 µL 以上', function () {
    var make = function (ul) {
      return GIBSON.reaction({ fragments: [
        { concentration: 166, bp: 5238, ul: ul }, { concentration: 100, bp: 308, ul: 1 }] });
    };
    var last = function (r) { return r.checks[r.checks.length - 1]; };
    eq(last(make(0.22)).ok, true);
    eq(last(make(0.21)).ok, false);
    eq(last(make(0.21)).value, '0.21 µL（Vector）');
  });

  test('おすすめ: 0.22 µL 未満の量は出さない', function () {
    // 600 ng/µL の Vector は 50 ng が 0.08 µL。0.22 µL（132 ng）では多すぎる
    var answer = GIBSON.recommend({ fragments: [
      { concentration: 600, bp: 5000 }, { concentration: 100, bp: 1500 }] });
    ok(!answer.satisfied, answer.message);
    if (answer.result) {
      answer.result.rows.forEach(function (row) { ok(row.ul >= 0.22, row.name + ' ' + row.ul); });
    }
  });

  test('反応: 0.5 µL 未満でも注意は出さない（不要とのこと）', function () {
    eq(GIBSON.reaction(BENCH).notes, []);     // Vector 0.35 µL
  });

  test('反応: Insert と Vector の合計が 10 µL を超えたら計算しない', function () {
    var refused = message(function () {
      GIBSON.reaction({ fragments: [
        { concentration: 100, bp: 5000, ul: 6 }, { concentration: 100, bp: 1000, ul: 4.01 }] });
    });
    ok(refused.indexOf('10.01 µL') !== -1 && refused.indexOf('計算しません') !== -1, refused);
    // ちょうど 10 µL は計算する（miliQ 0 µL）
    eq(GIBSON.reaction({ fragments: [
      { concentration: 100, bp: 5000, ul: 6 }, { concentration: 100, bp: 1000, ul: 4 }]
    }).waterUl, 0);
  });

  test('反応: 空欄は何が足りないかを言う', function () {
    ok(message(function () {
      GIBSON.reaction({ fragments: [
        { concentration: 166, bp: 5238, ul: '' }, { concentration: 100, bp: 308, ul: 1 }
      ] });
    }).indexOf('Vector の Volume') !== -1);
  });

  test('反応: 全量は常に 20 µL（入力しない）', function () {
    eq(GIBSON.TOTAL_UL, 20);
    eq(GIBSON.reaction({ totalUl: 10, fragments: BENCH.fragments }).totals.ul, 20);
  });

  // ---- Vector の量と倍数から計算する方式 ------------------------------

  test('倍数方式: Vector 0.35 µL（58.1 ng）・×2.5 → Insert 1.45 µL', function () {
    // 0.35 x 166 = 58.1 ng、58.1 x 2.5 = 145.25 ng / 100 = 1.4525 -> 1.45 µL
    var result = GIBSON.amountsFromRatio({ fragments: [
      { concentration: 166, bp: 5238, ul: 0.35 },
      { concentration: 100, bp: 1500, ratio: 2.5 }] });
    eq(result.vectorNg, 58.1);
    eq(result.rows.map(function (r) { return r.ul; }), [0.35, 1.45]);
    eq(result.targets[1].targetNg, 145.3);    // 145.25 -> 4 significant figures
    eq(result.targets[1].actualNg, 145);      // 1.45 µL x 100
    eq(result.targets[1].actualRatio, 2.5);   // 145 / 58.1 = 2.496 -> 2.5
    eq(result.checks.map(function (c) { return c.ok; }), [true, true, true, true]);
  });

  test('倍数方式: 結果は手入力の表と同じ判定になる', function () {
    var result = GIBSON.amountsFromRatio({ fragments: [
      { concentration: 166, bp: 5238, ul: '.35' }, { concentration: 100, bp: 308, ratio: 1.38 }] });
    var byHand = GIBSON.reaction({ fragments: [
      { concentration: 166, bp: 5238, ul: result.rows[0].ul },
      { concentration: 100, bp: 308, ul: result.rows[1].ul }] });
    eq(result.checks, byHand.checks);
  });

  test('倍数方式: Vector の Volume・倍数が無ければ何が足りないかを言う', function () {
    ok(message(function () {
      GIBSON.amountsFromRatio({ fragments: [
        { concentration: 166, bp: 5238, ul: '' }, { concentration: 100, bp: 308, ratio: 2 }] });
    }).indexOf('Vector の Volume') !== -1);
    ok(message(function () {
      GIBSON.amountsFromRatio({ fragments: [
        { concentration: 166, bp: 5238, ul: 0.3 }, { concentration: 100, bp: 308, ratio: '' }] });
    }).indexOf('Insert の倍数') !== -1);
  });

  test('倍数方式: 10 µL を超える量になったら計算しない', function () {
    ok(message(function () {
      GIBSON.amountsFromRatio({ fragments: [
        { concentration: 20, bp: 5000, ul: 5 }, { concentration: 20, bp: 1000, ratio: 3 }] });
    }).indexOf('計算しません') !== -1);
  });

  // ---- おすすめの Amount --------------------------------------------

  function allMet(result) {
    return result.checks.every(function (c) { return c.ok !== false; });
  }

  test('おすすめ: 満たせるときは、すべての束縛条件を満たす Amount を返す', function () {
    var answer = GIBSON.recommend({ fragments: [
      { concentration: 166, bp: 5238 }, { concentration: 100, bp: 1500 }] });
    eq([answer.satisfied, answer.relaxed], [true, false]);
    ok(allMet(answer.result), JSON.stringify(answer.result.checks));
  });

  test('おすすめ: 返した Amount は、表に入れて計算しても同じ判定になる', function () {
    var spec = { fragments: [
      { concentration: 166, bp: 5238 }, { concentration: 100, bp: 1500 }] };
    var answer = GIBSON.recommend(spec);
    var again = GIBSON.reaction({ fragments: spec.fragments.map(function (f, i) {
      return { concentration: f.concentration, bp: f.bp, ul: answer.result.rows[i].ul };
    }) });
    eq(again.checks, answer.result.checks);
  });

  test('おすすめ: 短い Insert で倍数を守れないときは、倍数だけを下げる', function () {
    // 実験ノートの例。Vector 50 ng・×2 でも 0.514 pmol で 0.5 を超える
    var answer = GIBSON.recommend({ fragments: [
      { concentration: 166, bp: 5238 }, { concentration: 100, bp: 308 }] });
    eq([answer.satisfied, answer.relaxed], [false, true]);
    var checks = answer.result.checks;
    eq(checks.map(function (c) { return c.ok; }), [true, false, true, true]);
    ok(answer.result.rows[1].ratio < 2, '倍数は下限より下');
    ok(answer.message.indexOf('倍数を目安の下限から') !== -1, answer.message);
  });

  test('おすすめ: 解が 2 通りしかない入力でも見つける（330 bp の Insert）', function () {
    // 0.01 µL 刻みの総当たりで、全条件を満たすのは 2 通りだけ（2026-10-10 確認）
    var answer = GIBSON.recommend({ fragments: [
      { concentration: 166, bp: 5238 }, { concentration: 100, bp: 330 }] });
    eq(answer.satisfied, true);
    ok(allMet(answer.result), JSON.stringify(answer.result.checks));
  });

  test('おすすめ: Insert が 3 本でも、それぞれの倍数を別々に合わせて満たす', function () {
    var answer = GIBSON.recommend({ fragments: [
      { concentration: 100, bp: 5000 }, { concentration: 50, bp: 1500 },
      { concentration: 200, bp: 600 }, { concentration: 80, bp: 2500 }] });
    eq(answer.satisfied, true);
    ok(allMet(answer.result), JSON.stringify(answer.result.checks));
  });

  test('おすすめ: どうしても満たせないときは、満たせない条件を言う', function () {
    var answer = GIBSON.recommend({ fragments: [
      { concentration: 5, bp: 5000 }, { concentration: 5, bp: 1500 }] });
    eq(answer.result, null);
    ok(answer.message.indexOf('Vector') !== -1, answer.message);
  });

  test('おすすめ: Amount は空欄でよいが、濃度は要る', function () {
    ok(message(function () {
      GIBSON.recommend({ fragments: [
        { concentration: '', bp: 5238 }, { concentration: 100, bp: 308 }] });
    }).indexOf('Vector の濃度') !== -1);
  });

  // ---- Tm ------------------------------------------------------------

  test('Gibson: Tm はエンジンの式そのもの（20 mer, GC 50% で 55℃）',
       function () {
    var row = GIBSON.tmOf('ATGCATGCATGCATGCATGC');
    eq(row.tm, 55);
    eq(row.tm, PB.tmTakara('ATGCATGCATGCATGCATGC'));
    eq(row.length, 20);
    eq(row.gcPercent, 50);
    eq(row.formulaApplies, true);
    eq(row.note, null);
  });

  test('Gibson: 小文字・空白・数字は取り除いて読む', function () {
    eq(GIBSON.tmOf(' 1 atg c ').sequence, 'ATGC');
    eq(GIBSON.tmOf('atgc').tm, 7);   // AT 2 -> 4, GC 2 -> 8, 4 + 8 - 5
  });

  test('Gibson: U は T として読む', function () {
    eq(GIBSON.tmOf('AUGC').sequence, 'ATGC');
  });

  test('Gibson: 25 mer を超えたら、値は出したうえで適用範囲外と言う',
       function () {
    var row = GIBSON.tmOf('ATGCATGCATGCATGCATGCATGCAT');   // 26 mer
    eq(row.length, 26);
    eq(row.formulaApplies, false);
    ok(row.note.indexOf('26 mer') !== -1, row.note);
    ok(typeof row.tm === 'number', 'Tm は隠さない');
  });

  test('Gibson: A/C/G/T 以外はエンジンのメッセージのまま断る', function () {
    ok(message(function () { GIBSON.tmOf('ATGN'); })
       .indexOf('A/C/G/T') !== -1);
  });

  test('Gibson: 空の配列は断る', function () {
    ok(message(function () { GIBSON.tmOf('   '); })
       .indexOf('配列') !== -1);
  });

  test('Gibson: 1 行 1 本で、空行は飛ばす', function () {
    var rows = GIBSON.tmRows('atgc\n\n  \nATGCATGC\n');
    eq(rows.map(function (r) { return r.length; }), [4, 8]);
  });

  // ---- PCR 条件（貼られたプロトコール, 2026-10-10） ------------------

  test('PCR: Tm 55℃ 以上は 5 sec、未満は 15 sec。全体は長いほう', function () {
    var rows = GIBSON.tmRows('ATGCATGCATGCATGCATGC\nATGCATGC');   // 55, 19
    var pcr = GIBSON.pcrConditions(rows);
    eq(pcr.perPrimerSec, [5, 15]);
    eq(pcr.anneal, { tempC: 55, sec: 15 });
    eq(pcr.mode, '3 step');
    eq(pcr.denature.tempC, 98);
  });

  test('PCR: 25 mer を超えるプライマーは 5 sec', function () {
    var pcr = GIBSON.pcrConditions(GIBSON.tmRows(
      'AAAAAAAAAAAAAAAAAAAAAAAAAA'));                  // 26 mer, Tm 47
    eq(pcr.perPrimerSec, [5]);
  });

  test('PCR: Tm 70℃ 以上（25 mer 以下）なら 2 step を勧める', function () {
    // 20 G = 75℃
    var pcr = GIBSON.pcrConditions(GIBSON.tmRows(
      'GGGGGGGGGGGGGGGGGGGG\nATGCATGCATGCATGCATGC'));
    eq(pcr.twoStepRows, [1]);
    ok(pcr.notes[0].indexOf('2 step') !== -1, pcr.notes[0]);
  });

  test('PCR: 25 mer 超は式の範囲外なので 70℃ の判定に使わない', function () {
    var pcr = GIBSON.pcrConditions(GIBSON.tmRows(new Array(41).join('GC')));
    eq(pcr.twoStepRows, []);
  });

  test('Gibson: 2 本の Tm 差を出す', function () {
    // ATGC = 7, ATGCATGC = 4*2 + 4*4 - 5 = 19
    eq(GIBSON.tmDifference('ATGC', 'ATGCATGC'), 12);
  });
})();
