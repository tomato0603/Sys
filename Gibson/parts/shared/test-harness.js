/* Self-test suite.
 *
 * Every expected value here was worked out by hand, not produced by the code
 * under test. That is the point: the Takara Tm formula, the PrimeSTAR HS
 * condition table and the junction equation are the places where a silent
 * error would emit unusable primers without ever failing.
 *
 * The whole suite runs in the browser from the ページ上のボタン, so the tool
 * can be checked without trusting anyone's word for it.
 */
var PBTESTS = (function () {
  'use strict';

  var cases = [];
  function test(name, fn) { cases.push({ name: name, fn: fn }); }

  function eq(actual, expected, what) {
    var a = JSON.stringify(actual), b = JSON.stringify(expected);
    if (a !== b) {
      throw new Error((what ? what + ': ' : '') + '期待 ' + b + ' / 実際 ' + a);
    }
  }
  function ok(value, what) {
    if (!value) { throw new Error((what || '条件') + ' が成り立ちません'); }
  }
  function throws(fn, what) {
    try { fn(); } catch (e) { return; }
    throw new Error((what || '呼び出し') + ' は例外になるべきです');
  }

  // A 60 bp circular test record; count the chunks of ten to check coordinates.
  // ATGAAACGTT GCGAATAAGG CCTTAAGGCC AATTGGCCAA TTCCGGAATT GGCCAATTGC
  var TEST_SEQ = 'ATGAAACGTTGCGAATAAGGCCTTAAGGCCAATTGGCCAATTCCGGAATTGGCCAATTGC';
  var TEST_GENBANK = [
    'LOCUS       TESTPLASMID               60 bp    DNA        circular    26-AUG-2026',
    'DEFINITION  synthetic record used only for tests.',
    'FEATURES             Location/Qualifiers',
    '     CDS             1..18',
    '                     /gene="testA"',
    '     CDS             complement(21..38)',
    '                     /gene="testB"',
    'ORIGIN',
    '        1 atgaaacgtt gcgaataagg ccttaaggcc aattggccaa ttccggaatt ggccaattgc',
    '//', ''].join('\n');
  var SMALL_CDS = 'ATGAAACGTTGCGAATAA';        // M K R C E *
  var SMALL_PROTEIN = 'MKRCE*';

  // ---- Tm (Takara Bio) ------------------------------------------------
  test('Tm: ATGC は 2*2 + 4*2 - 5 = 7', function () {
    eq(PB.tmTakara('ATGC'), 7);
  });
  test('Tm: 20 A は 40 - 5 = 35', function () {
    eq(PB.tmTakara(new Array(21).join('A')), 35);
  });
  test('Tm: 20 G は 80 - 5 = 75', function () {
    eq(PB.tmTakara(new Array(21).join('G')), 75);
  });
  test('Tm: A5T5G5C5 の 20mer はちょうど 55（時間規則の境界）', function () {
    eq(PB.tmTakara('ATGCATGCATGCATGCATGC'), 55);
  });
  test('Tm: 不確定塩基は例外（重みが定義されていないため）', function () {
    throws(function () { PB.tmTakara('ATGN'); });
    throws(function () { PB.tmTakara(''); });
  });
  test('Tm: 式の適用範囲は 25 mer まで', function () {
    eq(PB.tmFormulaApplies(25), true);
    eq(PB.tmFormulaApplies(26), false);
  });

  // ---- annealing / extension time ------------------------------------
  test('アニーリング時間: Tm>=55 は 5 秒、<55 は 15 秒、25mer 超は 5 秒', function () {
    eq(PB.annealingTimeSec(20, 55), 5);
    eq(PB.annealingTimeSec(20, 54), 15);
    eq(PB.annealingTimeSec(30, 40), 5);
  });
  test('アニーリング時間: 規則が矛盾するときだけ注意を出す', function () {
    ok(PB.annealingTimeNote(30, 40) !== null);
    eq(PB.annealingTimeNote(30, 60), null);
    eq(PB.annealingTimeNote(20, 40), null);
  });
  test('伸長時間: 1 分/kb と 30 秒単位の切り上げ', function () {
    eq(PB.extensionTimeSec(5000), 300);
    eq(PB.practicalExtensionTimeSec(5000), 300);
    eq(PB.practicalExtensionTimeSec(2686), 180);   // 161.16 -> 180
    eq(PB.practicalExtensionTimeSec(5500), 330);
    throws(function () { PB.extensionTimeSec(0); });
  });
  test('時間の表示', function () {
    eq(PB.formatSeconds(300), '5 分');
    eq(PB.formatSeconds(330), '5 分 30 秒');
    eq(PB.formatSeconds(45), '45 秒');
  });

  // ---- genetic code ---------------------------------------------------
  test('遺伝暗号表: 64 コドン、終止は TAA/TAG/TGA', function () {
    eq(Object.keys(PB.GENETIC_CODE).length, 64);
    eq(PB.GENETIC_CODE.ATG, 'M');
    eq(PB.GENETIC_CODE.TGG, 'W');
    var stops = Object.keys(PB.GENETIC_CODE).filter(function (c) {
      return PB.GENETIC_CODE[c] === '*';
    }).sort();
    eq(stops, ['TAA', 'TAG', 'TGA']);
  });
  test('翻訳: ORF と端数の報告', function () {
    eq(PB.translate('ATGAAACGTTGCGAATAA').protein, 'MKRCE*');
    eq(PB.translate('ATGAA'), { protein: 'M', trailing: 2 });
    eq(PB.translate('ATGNNN').protein, 'MX');
  });
  test('逆相補', function () {
    eq(PB.reverseComplement('ATGC'), 'GCAT');
    eq(PB.reverseComplement('CCTTAAGGCCAATTGGCC'), 'GGCCAATTGGCCTTAAGG');
  });

  // ---- occurrence counting -------------------------------------------
  test('出現数: 回文は 1 回として数える（制限酵素サイトは全部回文）', function () {
    eq(PB.countOccurrences('AAGGATCCAA', 'GGATCC'), 1);
  });
  test('出現数: 非回文は両鎖を数える／重なりも数える', function () {
    eq(PB.countOccurrences('GGGAAACCC', 'GGG'), 2);
    eq(PB.countOccurrences('AAAA', 'AA'), 3);
  });
  test('出現数: 原点をまたぐ一致', function () {
    eq(PB.countOccurrences('TCCAAAAAAGGA', 'GGATCC'), 0);
    eq(PB.countOccurrences('TCCAAAAAAGGA', 'GGATCC', true), 1);
  });

  // ---- circular coordinates ------------------------------------------
  test('環状座標: 原点をまたぐ切り出し（58..3 は T,G,C,A,T,G）', function () {
    var circular = new PB.Record({ seq: TEST_SEQ, isCircular: true });
    var linear = new PB.Record({ seq: TEST_SEQ, isCircular: false });
    eq(circular.subseq(1, 4), 'ATGA');
    eq(circular.subseq(58, 3), 'TGCATG');
    throws(function () { linear.subseq(58, 3); }, '直鎖でのまたぎ');
    throws(function () { circular.subseq(0, 5); });
  });
  test('環状座標: 原点をまたぐ検索と GC 含量', function () {
    var circular = new PB.Record({ seq: TEST_SEQ, isCircular: true });
    var linear = new PB.Record({ seq: TEST_SEQ, isCircular: false });
    eq(circular.findAll('GCATG'), [59]);
    eq(linear.findAll('GCATG'), []);
    eq(linear.findAll('AATT'), [31, 39, 47, 55]);
    eq(circular.gcContent().toFixed(6), (28 / 60).toFixed(6));
  });

  // ---- location parsing ----------------------------------------------
  test('location: 範囲・complement・join・あいまい記号', function () {
    eq(PB.parseLocation('1..18'), { segments: [[1, 18]], strand: 1 });
    eq(PB.parseLocation('complement(21..38)'),
       { segments: [[21, 38]], strand: -1 });
    eq(PB.parseLocation('join(2600..2686,1..120)'),
       { segments: [[2600, 2686], [1, 120]], strand: 1 });
    eq(PB.parseLocation('complement(join(10..20,30..40))'),
       { segments: [[10, 20], [30, 40]], strand: -1 });
    eq(PB.parseLocation('<10..>30'), { segments: [[10, 30]], strand: 1 });
    eq(PB.parseLocation('42'), { segments: [[42, 42]], strand: 1 });
    throws(function () { PB.parseLocation('10^11'); });
  });

  // ---- GenBank / ApE --------------------------------------------------
  test('GenBank: LOCUS・CDS・逆鎖の読み取り', function () {
    var record = PB.parseGenbank(TEST_GENBANK);
    eq(record.name, 'TESTPLASMID');
    eq(record.length(), 60);
    eq(record.isCircular, true);
    eq(record.seq, TEST_SEQ);
    var cds = record.cdsFeatures();
    eq(cds.length, 2);
    eq(cds[0].segments, [[1, 18]]);
    eq(cds[0].qualifiers.gene, 'testA');
    eq(cds[1].strand, -1);
    eq(record.featureSequence(cds[0]), SMALL_CDS);
    eq(record.featureSequence(cds[1]), 'GGCCAATTGGCCTTAAGG');
  });
  test('入力形式: .ape は GenBank として読む／GenBank でない .ape は拒否', function () {
    eq(PB.detectFormat(TEST_GENBANK, 'x.ape'), 'genbank');
    throws(function () { PB.detectFormat('not genbank\n', 'x.ape'); });
  });
  test('入力形式: FASTA は対象外として拒否する', function () {
    throws(function () { PB.detectFormat('>seq\nACGT\n', 'x.txt'); }, 'FASTA');
    throws(function () { PB.loadFromText('>seq\nACGTACGT\n', 'x.txt'); });
  });
  test('入力形式: タンパク質配列のファイルは拒否する', function () {
    // 文字種の検査だけでは通ってしまうので、A/C/G/T の割合で判定する
    throws(function () {
      PB.loadFromText('MDIPLPLRIREVLADLGLNVVRSVMFGDSATALAK', 'x.txt');
    }, 'タンパク質配列');
    throws(function () {
      PB.requireDna('MDIPLPLRIREVLADLGLNVVRSV', 'テスト');
    });
  });
  test('入力形式: プレーンテキストの塩基配列は受け付ける', function () {
    var record = PB.loadFromText('atgaaacgttgcgaataa' + 'ACGT'.repeat(50), 'x.txt');
    eq(record.seq.slice(0, 18), 'ATGAAACGTTGCGAATAA');
    eq(record.isCircular, true);
    ok(record.warnings.length > 0, '注釈がない旨の注意');
  });
  test('ApE: LOCUS 行の桁位置（ApE の出力から実測）', function () {
    eq(PB.formatLocusLine('pET_LgBit-his', 5715, false, '28-AUG-2025'),
       'LOCUS       pET_LgBit-his           5715 bp    DNA'
       + '        linear       28-AUG-2025');
    var line = PB.formatLocusLine('pTEST', 900, true, '26-AUG-2026');
    eq(line.substr(PB.LOCUS_TOPOLOGY_COL, 8), 'circular');
    eq(line.slice(PB.LOCUS_DATE_COL), '26-AUG-2026');
  });
  test('ApE: 修飾子の折り返しは空白位置のみ（読み直しで復元できる）', function () {
    var value = 'arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0';
    var lines = PB.formatQualifier('ApEinfo_graphicformat', value);
    eq(lines.length, 2);
    lines.forEach(function (l) { ok(l.length <= PB.GENBANK_LINE_WIDTH); });
    eq(lines[0].trim() + ' ' + lines[1].trim(),
       '/ApEinfo_graphicformat="' + value + '"');
  });
  test('ApE: 書き出して読み直すと同じ内容になる', function () {
    var record = PB.parseGenbank(TEST_GENBANK);
    var again = PB.parseGenbank(PB.formatGenbank(record, '26-AUG-2026'));
    eq(again.seq, record.seq);
    eq(again.isCircular, record.isCircular);
    eq(again.features.map(function (f) { return f.segments; }),
       record.features.map(function (f) { return f.segments; }));
    eq(again.features.map(function (f) { return f.qualifiers; }),
       record.features.map(function (f) { return f.qualifiers; }));
  });
  test('ApE: ORIGIN の桁組み', function () {
    var lines = PB.formatOrigin(new Array(131).join('A'));
    eq(lines.length, 3);
    eq(lines[0], '        1 ' + new Array(7).join('aaaaaaaaaa ').trim());
    ok(lines[2].indexOf('      121 ') === 0);
  });

  // ---- edits ----------------------------------------------------------
  test('変異表記: 4 種類の書き方', function () {
    eq(PB.parseEdits('K42A')[0].kind(), 'substitution');
    var ins = PB.parseEdits('ins@50:GGSGGS')[0];
    eq([ins.start, ins.end], [51, 50]);       // 空区間 = 残基 51 の直前
    eq(ins.kind(), 'insertion');
    eq(PB.parseEdits('del45-50')[0].kind(), 'deletion');
    eq(PB.parseEdits('del45')[0].kind(), 'deletion');
    eq(PB.parseEdits('rep45-50:GS')[0].kind(), 'replacement');
    eq(PB.parseEdits('del45-47, K10A, ins@60:WW').map(function (e) {
      return e.start;
    }), [10, 45, 61]);
    throws(function () { PB.parseEdits('K42'); });
    throws(function () { PB.parseEdits('del50-45'); });
    throws(function () { PB.parseEdits('   '); });
  });
  test('検証: 元の残基が違えば止まる（残基番号ずれの検出）', function () {
    var edits = PB.parseEdits('del2-4');
    PB.validateEdits(edits, SMALL_PROTEIN);
    eq(edits[0].original, 'KRC');
    throws(function () {
      PB.validateEdits(PB.parseEdits('A2G'), SMALL_PROTEIN);
    }, '残基 2 は K なので A2G');
    throws(function () {
      PB.validateEdits(PB.parseEdits('K99A'), SMALL_PROTEIN);
    });
    throws(function () {
      PB.validateEdits(PB.parseEdits('del2-4,rep3-5:WW'), SMALL_PROTEIN);
    }, '範囲の重なり');
  });
  test('適用: 置換・挿入・欠失・一般置換・複合', function () {
    function apply(notation) {
      var edits = PB.parseEdits(notation);
      PB.validateEdits(edits, SMALL_PROTEIN);
      return PB.applyEdits(SMALL_PROTEIN, edits);
    }
    eq(apply('K2A'), 'MARCE*');
    eq(apply('ins@2:WW'), 'MKWWRCE*');
    eq(apply('del2-4'), 'ME*');
    eq(apply('rep2-4:Y'), 'MYE*');
    eq(apply('K2A,ins@4:WW'), 'MARCWWE*');
  });
  test('アライメント: 差分を適用すると必ず目的の配列になる', function () {
    [['MKRCEAAGGWY', 'MKACEAAGGWY'],
     ['MKRCEAAGGWY', 'MKRCEAAQQGGWY'],
     ['MKRCEAAGGWY', 'MKRCEGGWY'],
     ['MKRCEAAGGWY', 'MKACEAAGGWWWY'],
     ['MKRCEAAGGWY', 'MKRCEAAGGWY']].forEach(function (pair) {
      var edits = PB.diffProteins(pair[0], pair[1]);
      eq(PB.applyEdits(pair[0], edits), pair[1], pair[1]);
    });
    eq(PB.diffProteins('MKRCE', 'MKRCE'), []);
  });

  // ---- codon selection ------------------------------------------------
  test('コドン選択: 障害がなければ最優先コドン', function () {
    var choice = PB.chooseCodons('MKRCE');
    eq(choice.dna, 'ATGAAACGCTGCGAA');
    eq(choice.violations, 0);
  });
  test('コドン選択: GS リンカーも制約違反なしで通る（貪欲法では破綻する）', function () {
    var choice = PB.chooseCodons('GGSGGSGGS');
    eq(choice.violations, 0);
    ok(PB.longestGcRun(choice.dna) <= PB.MAX_GC_RUN, 'GC 連続');
    ok(PB.longestHomopolymerRun(choice.dna) <= PB.MAX_HOMOPOLYMER_RUN, '同一塩基連続');
    eq(PB.translate(choice.dna).protein, 'GGSGGSGGS');
  });
  test('コドン選択: 翻訳は必ず指定したペプチドに戻る', function () {
    ['HHHHHH', 'GGSGGS', 'MKRCEAAWY', 'PPPP', 'AAAA', 'RRRR'].forEach(function (p) {
      eq(PB.translate(PB.chooseCodons(p).dna).protein, p, p);
    });
  });
  test('コドン選択: 避けられない違反は隠さず数を報告する', function () {
    var choice = PB.chooseCodons('GG', { upstream: 'GTTGGG' });
    ok(choice.violations > 0, '違反数');
    ok(choice.level.indexOf('制約違反') === 0);
  });
  test('コドン選択: 回避配列を作らない／レアコドンを避ける', function () {
    var choice = PB.chooseCodons('EF', { avoid: ['GAATTC'] });
    eq(choice.violations, 0);
    ok(choice.dna.indexOf('GAATTC') === -1);
    PB.chooseCodons('RRRR').picks.forEach(function (pick) {
      ok(PB.ECOLI_RARE_CODONS.indexOf(pick.codon) === -1, pick.codon);
    });
  });
  test('コドン選択: 未知のアミノ酸・未対応の宿主は例外', function () {
    throws(function () { PB.chooseCodons('MKZ'); });
    throws(function () { PB.chooseCodons('M', { host: 'yeast' }); });
  });

  // ---- grouping -------------------------------------------------------
  test("5' テールの上限 = 2 x (60 - 18) = 84", function () {
    eq(PB.tailBudgetNt(), 84);
    eq(PB.tailBudgetNt(40, 20), 40);
    throws(function () { PB.tailBudgetNt(20, 20); });
  });
  test('区間: 変更しない残基は元のコドンを流用する', function () {
    var edits = PB.parseEdits('K2A,E5W');
    PB.validateEdits(edits, SMALL_PROTEIN);
    var group = PB.buildEditGroup(SMALL_PROTEIN, SMALL_CDS, edits, {});
    eq([group.spanStart, group.spanEnd], [2, 5]);
    eq(group.newPeptide, 'ARCW');
    eq(group.newDna, 'GCGCGTTGCTGG');      // GCG + CGTTGC(流用) + TGG
    eq(group.keptNt, 6);
    eq(group.chosenNt, 6);
    eq(PB.translate(group.newDna).protein, 'ARCW');
  });
  test('区間: 欠失はテールを必要としない', function () {
    var edits = PB.parseEdits('del2-4');
    PB.validateEdits(edits, SMALL_PROTEIN);
    var group = PB.buildEditGroup(SMALL_PROTEIN, SMALL_CDS, edits, {});
    eq(group.newDna, '');
    eq(group.requiredTailNt, 0);
    eq(group.removedNt, 9);
  });
  test('ラウンド分け: 近い編集はまとめ、離れた編集は分ける', function () {
    var protein = 'M' + new Array(201).join('A') + '*';
    var cds = 'ATG' + new Array(201).join('GCG') + 'TAA';
    var edits = PB.parseEdits('A2W,A150W');
    PB.validateEdits(edits, protein);
    eq(PB.groupEdits(protein, edits, { cdsDna: cds }).length, 2);
    var near = PB.parseEdits('A2W,A12W');
    PB.validateEdits(near, protein);
    eq(PB.groupEdits(protein, near, { cdsDna: cds, budget: 33 }).length, 1);
    eq(PB.groupEdits(protein, near, { cdsDna: cds, budget: 32 }).length, 2);
  });

  // ---- CDS coordinates -------------------------------------------------
  test('CDS 座標: 順鎖・逆鎖・原点をまたぐ join', function () {
    var forward = new PB.Feature('CDS', [[1, 18]], 1);
    var map = PB.cdsPositionMap(forward);
    eq(PB.residueToPlasmid(map, 1), [1, 2, 3]);
    eq(PB.residueToPlasmid(map, 6), [16, 17, 18]);
    var reverse = new PB.Feature('CDS', [[21, 38]], -1);
    eq(PB.cdsPositionMap(reverse).slice(0, 3), [38, 37, 36]);
    var wrapped = new PB.Feature('CDS', [[58, 60], [1, 3]], 1);
    eq(PB.cdsPositionMap(wrapped), [58, 59, 60, 1, 2, 3]);
    throws(function () { PB.residueToPlasmid(map, 7); });
  });

  return { cases: cases, TEST_SEQ: TEST_SEQ, TEST_GENBANK: TEST_GENBANK,
           SMALL_CDS: SMALL_CDS, SMALL_PROTEIN: SMALL_PROTEIN,
           test: test, eq: eq, ok: ok, throws: throws };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = PBTESTS; }

