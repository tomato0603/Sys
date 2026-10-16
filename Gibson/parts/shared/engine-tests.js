/* Part 2: design, verification and diagnostic sites.
 *
 * These need a template whose 18-mers are all unique on both strands, so the
 * synthetic example plasmid is embedded here. It is generated, not real: no
 * one's construct is in this file.
 */
(function () {
  'use strict';
  var test = PBTESTS.test, eq = PBTESTS.eq, ok = PBTESTS.ok,
      throws = PBTESTS.throws;

  var EXAMPLE_GENBANK = "LOCUS       pTEST-INV                900 bp    DNA     circular SYN 26-AUG-2026\nDEFINITION  synthetic example plasmid for Primer Builder. Not a real construct.\nACCESSION   .\nKEYWORDS    .\nSOURCE      synthetic DNA construct\n  ORGANISM  synthetic DNA construct\nCOMMENT     Generated deterministically; every 18-mer is unique on both strands.\nFEATURES             Location/Qualifiers\n     source          1..900\n                     /organism=\"synthetic DNA construct\"\n                     /mol_type=\"other DNA\"\n     misc_feature    60..95\n                     /label=\"promoter (placeholder)\"\n     CDS             101..700\n                     /gene=\"targetA\"\n                     /product=\"synthetic target protein\"\n                     /note=\"edit target for inverse PCR\"\nORIGIN\n        1 accctcagac gctactaaag cacgactgct tagataagtc catcatcgac actgggggcg\n       61 gcatggtggg ccgacttcga gctttcgcct cgctgtattg atggacattc ccttacccct\n      121 tcgtataagg gaagttctag ctgatctagg cttaaatgtt gtacgatcag taatgtttgg\n      181 ggatagcgcc acggcattgg ctaaagctgg ctgccacgtg ggacaggtaa cctcccgtgc\n      241 tctccgctcg caagacgcgc ctgcccagtc ggctcgttcc ggtggcaatc gaagtagatt\n      301 tgcaccaata cagaaacatg atctctcatt tgccaggccg cagtactgcc cgagaccgcg\n      361 aacgaccgcc gcgtggtacg tcgacgaaaa gtcctggccc tgcctgaaat acgggttgtt\n      421 caagggtgaa ttcctccaca tctctctatt tgccagtcca tgcgatcggg cgatacccct\n      481 acctcctcga aattcgtcta tccccacccg gataatatat ctgcaacgtc tgctttatct\n      541 atacgaggat gtgcgaacct atctatgctc attaccaccg cacatacaac tgtccgtgtt\n      601 cattgcgcat gtacgtttct atcccgaatg tatttcaaac agcagagccg aaccccgaag\n      661 cgacttcgcg cagactcacc cctgcgaaat tgcatcataa ctctagccca tggcgattaa\n      721 tttaaccgtt tcctttccgt gtcataaact agctcgcttc cgctgcaagc aggttgtggg\n      781 tctctacagc ttccaaactt tgagagctct tttccgatta cagcggcttt aatactggaa\n      841 atttttccgt gatgttgttc gaatataatc acgttcaatc aataactatt tcttttcaag\n//\n";
  PBTESTS.EXAMPLE_GENBANK = EXAMPLE_GENBANK;

  function example() {
    var record = PB.parseGenbank(EXAMPLE_GENBANK);
    var feature = record.cdsFeatures()[0];
    return { record: record, feature: feature,
             cds: record.featureSequence(feature),
             protein: PB.translate(record.featureSequence(feature)).protein };
  }
  PBTESTS.example = example;

  function oneGroup(context, notation) {
    var edits = PB.parseEdits(notation);
    PB.validateEdits(edits, context.protein);
    var groups = PB.groupEdits(context.protein, edits, { cdsDna: context.cds });
    if (groups.length !== 1) { throw new Error(notation + ' が 1 ラウンドになりません'); }
    return groups[0];
  }

  function onePair(context, notation) {
    var group = oneGroup(context, notation);
    var site = PB.editSite(context.record, context.feature, group);
    var pairs = PB.designPrimerPairs(context.record, site, {});
    if (!pairs.length) { throw new Error(notation + ' の候補がありません'); }
    return { group: group, site: site, pair: pairs[0], pairs: pairs };
  }

  test('例のプラスミド: 全 18-mer が両鎖で一意（設計の前提）', function () {
    var context = example();
    var seq = context.record.seq;
    var circular = seq + seq.slice(0, 24);
    var seen = {};
    for (var i = 0; i < seq.length; i++) {
      var window = circular.substr(i, 18);
      seen[window] = (seen[window] || 0) + 1;
      var complement = PB.reverseComplement(window);
      if (complement in seen) { seen[complement] += 1; }
    }
    var worst = 0;
    Object.keys(seen).forEach(function (k) { worst = Math.max(worst, seen[k]); });
    eq(worst, 1);
  });

  test('編集部位: 置換は残基 2 が 104..106（CDS は 101..700）', function () {
    var context = example();
    var group = oneGroup(context, 'D2A');
    var site = PB.editSite(context.record, context.feature, group);
    eq(site.leftPos, 103);
    eq(site.rightPos, 107);
    eq(site.removed, context.record.subseq(104, 106));
    eq(site.insert.length, 3);
  });

  test('編集部位: 挿入は何も除去しない／欠失は何も挿入しない', function () {
    var context = example();
    var insertion = PB.editSite(context.record, context.feature,
                                oneGroup(context, 'ins@6:GGSGGS'));
    eq(insertion.removed, '');
    eq(insertion.rightPos, insertion.leftPos + 1);
    eq(insertion.insert.length, 18);
    var deletion = PB.editSite(context.record, context.feature,
                               oneGroup(context, 'del5-8'));
    eq(deletion.insert, '');
    eq(deletion.removed.length, 12);
  });

  test('編集部位: 逆鎖 CDS では区間と挿入配列が反転する', function () {
    var record = new PB.Record({ seq: PBTESTS.TEST_SEQ, isCircular: true });
    var feature = new PB.Feature('CDS', [[21, 38]], -1);
    var protein = PB.translate(record.featureSequence(feature)).protein;
    var cds = record.featureSequence(feature);
    var edits = PB.parseEdits(protein[1] + '2W');
    PB.validateEdits(edits, protein);
    var group = PB.groupEdits(protein, edits, { cdsDna: cds })[0];
    var site = PB.editSite(record, feature, group);
    eq(site.removed, record.subseq(33, 35));
    eq(site.leftPos, 32);
    eq(site.rightPos, 36);
    eq(site.insert, PB.reverseComplement(group.newDna));
  });

  test("★ テールが挿入配列を復元する（設計全体が乗る等式）", function () {
    var context = example();
    ['D2A', 'ins@6:GGSGGS', 'del5-8', 'rep5-8:WW'].forEach(function (notation) {
      var result = onePair(context, notation);
      eq(PB.reverseComplement(result.pair.reverse.tail) + result.pair.forward.tail,
         result.site.insert, notation);
    });
  });

  test('アニーリング領域は鋳型と完全一致（3\' 末端のミスマッチなし）', function () {
    var context = example();
    var length = context.record.seq.length;
    ['D2A', 'ins@6:GGSGGS', 'del5-8'].forEach(function (notation) {
      var r = onePair(context, notation);
      var forward = context.record.subseq(r.site.rightPos,
        PB.wrapPosition(r.site.rightPos + r.pair.forward.anneal.length - 1, length));
      eq(r.pair.forward.anneal, forward, notation + ' の F');
      var reverse = PB.reverseComplement(context.record.subseq(
        PB.wrapPosition(r.site.leftPos - r.pair.reverse.anneal.length + 1, length),
        r.site.leftPos));
      eq(r.pair.reverse.anneal, reverse, notation + ' の R');
    });
  });

  test('産物長 = プラスミド長 - 除去 + 挿入', function () {
    var context = example();
    var r = onePair(context, 'ins@6:GGSGGS');
    eq(r.pair.productLength,
       context.record.seq.length - r.site.removedLength + r.site.insertLength);
  });

  test('プライマー長は上限内／アニール域は Tm 式の適用範囲内', function () {
    var context = example();
    var r = onePair(context, 'ins@6:GGSGGS');
    ok(r.pair.forward.length <= PB.DEFAULT_MAX_PRIMER_LEN, 'F の長さ');
    ok(r.pair.reverse.length <= PB.DEFAULT_MAX_PRIMER_LEN, 'R の長さ');
    ok(r.pair.forward.anneal.length <= PB.TM_FORMULA_MAX_LEN, 'F のアニール域');
    ok(r.pair.reverse.anneal.length <= PB.TM_FORMULA_MAX_LEN, 'R のアニール域');
  });

  test('候補はスコア順に並ぶ／この鋳型では結合部位が一意', function () {
    var context = example();
    var r = onePair(context, 'D2A');
    var scores = r.pairs.map(function (p) { return p.score; });
    eq(scores.slice().sort(function (a, b) { return a - b; }), scores);
    ok(r.pair.usable, '使用可');
    eq(r.pair.forward.bindingSites, 1);
    eq(r.pair.reverse.bindingSites, 1);
  });

  test('繰り返し配列の鋳型では非特異結合として検出される', function () {
    var unit = 'ATGCGGATTACCGTTAAGCCTTGACAGTCAGGTTACCTGA';
    var record = new PB.Record({ seq: new Array(11).join(unit), isCircular: true });
    var site = { leftPos: 100, rightPos: 104, removed: record.subseq(101, 103),
                 insert: 'GCG', removedLength: 3, insertLength: 3 };
    var pairs = PB.designPrimerPairs(record, site, {});
    ok(pairs.length > 0);
    eq(pairs[0].usable, false);
    ok(pairs[0].score > PB.PENALTY_NON_UNIQUE, 'スコアに失格分の加点');
  });

  test('直鎖の鋳型は拒否される（逆PCRは環状が前提）', function () {
    var record = new PB.Record({ seq: new Array(101).join('ACGT'),
                                 isCircular: false });
    throws(function () {
      PB.designPrimerPairs(record, { leftPos: 10, rightPos: 14, removed: 'ACG',
                                     insert: 'TTT', removedLength: 3,
                                     insertLength: 3 }, {});
    });
  });

  test('テールの分配は半分に近く、収まらなければ null', function () {
    eq(PB.splitPoint(40, 42, 42), 20);
    eq(PB.splitPoint(40, 10, 42), 30);
    eq(PB.splitPoint(100, 20, 20), null);
  });

  test('★ 検算: 4 種の編集すべてで目的のタンパク質に一致する', function () {
    var context = example();
    ['D2A', 'ins@6:GGSGGS', 'del5-8', 'rep5-8:WW', 'D2A,ins@6:GGSGGS']
      .forEach(function (notation) {
        var r = onePair(context, notation);
        var v = PB.verifyRound(context.record, context.feature, context.protein,
                               r.group, r.pair);
        eq(v.problems, [], notation);
        eq(v.observedProtein, v.expectedProtein, notation);
      });
  });

  test('検算: 長さの計算が合う', function () {
    var context = example();
    [['D2A', 0], ['ins@6:GGSGGS', 18], ['del5-8', -12], ['rep5-8:WW', -6]]
      .forEach(function (entry) {
        var r = onePair(context, entry[0]);
        var v = PB.verifyRound(context.record, context.feature, context.protein,
                               r.group, r.pair);
        eq(v.newRecord.length(), context.record.seq.length + entry[1], entry[0]);
      });
  });

  test('検算: 原点が保たれる／source は新しい全長に作り直される', function () {
    var context = example();
    var r = onePair(context, 'D2A');
    var v = PB.verifyRound(context.record, context.feature, context.protein,
                           r.group, r.pair);
    eq(v.newRecord.seq.slice(0, 50), context.record.seq.slice(0, 50));
    var inserted = onePair(context, 'ins@6:GGSGGS');
    var v2 = PB.verifyRound(context.record, context.feature, context.protein,
                            inserted.group, inserted.pair);
    var sources = v2.newRecord.features.filter(function (f) {
      return f.key.toLowerCase() === 'source';
    });
    eq(sources.length, 1);
    eq(sources[0].segments, [[1, v2.newRecord.length()]]);
  });

  test('★ 検算: テールを 1 塩基壊すと必ず検出される', function () {
    var context = example();
    var r = onePair(context, 'ins@6:GGSGGS');
    var broken = r.pair.forward.tail;
    var swapped = (broken[0] === 'A' ? 'C' : 'A') + broken.slice(1);
    var tampered = { site: r.pair.site, split: r.pair.split,
      forward: Object.assign({}, r.pair.forward,
        { tail: swapped, sequence: swapped + r.pair.forward.anneal }),
      reverse: r.pair.reverse };
    var v = PB.verifyRound(context.record, context.feature, context.protein,
                           r.group, tampered);
    eq(v.ok, false);
  });

  test('検算: PCR 産物の長さ', function () {
    var context = example();
    var r = onePair(context, 'ins@6:GGSGGS');
    eq(PB.pcrProduct(context.record, r.pair).length, r.pair.productLength);
  });

  test('検算: 逆鎖の CDS でも通る', function () {
    var original = PB.parseGenbank(EXAMPLE_GENBANK);
    var flipped = original.seq.slice(0, 100)
      + PB.reverseComplement(original.seq.slice(100, 700))
      + original.seq.slice(700);
    var record = new PB.Record({ name: 'pTEST-REV', seq: flipped,
                                 isCircular: true });
    var feature = new PB.Feature('CDS', [[101, 700]], -1);
    var protein = PB.translate(record.featureSequence(feature)).protein;
    ok(protein[0] === 'M' && protein.slice(-1) === '*', '逆鎖が ORF であること');
    eq(protein.slice(0, -1).split('*').length - 1, 0);
    var cds = record.featureSequence(feature);
    var edits = PB.parseEdits(protein[2] + '3W');
    PB.validateEdits(edits, protein);
    var group = PB.groupEdits(protein, edits, { cdsDna: cds })[0];
    var site = PB.editSite(record, feature, group);
    var pair = PB.designPrimerPairs(record, site, {})[0];
    var v = PB.verifyRound(record, feature, protein, group, pair);
    eq(v.problems, []);
    eq(v.observedProtein, v.expectedProtein);
  });

  test('診断用サイト: Thr-Gly は AgeI を同義コドンで書ける', function () {
    var preference = PB.getCodonPreference('ecoli').preference;
    var cds = 'ATG' + 'ACAGGA' + 'TAA';
    var placed = PB.placeSite(cds, 'ACCGGT', 3, 3, 8, preference);
    ok(placed !== null);
    eq(placed.cds, 'ATGACCGGTTAA');
    eq(PB.translate(placed.cds).protein, PB.translate(cds).protein);
    eq(placed.changed.length, 2);
  });

  test('診断用サイト: 区間の外の塩基は変えない', function () {
    var preference = PB.getCodonPreference('ecoli').preference;
    eq(PB.placeSite('ATGACAGGATAA', 'ACCGGT', 3, 3, 5, preference), null);
  });

  test('★ 診断用サイト: サイレントで、親になく、改変後に 1 箇所', function () {
    var context = example();
    var group = oneGroup(context, 'ins@6:GGSGGS');
    var site = PB.editSite(context.record, context.feature, group);
    var found = PB.findDiagnosticSites(context.record, context.feature, group,
                                       site, {});
    ok(found.length > 0, '候補が見つかること');
    var peptide = PB.translate(group.newDna).protein;
    found.forEach(function (candidate) {
      eq(PB.translate(candidate.newDna).protein, peptide, candidate.enzyme);
      eq(candidate.newDna.length, group.newDna.length);
      eq(PB.countOccurrences(context.record.seq, candidate.recognition, true), 0,
         candidate.enzyme + ' が親に無いこと');
      var rebuilt = context.record.subseq(site.rightPos, site.leftPos)
        + candidate.newDna;
      eq(PB.countOccurrences(rebuilt, candidate.recognition, true), 1,
         candidate.enzyme + ' が改変後に 1 箇所');
    });
  });

  test('診断用サイト: 挿入が無ければ候補は出ない', function () {
    var context = example();
    var group = oneGroup(context, 'del5-8');
    var site = PB.editSite(context.record, context.feature, group);
    eq(PB.findDiagnosticSites(context.record, context.feature, group, site, {}),
       []);
  });

  test('制限酵素表: 組み込みの認識配列はすべて回文', function () {
    Object.keys(PB.BUILTIN_ENZYMES).forEach(function (name) {
      var site = PB.BUILTIN_ENZYMES[name];
      eq(PB.reverseComplement(site), site, name);
    });
  });

  test('制限酵素表: ファイル読み込みと不正な塩基の拒否', function () {
    eq(PB.parseEnzymeTable('# custom\nMyI GGATCC\nYourI AAGCTT\n'),
       { MyI: 'GGATCC', YourI: 'AAGCTT' });
    throws(function () { PB.parseEnzymeTable('MyI GGNTCC\n'); });
  });

  test('コドン表: 大腸菌以外の宿主は断る', function () {
    // 2026-09-28: 測定値ファイルを渡す口は廃止した。組み込みの大腸菌だけ。
    throws(function () { PB.getCodonPreference('yeast'); }, 'ecoli のみ');
    eq(PB.getCodonPreference('ecoli').preference.L[0], 'CTG');
  });

  test('★ 全体: 設計から検算まで通り、診断用サイトを採用しても通る', function () {
    var context = example();
    var run = PB.designRun(context.record, context.feature, {
      mutations: 'D2A,ins@6:GGSGGS', constructName: 'pTEST_v2' });
    eq(run.ok, true);
    eq(run.rounds.length, 1);
    ok(run.rounds[0].diagnostics.length > 0, '診断用サイトの候補');
    var enzyme = run.rounds[0].diagnostics[0].enzyme;
    var adopted = PB.designRun(context.record, context.feature, {
      mutations: 'D2A,ins@6:GGSGGS', constructName: 'pTEST_v2',
      diagnosticEnzyme: enzyme });
    eq(adopted.ok, true);
    eq(adopted.rounds[0].adoptedDiagnostic.enzyme, enzyme);
    ok(adopted.rounds[0].site.insert !== run.rounds[0].site.insert,
       '挿入配列が書き換わること');
    eq(PB.countOccurrences(adopted.rounds[0].verification.newRecord.seq,
       adopted.rounds[0].adoptedDiagnostic.recognition, true), 1);
  });

  test('全体: 上限を超えるラウンドはプライマーを出さない', function () {
    var context = example();
    var run = PB.designRun(context.record, context.feature, {
      mutations: 'ins@6:' + new Array(41).join('W') });
    eq(run.ok, false);
    eq(run.rounds[0].pair, null);
  });

  test('全体: 知らない酵素名を指定したら報告する', function () {
    var context = example();
    var run = PB.designRun(context.record, context.feature, {
      mutations: 'D2A,ins@6:GGSGGS', diagnosticEnzyme: 'NotAnEnzyme' });
    eq(run.ok, false);
    ok(run.rounds[0].problems.join(' ').indexOf('NotAnEnzyme') !== -1);
  });

  test('全体: 目的タンパク質を渡しても差分から設計できる', function () {
    var context = example();
    var body = context.protein.replace(/\*+$/, '');
    var target = body.slice(0, 9) + 'GGSGGS' + body.slice(9);
    var run = PB.designRun(context.record, context.feature, {
      targetProtein: target });
    eq(run.ok, true);
    eq(run.newProtein.replace(/\*+$/, ''), target);
  });

  // ---- replacement mode ------------------------------------------------

  test('置換: 座標で指定した領域が正しく EditSite になる', function () {
    var context = example();
    var site = PB.replacementSite(context.record, 104, 115, 'GCG');
    eq(site.leftPos, 103);
    eq(site.rightPos, 116);
    eq(site.removed, context.record.subseq(104, 115));
    eq(site.removedLength, 12);
    eq(site.insert, 'GCG');
  });

  test('置換: 終了 = 開始 - 1 は「除去なしの挿入」', function () {
    var context = example();
    var site = PB.replacementSite(context.record, 118, 117, 'GGT');
    eq(site.removed, '');
    eq(site.removedLength, 0);
    eq(site.rightPos, 118);
    eq(site.leftPos, 117);
  });

  test('置換: 空の指定と不正な塩基は拒否する', function () {
    var context = example();
    throws(function () { PB.replacementSite(context.record, 118, 117, ''); },
           '除去も挿入もない');
    throws(function () { PB.replacementSite(context.record, 104, 115, 'GCX'); },
           '不正な塩基');
    throws(function () { PB.replacementSite(context.record, 0, 10, 'A'); },
           '範囲外');
  });

  test('置換: 貼り付けた配列から位置を特定する', function () {
    var context = example();
    var located = PB.locateSelection(context.record,
                                     context.record.subseq(300, 331));
    eq([located.start, located.end, located.length], [300, 331, 32]);
    eq(located.onReverseStrand, false);
  });

  test('置換: 逆鎖に一致した場合も順方向の座標で返す', function () {
    var context = example();
    var located = PB.locateSelection(context.record,
      PB.reverseComplement(context.record.subseq(500, 531)));
    eq([located.start, located.end], [500, 531]);
    eq(located.onReverseStrand, true);
  });

  test('置換: 空白や数字が混じっていても位置を特定できる', function () {
    var context = example();
    var picked = context.record.subseq(300, 331);
    var messy = '  1 ' + picked.slice(0, 10) + '\n 11 ' + picked.slice(10) + ' ';
    eq(PB.locateSelection(context.record, messy).start, 300);
  });

  test('置換: 一意でない選択・見つからない選択は拒否する', function () {
    var context = example();
    throws(function () { PB.locateSelection(context.record, 'ACGT'); }, '複数一致');
    throws(function () {
      PB.locateSelection(context.record, 'TTTTTTTTTTTTTTTTTTTTTTTTGGGG');
    }, '不一致');
    throws(function () { PB.locateSelection(context.record, ''); }, '空');
  });

  test('★ 置換: 4 通りの操作すべてで、期待どおりの環状配列になる', function () {
    var context = example();
    var record = context.record;
    [['置換', { start: 104, end: 115, newSequence: 'GCGGAACTGGCGACC' }],
     ['欠失', { start: 104, end: 130, newSequence: '' }],
     ['挿入', { start: 118, end: 117, newSequence: 'GGTGGCAGCGGTGGCAGC' }],
     ['同長置換', { start: 104, end: 106, newSequence: 'GCG' }]
    ].forEach(function (entry) {
      var result = PB.designReplacement(record, entry[1]);
      eq(result.problems, [], entry[0]);
      eq(result.verification.problems, [], entry[0] + ' の検算');
      // 期待する環状配列を、プライマーとは独立に組み立てて比較する
      var expected = record.subseq(result.site.rightPos, result.site.leftPos)
        + result.site.insert;
      eq(result.verification.expectedSequence, expected, entry[0]);
      eq(result.verification.newRecord.length(),
         record.seq.length - result.site.removedLength
         + result.site.insertLength, entry[0] + ' の長さ');
    });
  });

  test('★ 置換: テールが新しい配列を復元する', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 115, newSequence: 'GCGGAACTGGCGACC' });
    eq(PB.reverseComplement(result.pair.reverse.tail) + result.pair.forward.tail,
       result.site.insert);
  });

  test('★ 置換: テールを 1 塩基壊すと検算が落ちる', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 115, newSequence: 'GCGGAACTGGCGACC' });
    var tail = result.pair.forward.tail;
    var swapped = (tail[0] === 'A' ? 'C' : 'A') + tail.slice(1);
    var tampered = { site: result.pair.site, split: result.pair.split,
      forward: Object.assign({}, result.pair.forward,
        { tail: swapped, sequence: swapped + result.pair.forward.anneal }),
      reverse: result.pair.reverse };
    eq(PB.verifyReplacement(context.record, result.site, tampered).ok, false);
  });

  test('置換: 選択範囲の外は 1 塩基も変わらない', function () {
    var context = example();
    var record = context.record;
    var result = PB.designReplacement(record,
      { start: 104, end: 115, newSequence: 'GCGGAACTGGCGACC' });
    var after = result.verification.newRecord;
    // 編集前の 1..103 は、新しい記録でもそのまま先頭にある
    eq(after.subseq(1, 103), record.subseq(1, 103));
    // 編集後の残りも、元の 116 以降と一致する
    eq(after.subseq(104 + 15, after.length()),
       record.subseq(116, record.seq.length));
  });

  test('置換: アミノ酸で指定すると翻訳が一致する DNA になる', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 115, newSequence: 'AELATGGS', asProtein: true });
    eq(result.problems, []);
    eq(PB.translate(result.site.insert).protein, 'AELATGGS');
    eq(result.site.insertLength, 24);
    ok(result.notes.join(' ').indexOf('アミノ酸 8 残基') !== -1, '変換の注記');
  });

  test('置換: アミノ酸指定で回避配列を作らない', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 115, newSequence: 'EFAAKK', asProtein: true,
        avoid: ['GAATTC'] });
    eq(result.problems, []);
    ok(result.site.insert.indexOf('GAATTC') === -1, 'EcoRI を作らない');
    eq(PB.translate(result.site.insert).protein, 'EFAAKK');
  });

  test('置換: 読み枠に合わない選択は注意を出す', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 114, newSequence: 'AA', asProtein: true });
    ok(result.notes.join(' ').indexOf('3 の倍数ではありません') !== -1);
  });

  test('置換: 長すぎる新配列は上限として報告する', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 106, newSequence: new Array(101).join('A') });
    eq(result.ok, false);
    eq(result.pair, null);
    ok(result.problems.join(' ').indexOf('上限') !== -1);
  });

  test('置換: プライマーは上限内、アニール域は Tm 式の範囲内', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 115, newSequence: 'GCGGAACTGGCGACC' });
    ok(result.pair.forward.length <= PB.DEFAULT_MAX_PRIMER_LEN);
    ok(result.pair.reverse.length <= PB.DEFAULT_MAX_PRIMER_LEN);
    ok(result.pair.forward.anneal.length <= PB.TM_FORMULA_MAX_LEN);
    ok(result.pair.reverse.anneal.length <= PB.TM_FORMULA_MAX_LEN);
  });


  // ---- reading frame ---------------------------------------------------

  test('★ 読み枠: コドンの境界からずれた挿入を拒否する', function () {
    var context = example();
    var cds = context.record.cdsFeatures()[0];          // 101..700
    // 104 は残基 2 の先頭。105 はその 2 塩基目。
    var ok = PB.designReplacement(context.record,
      { start: 104, end: 103, newSequence: 'GGS', asProtein: true, cds: cds });
    eq(PB.translate(ok.site.insert).protein, 'GGS');
    throws(function () {
      PB.designReplacement(context.record,
        { start: 105, end: 104, newSequence: 'GGS', asProtein: true, cds: cds });
    }, 'コドンの境界からずれた挿入');
  });

  test('読み枠: CDS の外は拒否／3 の倍数でない範囲も拒否', function () {
    var context = example();
    var cds = context.record.cdsFeatures()[0];
    throws(function () {
      PB.designReplacement(context.record,
        { start: 50, end: 49, newSequence: 'GGS', asProtein: true, cds: cds });
    }, 'CDS の外');
    throws(function () {
      PB.designReplacement(context.record,
        { start: 104, end: 114, newSequence: 'GG', asProtein: true, cds: cds });
    }, '3 の倍数でない範囲');
  });

  test('読み枠: CDS を渡さなければ確認できない旨を伝える', function () {
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 105, end: 104, newSequence: 'GGS', asProtein: true });
    // 拒否はできないが、黙って通してもいけない
    eq(result.problems, []);
    ok(result.notes.join(' ').indexOf('確認できませんでした') !== -1,
       '確認できない旨の注意');
  });

  test('読み枠: 逆鎖の CDS では高い側の端でコドン境界を見る', function () {
    var context = example();
    var minus = new PB.Feature('CDS', [[101, 700]], -1, {});
    var result = PB.designReplacement(context.record,
      { start: 104, end: 115, newSequence: 'AELA', asProtein: true,
        cds: minus });
    eq(result.problems, []);
    ok(result.notes.join(' ').indexOf('読み枠を確認しました') !== -1);
  });

  test('★ 読み枠: アミノ酸での挿入が実際に動く（以前は落ちていた）', function () {
    // 文脈取得のための仮サイトを、空の挿入として組み立てていたため、
    // アミノ酸 + 挿入 の組み合わせが例外になっていた。
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 103, newSequence: 'GGSGGS', asProtein: true });
    eq(result.problems, []);
    eq(PB.translate(result.site.insert).protein, 'GGSGGS');
    eq(result.site.removedLength, 0);
    eq(result.verification.problems, []);
  });

  test('回帰: 報告された入力（104 の直前に WW を挿入）', function () {
    // 2026-08-28 報告: asProtein + 挿入（end = start - 1）が
    // 「何も除去せず、何も挿入しない指定です」で落ちていた。
    // 報告された入力をそのまま残して固定する。
    var context = example();
    var result = PB.designReplacement(context.record,
      { start: 104, end: 103, newSequence: 'WW', asProtein: true });
    eq(result.problems, []);
    eq(result.site.insertLength, 6);
    eq(PB.translate(result.site.insert).protein, 'WW');
    eq(result.verification.ok, true);
  });

  test('spanBounds: 挿入でも例外にならず、範囲だけを返す', function () {
    var context = example();
    var empty = PB.spanBounds(context.record, 118, 117);
    eq(empty.empty, true);
    eq(empty.removedLength, 0);
    eq([empty.leftPos, empty.rightPos], [117, 118]);
    var span = PB.spanBounds(context.record, 104, 115);
    eq(span.empty, false);
    eq(span.removedLength, 12);
    eq([span.leftPos, span.rightPos], [103, 116]);
  });

  /* Reading direction, fixed 2026-09-11.
   *
   * W has exactly one codon, TGG. A CDS on the minus strand reads the plus
   * strand backwards, so the plus strand has to carry CCA CCA for that CDS to
   * read WW. Before the fix the tool inserted TGGTGG whatever the strand,
   * which a minus-strand CDS reads as PP. */

  test('逆鎖 CDS: アミノ酸は読み枠の向きでコドンにする', function () {
    var record = PB.parseGenbank(EXAMPLE_GENBANK);
    var cds = new PB.Feature('CDS', [[101, 700]], -1, {});
    var map = PB.cdsPositionMap(cds);
    // residue 5 is 688,687,686 and residue 6 is 685,684,683
    var span = PB.residueToPlasmid(map, 5).concat(PB.residueToPlasmid(map, 6));
    var result = PB.designReplacement(record, {
      start: Math.min.apply(null, span), end: Math.max.apply(null, span),
      newSequence: 'WW', asProtein: true, cds: cds });
    eq(result.site.insert, 'CCACCA');
    eq(PB.translate(result.verification.newRecord.featureSequence(cds))
       .protein.slice(4, 6), 'WW');
  });

  test('順鎖 CDS: これまでどおり前向きにコドンにする', function () {
    var record = PB.parseGenbank(EXAMPLE_GENBANK);
    var cds = new PB.Feature('CDS', [[101, 700]], 1, {});
    var map = PB.cdsPositionMap(cds);
    var span = PB.residueToPlasmid(map, 5).concat(PB.residueToPlasmid(map, 6));
    var result = PB.designReplacement(record, {
      start: Math.min.apply(null, span), end: Math.max.apply(null, span),
      newSequence: 'WW', asProtein: true, cds: cds });
    eq(result.site.insert, 'TGGTGG');
    eq(PB.translate(result.verification.newRecord.featureSequence(cds))
       .protein.slice(4, 6), 'WW');
  });

  test('逆鎖 CDS: 挿入も読み枠の向きになる', function () {
    var record = PB.parseGenbank(EXAMPLE_GENBANK);
    var cds = new PB.Feature('CDS', [[101, 700]], -1, {});
    var map = PB.cdsPositionMap(cds);
    var start = Math.min.apply(null, PB.residueToPlasmid(map, 5));
    var result = PB.designReplacement(record, {
      start: start, end: start - 1, newSequence: 'WW', asProtein: true,
      cds: cds });
    eq(result.site.removedLength, 0);
    eq(result.site.insert, 'CCACCA');
    var grown = new PB.Feature('CDS', [[101, 706]], -1, {});
    ok(PB.translate(result.verification.newRecord.featureSequence(grown))
       .protein.indexOf('WW') !== -1);
  });
})();
