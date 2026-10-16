
/* Tests for the sequence view: T7 recognition, the display window and the
 * ORF between promoter and terminator. The plasmid below is generated, with
 * the standard T7 promoter and terminator spliced into the synthetic example.
 */
(function () {
  'use strict';
  var test = PBTESTS.test, eq = PBTESTS.eq, ok = PBTESTS.ok;
  var T7_EXAMPLE = "LOCUS       pT7-EXAMPLE             966 bp    DNA     circular SYN 27-AUG-2026\nDEFINITION  synthetic example with a T7 promoter and terminator. Not a real construct.\nCOMMENT     Generated for testing; every 18-mer of the parent was unique.\nFEATURES             Location/Qualifiers\n     source          1..966\n                     /organism=\"synthetic DNA construct\"\n     promoter        101..118\n                     /label=\"T7 promoter\"\n     terminator      789..836\n                     /label=\"T7 terminator\"\nORIGIN\n        1 accctcagac gctactaaag cacgactgct tagataagtc catcatcgac actgggggcg\n       61 gcatggtggg ccgacttcga gctttcgcct cgctgtattg taatacgact cactatagct\n      121 ctagcccatg gcgattaatt taaccgtttc ctttccgtat ggacattccc ttaccccttc\n      181 gtataaggga agttctagct gatctaggct taaatgttgt acgatcagta atgtttgggg\n      241 atagcgccac ggcattggct aaagctggct gccacgtggg acaggtaacc tcccgtgctc\n      301 tccgctcgca agacgcgcct gcccagtcgg ctcgttccgg tggcaatcga agtagatttg\n      361 caccaataca gaaacatgat ctctcatttg ccaggccgca gtactgcccg agaccgcgaa\n      421 cgaccgccgc gtggtacgtc gacgaaaagt cctggccctg cctgaaatac gggttgttca\n      481 agggtgaatt cctccacatc tctctatttg ccagtccatg cgatcgggcg atacccctac\n      541 ctcctcgaaa ttcgtctatc cccacccgga taatatatct gcaacgtctg ctttatctat\n      601 acgaggatgt gcgaacctat ctatgctcat taccaccgca catacaactg tccgtgttca\n      661 ttgcgcatgt acgtttctat cccgaatgta tttcaaacag cagagccgaa ccccgaagcg\n      721 acttcgcgca gactcacccc tgcgaaattg catcataagt cataaactag ctcgcttccg\n      781 ctgcaagcct agcataaccc cttggggcct ctaaacgggt cttgaggggt tttttgaggt\n      841 tgtgggtctc tacagcttcc aaactttgag agctcttttc cgattacagc ggctttaata\n      901 ctggaaattt ttccgtgatg ttgttcgaat ataatcacgt tcaatcaata actatttctt\n      961 ttcaag\n//\n";
  PBTESTS.T7_EXAMPLE = T7_EXAMPLE;

  var PROMOTER = [101, 118], TERMINATOR = [789, 836], CDS = [159, 758];

  function loaded() { return PB.parseGenbank(T7_EXAMPLE); }

  test('T7: promoter と terminator を見つける', function () {
    var record = loaded();
    var promoter = GEVIEW.findMotif(record, GEVIEW.T7_PROMOTER);
    var terminator = GEVIEW.findMotif(record, GEVIEW.T7_TERMINATOR);
    eq([promoter.start, promoter.end], PROMOTER);
    eq([terminator.start, terminator.end], TERMINATOR);
    eq(promoter.strand, 1);
    eq(promoter.hits, 1);
  });

  test('T7: 見つからない配列は null を返す', function () {
    eq(GEVIEW.findMotif(loaded(), 'GGGGGGGGGGGGGGGGGGGGGGGG'), null);
  });

  test('T7: 表示窓は promoter の N 前から terminator の N 後まで', function () {
    var record = loaded();
    // 前後 40 塩基なら窓はプラスミドより短く収まる
    var w = GEVIEW.expressionWindow(record, null, null, 40);
    eq(w.found, true);
    eq(w.wholePlasmid, false);
    eq(w.from, PROMOTER[0] - 40);
    eq(w.to, TERMINATOR[1] + 40);
    var positions = GEVIEW.windowPositions(record, w.from, w.to);
    eq(positions.length, TERMINATOR[1] + 40 - (PROMOTER[0] - 40) + 1);
    eq(positions[0], w.from);
    eq(positions[positions.length - 1], w.to);
    // 窓の中に promoter と terminator が入っていること
    ok(positions.indexOf(PROMOTER[0]) !== -1, 'promoter が窓の中');
    ok(positions.indexOf(TERMINATOR[1]) !== -1, 'terminator が窓の中');
  });

  test('T7: 窓がプラスミドより長くなる場合は全長を表示する', function () {
    var record = loaded();
    // 前後 120 塩基だと 976 塩基必要で、966 bp のプラスミドに収まらない
    var w = GEVIEW.expressionWindow(record);
    eq(w.wholePlasmid, true);
    eq([w.from, w.to], [1, record.length()]);
    var positions = GEVIEW.windowPositions(record, w.from, w.to);
    eq(positions.length, record.length());
    ok(positions.indexOf(PROMOTER[0]) !== -1, 'promoter が窓の中');
    ok(positions.indexOf(TERMINATOR[1]) !== -1, 'terminator が窓の中');
  });

  test('T7: promoter と terminator に挟まれた領域', function () {
    var w = GEVIEW.expressionWindow(loaded());
    eq(w.between.start, PROMOTER[1] + 1);
    eq(w.between.end, TERMINATOR[0] - 1);
  });

  test('T7: 挟まれた領域から ORF を見つけ、元の CDS と一致する', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record);
    var orf = GEVIEW.findOrf(record, w.between, w.strand);
    ok(orf !== null, 'ORF が見つかること');
    eq(orf.segments, [CDS]);
    var protein = PB.translate(record.featureSequence(orf)).protein;
    eq(protein[0], 'M');
    eq(protein.slice(-1), '*');
    eq(protein.slice(0, -1).split('*').length - 1, 0);
  });

  test('T7: 表示窓が原点をまたぐ場合も長さが合う', function () {
    var record = loaded();
    var n = record.length();
    var positions = GEVIEW.windowPositions(record, n - 10, 10);
    eq(positions.length, 21);
    eq(positions[0], n - 10);
    eq(positions[10], n);
    eq(positions[11], 1);
    eq(positions[20], 10);
  });

  test('表示: 塩基ごとに data-pos が付き、位置と文字が対応する', function () {
    var record = loaded();
    var html = GEVIEW.renderDna(record, GEVIEW.windowPositions(record, 1, 60), []);
    var matches = html.match(/data-pos="(\d+)">([ACGT])</g) || [];
    eq(matches.length, 60);
    var first = /data-pos="1">([ACGT])</.exec(html);
    eq(first[1], record.seq[0]);
    var last = /data-pos="60">([ACGT])</.exec(html);
    eq(last[1], record.seq[59]);
  });

  test('表示: 残基ごとに data-res と data-codon が付く', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record);
    var orf = GEVIEW.findOrf(record, w.between, w.strand);
    var rendered = GEVIEW.renderProtein(record, orf);
    var first = /data-res="1" data-codon="([\d,]+)"/.exec(rendered.html);
    eq(first[1].split(',').map(Number), PB.residueToPlasmid(rendered.map, 1));
    var third = /data-res="3" data-codon="([\d,]+)"/.exec(rendered.html);
    eq(third[1].split(',').map(Number), PB.residueToPlasmid(rendered.map, 3));
  });

  test('表示: promoter と terminator に色のクラスが付く', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record);
    var html = GEVIEW.renderDna(record,
      GEVIEW.windowPositions(record, w.from, w.to),
      [{ start: w.promoter.start, end: w.promoter.end, cls: 'prom' },
       { start: w.terminator.start, end: w.terminator.end, cls: 'term' }]);
    eq((html.match(/class="b prom"/g) || []).length, GEVIEW.T7_PROMOTER.length);
    eq((html.match(/class="b term"/g) || []).length, GEVIEW.T7_TERMINATOR.length);
  });

  test('置換モード: 選択範囲と新配列からプライマーが出る（別エージェントのエンジン）', function () {
    var record = loaded();
    var result = PB.designReplacement(record,
      { start: 162, end: 164, newSequence: 'GCG' });
    eq(result.verification.problems, []);
    eq(result.site.removedLength, 3);
    eq(result.site.insertLength, 3);
  });

  test('読み枠: ORF を cds として渡すと「確認しました」になる', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record);
    var orf = GEVIEW.findOrf(record, w.between, w.strand);
    var map = PB.cdsPositionMap(orf);
    var first = PB.residueToPlasmid(map, 5), last = PB.residueToPlasmid(map, 8);
    var start = Math.min.apply(null, first.concat(last));
    var end = Math.max.apply(null, first.concat(last));

    // cds なしだと、確認できなかったという注記が出る
    var without = PB.designReplacement(record,
      { start: start, end: end, newSequence: 'WWWW', asProtein: true });
    ok(without.notes.join(' ').indexOf('確認できませんでした') !== -1,
       'cds なしでは確認できない旨の注記');

    // Feature を渡すと確認済みになる
    var withFeature = PB.designReplacement(record,
      { start: start, end: end, newSequence: 'WWWW', asProtein: true, cds: orf });
    ok(withFeature.notes.join(' ').indexOf('読み枠を確認しました') !== -1,
       'Feature を渡すと確認済み');
    ok(withFeature.notes.join(' ').indexOf('確認できませんでした') === -1,
       '確認できない旨の注記が消えている');

    // 位置の文字列でも同じ
    var withText = PB.designReplacement(record,
      { start: start, end: end, newSequence: 'WWWW', asProtein: true,
        cds: orf.locationText() });
    ok(withText.notes.join(' ').indexOf('読み枠を確認しました') !== -1,
       '文字列でも確認済み');
  });

  test('読み枠: アミノ酸での挿入のみも cds があれば通る', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record);
    var orf = GEVIEW.findOrf(record, w.between, w.strand);
    var start = PB.residueToPlasmid(PB.cdsPositionMap(orf), 5)[0];
    var result = PB.designReplacement(record,
      { start: start, end: start - 1, newSequence: 'WW', asProtein: true,
        cds: orf });
    eq(result.site.removedLength, 0);
    eq(result.site.insertLength, 6);
    eq(PB.translate(result.site.insert).protein, 'WW');
    ok(result.notes.join(' ').indexOf('読み枠を確認しました') !== -1);
  });

  test('読み枠: コドンの境界でない位置は拒否される', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record);
    var orf = GEVIEW.findOrf(record, w.between, w.strand);
    var start = PB.residueToPlasmid(PB.cdsPositionMap(orf), 5)[0];
    PBTESTS.throws(function () {
      PB.designReplacement(record,
        { start: start + 1, end: start + 6, newSequence: 'WW',
          asProtein: true, cds: orf });
    }, 'コドン境界でない開始位置');
  });

  test('置換モード: アミノ酸で指定してもコドンに変換される', function () {
    var record = loaded();
    var result = PB.designReplacement(record,
      { start: 162, end: 167, newSequence: 'WW', asProtein: true });
    eq(result.site.insertLength, 6);
    eq(PB.translate(result.site.insert).protein, 'WW');
    eq(result.verification.problems, []);
  });

  /* Which ORF is shown decides every residue number the mutation-notation
   * mode uses, so the three ways of getting it wrong are pinned here. */

  test('ORF: ファイルの CDS 注釈があればそれを使う', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record, GEVIEW.T7_PROMOTER,
                                    GEVIEW.T7_TERMINATOR);
    // The guessed ORF for this plasmid runs 159..758.
    eq(GEVIEW.findOrf(record, w.between, w.strand).segments, [CDS]);
    // Declaring a shorter CDS must win over the longest-ORF guess.
    record.features.push(new PB.Feature('CDS', [[201, 500]], 1,
                                        { label: 'declared' }));
    var chosen = GEVIEW.findOrf(record, w.between, w.strand);
    eq(chosen.segments, [[201, 500]]);
    eq(chosen.label(), 'declared');
  });

  test('ORF: 同じ終止コドンを共有する開始候補を長い順に返す', function () {
    var record = loaded();
    var w = GEVIEW.expressionWindow(record, GEVIEW.T7_PROMOTER,
                                    GEVIEW.T7_TERMINATOR);
    var candidates = GEVIEW.orfCandidates(record, w.between, w.strand);
    ok(candidates.length > 1);
    eq(candidates[0].segments, [CDS]);
    for (var i = 1; i < candidates.length; i++) {
      ok(candidates[i - 1].residues >= candidates[i].residues);
    }
  });

  /* A plasmid whose only ORF stops inside the terminator. Built here rather
   * than pasted, so the numbers below can be counted by hand: 60 bases of
   * filler, the 18-base promoter, ATG, 20 GCT codons, then the terminator,
   * whose first codons read CTA GCA TAA -- the stop. */
  function stopPastTerminator() {
    var filler = '';
    while (filler.length < 60) { filler += 'CTT'; }
    var body = 'ATG';
    for (var i = 0; i < 20; i++) { body += 'GCT'; }
    var seq = filler + GEVIEW.T7_PROMOTER + body + GEVIEW.T7_TERMINATOR + filler;
    return new PB.Record({ name: 'stop-past-term', seq: seq, isCircular: true,
                           features: [], sourceFormat: 'test' });
  }

  test('ORF: 終止コドンが terminator の内側にあっても見つける', function () {
    var record = stopPastTerminator();
    var w = GEVIEW.expressionWindow(record, GEVIEW.T7_PROMOTER,
                                    GEVIEW.T7_TERMINATOR);
    ok(w.found);
    var orf = GEVIEW.findOrf(record, w.between, w.strand);
    ok(!!orf);
    // M + 20 A + L + A, then the stop. Counted from the construction above.
    eq(PB.translate(record.featureSequence(orf)).protein, 'M' +
       'AAAAAAAAAAAAAAAAAAAA' + 'LA*');
    ok(orf.beyondTerminator);
  });
})();
