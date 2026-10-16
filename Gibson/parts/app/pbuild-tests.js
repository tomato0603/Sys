
/* Tests for reading a Protein Builder backup. The file below is written out
 * here rather than loaded, so the page can test itself with no file and no
 * network. Every sequence in it is made up. */
(function () {
  'use strict';
  var test = PBTESTS.test, eq = PBTESTS.eq, ok = PBTESTS.ok;

  // A protein to search in. HeadPiece sits at residues 11..15.
  var PROTEIN = 'MKVLATGDNS' + 'ACDEF' + 'GHIK' + 'LMNPQ' + 'STVWYRRQ*';
  //             1........10   11..15   16.19   20..24   25...

  var BACKUP = JSON.stringify({
    format: 'protein-domain-builder',
    version: 2,
    library: [
      { id: 'lib-a', name: 'HeadPiece', seq: 'acdef' },
      { id: 'lib-b', name: 'TailPiece', seq: 'LMN PQ' },
      { id: 'lib-gone', name: 'ElsewhereOnly', seq: 'WWWWWWWW' },
      { id: 'lib-empty', name: 'NoSequence', seq: '' }
    ],
    // The builds are deliberately not offered; their presence must not
    // change what comes out.
    constructs: [{ id: 'c1', name: 'BUILD-1', elements: [
      { type: 'block', libId: 'lib-a', start: 1, end: 5 }] }]
  });

  test('P-Build: JSON でないものは理由を言って断る', function () {
    var message = '';
    try { PBUILD.parseBackup('<html>'); }
    catch (error) { message = error.message; }
    ok(message.indexOf('JSON') !== -1, message);
  });

  test('P-Build: format が違うファイルは断る', function () {
    var message = '';
    try { PBUILD.parseBackup('{"format":"something-else"}'); }
    catch (error) { message = error.message; }
    ok(message.indexOf(PBUILD.FORMAT) !== -1, message);
  });

  test('P-Build: ライブラリの配列を名前つきで取り出す', function () {
    eq(PBUILD.parseBackup(BACKUP).sequences.map(function (s) {
      return s.name + ' ' + s.sequence;
    }), ['HeadPiece ACDEF', 'TailPiece LMNPQ', 'ElsewhereOnly WWWWWWWW']);
  });

  test('P-Build: 配列が空の項目は飛ばして、その旨を残す', function () {
    var parsed = PBUILD.parseBackup(BACKUP);
    eq(parsed.sequences.length, 3);
    ok(parsed.warnings.join(' ').indexOf('1 件') !== -1,
       parsed.warnings.join(' '));
  });

  test('P-Build: 配列が 1 つも無ければ断る', function () {
    var message = '';
    try {
      PBUILD.parseBackup(JSON.stringify({
        format: PBUILD.FORMAT, library: [], constructs: [] }));
    } catch (error) { message = error.message; }
    ok(message.indexOf('タンパク質配列') !== -1, message);
  });

  test('P-Build: 1 箇所だけにある配列は残基番号を返す', function () {
    eq(PBUILD.findIn(PROTEIN, 'ACDEF'), { count: 1, resStart: 11, resEnd: 15 });
  });

  test('P-Build: 無い配列は count 0 で、位置を返さない', function () {
    eq(PBUILD.findIn(PROTEIN, 'WWWWWWWW'), { count: 0 });
  });

  test('P-Build: 2 箇所にある配列は位置を決めない', function () {
    eq(PBUILD.findIn('ACDEFQQACDEF', 'ACDEF'), { count: 2 });
  });

  test('P-Build: 空の配列を探しには行かない', function () {
    eq(PBUILD.findIn(PROTEIN, ''), { count: 0 });
  });

  test('P-Build: 小文字と空白は取り除いて比べる', function () {
    eq(PBUILD.cleanSeq(' ac d\nef1 '), 'ACDEF');
  });

  test('P-Build: 同じ名前と配列の組は二重に貯めない', function () {
    var merged = PBUILD.merge(
      [{ name: 'A', sequence: 'ACDEF' }],
      [{ name: 'A', sequence: 'ACDEF' }, { name: 'B', sequence: 'LMNPQ' }]);
    eq(merged.added, 1);
    eq(merged.list.map(function (e) { return e.name; }), ['A', 'B']);
  });

  test('P-Build: 名前が違えば同じ配列でも別の項目として貯める', function () {
    var merged = PBUILD.merge(
      [{ name: 'A', sequence: 'ACDEF' }],
      [{ name: 'AnotherName', sequence: 'ACDEF' }]);
    eq(merged.added, 1);
    eq(merged.list.length, 2);
  });

  test('P-Build: 一致したものを上に、プラスミドに並ぶ順で出す', function () {
    var rows = PBUILD.arrange(PROTEIN, [
      { name: 'Missing', sequence: 'WWWWWWWW' },   // なし
      { name: 'Tail', sequence: 'LMNPQ' },         // 残基 20..24
      { name: 'Head', sequence: 'ACDEF' }          // 残基 11..15
    ]);
    eq(rows.map(function (r) { return r.entry.name; }),
       ['Head', 'Tail', 'Missing']);
    // 元の位置は残す。押したときにどれを指すか決まらなくなるため
    eq(rows.map(function (r) { return r.index; }), [2, 1, 0]);
  });

  test('P-Build: 決められないものは、見つからないものより上に置く',
       function () {
    var rows = PBUILD.arrange('AAQQAAZZ', [
      { name: 'Missing', sequence: 'WWWW' },
      { name: 'Twice', sequence: 'AA' },
      { name: 'Once', sequence: 'ZZ' }
    ]);
    eq(rows.map(function (r) { return r.entry.name; }),
       ['Once', 'Twice', 'Missing']);
  });

  test('P-Build: 読み枠が無ければ、どれも一致なしとして並べる', function () {
    var rows = PBUILD.arrange('', [{ name: 'Head', sequence: 'ACDEF' }]);
    eq(rows[0].found.count, 0);
  });
})();
