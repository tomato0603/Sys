/* ===================================================================
 * Reading a Protein Builder backup.
 *
 * Protein Builder (P-Build/protein_builder.html) is the sibling tool that
 * holds the user's protein sequences. They live in that page's localStorage,
 * not inside its HTML, so the file that comes over here is the one its
 * "SYSTEM SAVE (.json)" button writes.
 *
 * Only the library is read -- the named protein sequences themselves. That
 * file also carries the builds assembled out of them, plus requests and
 * colour rules; on 2026-09-27 the user said the builds are not what they
 * want here, so they are not offered. (They sit under `constructs`, each
 * with an `elements` list naming library ids and slices, if that changes.)
 * =================================================================== */
var PBUILD = (function () {
  'use strict';

  var FORMAT = 'protein-domain-builder';

  // Protein Builder's own cleanSeq: letters only, upper case.
  function cleanSeq(text) {
    return String(text === undefined || text === null ? '' : text)
      .toUpperCase().replace(/[^A-Z]/g, '');
  }

  function parseBackup(text) {
    var data;
    try { data = JSON.parse(text); }
    catch (error) {
      throw new Error('JSON として読めませんでした。Protein Builder の'
        + '「SYSTEM SAVE (.json)」で保存したファイルを選んでください。');
    }
    if (!data || data.format !== FORMAT) {
      throw new Error('Protein Builder の保存ファイルではないようです'
        + '（format が ' + FORMAT + ' ではありません）。');
    }
    var sequences = [], empty = 0;
    (data.library || []).forEach(function (item) {
      if (!item) { return; }
      var sequence = cleanSeq(item.seq);
      if (!sequence) { empty++; return; }
      sequences.push({ name: String(item.name || '名前なし'),
                       sequence: sequence });
    });
    var warnings = [];
    if (empty) {
      warnings.push('配列が空の項目 ' + empty + ' 件は飛ばしました。');
    }
    if (!sequences.length) {
      throw new Error('タンパク質配列が入っていませんでした。');
    }
    return { sequences: sequences, warnings: warnings };
  }

  /* Where a sequence sits in a translated CDS, as residue numbers.
   *
   * count is what the caller acts on: 0 means it is not in this reading
   * frame, and more than 1 means the position cannot be decided. Neither is
   * guessed at, because a label put on the wrong repeat is worse than no
   * label at all. */
  function findIn(protein, sequence) {
    if (!sequence) { return { count: 0 }; }
    var count = 0, at = -1, from = 0, next;
    while ((next = protein.indexOf(sequence, from)) !== -1) {
      if (at === -1) { at = next; }
      count++;
      from = next + 1;
    }
    return count === 1
      ? { count: 1, resStart: at + 1, resEnd: at + sequence.length }
      : { count: count };
  }

  /* Add newly read sequences to the ones already kept, without duplicating.
   *
   * A backup is read again whenever the user exports a new one, and most of
   * it will be what is already here. Name and sequence together are the
   * identity: the same peptide under two names is two entries, because the
   * name is what ends up on the plasmid. Protein Builder's own ids are not
   * used -- they do not survive its merges. */
  function merge(stored, incoming) {
    var list = (stored || []).slice(), seen = {}, added = 0;
    list.forEach(function (entry) {
      seen[entry.name + '\t' + entry.sequence] = true;
    });
    (incoming || []).forEach(function (entry) {
      var key = entry.name + '\t' + entry.sequence;
      if (seen[key]) { return; }
      seen[key] = true;
      list.push({ name: entry.name, sequence: entry.sequence });
      added++;
    });
    return { list: list, added: added };
  }

  /* The order to show them in: the ones that are on this plasmid first, in
   * the order they sit along it, then the ambiguous ones, then the rest.
   *
   * A library grows, and once it is kept between visits the rows worth
   * acting on are the few that match. Leaving them in file order pushes
   * them off the screen (asked for 2026-09-27). */
  function arrange(protein, sequences) {
    var rows = (sequences || []).map(function (entry, index) {
      return { entry: entry, index: index,
               found: protein ? findIn(protein, entry.sequence)
                              : { count: 0 } };
    });
    function rank(row) {
      return row.found.count === 1 ? 0 : row.found.count > 1 ? 1 : 2;
    }
    rows.sort(function (a, b) {
      if (rank(a) !== rank(b)) { return rank(a) - rank(b); }
      if (rank(a) === 0) { return a.found.resStart - b.found.resStart; }
      return a.index - b.index;
    });
    return rows;
  }

  return { FORMAT: FORMAT, cleanSeq: cleanSeq, parseBackup: parseBackup,
           findIn: findIn, merge: merge, arrange: arrange };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = PBUILD; }
