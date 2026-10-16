
/* ===================================================================
 * Gene Editor: the sequence view and range selection.
 *
 * Everything the user does starts here: they pick a region by looking at
 * the sequence, not by typing coordinates. Positions are 1-based, matching
 * ApE, so what is on screen matches what the sequence editor shows.
 * =================================================================== */
var GEVIEW = (function () {
  'use strict';

  /* T7 promoter and terminator, as used in the pET vectors.
   *
   * PROVENANCE: standard vector sequences. The promoter core is the widely
   * quoted 18-mer; the terminator is the pET T7 terminator. Both are shown
   * in the interface and can be edited there, because a search motif that
   * cannot be checked should not be buried in the code. */
  var T7_PROMOTER = 'TAATACGACTCACTATAG';
  var T7_TERMINATOR = 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG';
  var FLANK_BASES = 120;

  /* Find a motif on either strand of a circular sequence. */
  function findMotif(record, motif) {
    var forward = record.findAll(motif);
    if (forward.length) {
      return { start: forward[0], end: forward[0] + motif.length - 1,
               strand: 1, hits: forward.length };
    }
    var reverse = record.findAll(PB.reverseComplement(motif));
    if (reverse.length) {
      return { start: reverse[0], end: reverse[0] + motif.length - 1,
               strand: -1, hits: reverse.length };
    }
    return null;
  }

  /* The window to show: 120 bases before the promoter through 120 bases
   * after the terminator, as one continuous stretch.
   *
   * If that span is as long as the plasmid, the whole plasmid is shown
   * instead. Wrapping it further would make the window bite its own tail and
   * silently leave the promoter outside the view, which is worse than showing
   * a little extra.
   */
  function expressionWindow(record, promoterMotif, terminatorMotif, flank) {
    var n = record.seq.length;
    flank = (flank === undefined || flank === null) ? FLANK_BASES : flank;
    var promoter = findMotif(record, promoterMotif || T7_PROMOTER);
    var terminator = findMotif(record, terminatorMotif || T7_TERMINATOR);
    if (!promoter || !terminator) {
      return { found: false, promoter: promoter, terminator: terminator };
    }
    // On the minus strand the promoter sits downstream in plasmid
    // coordinates, so the window runs the other way.
    var forward = promoter.strand > 0;
    var rawFrom = forward ? promoter.start - flank : terminator.start - flank;
    var rawTo = forward ? terminator.end + flank : promoter.end + flank;
    var from, to, wholePlasmid = false;
    if (rawTo - rawFrom + 1 >= n) {
      from = 1; to = n; wholePlasmid = true;
    } else {
      from = PB.wrapPosition(rawFrom, n);
      to = PB.wrapPosition(rawTo, n);
    }
    var between = forward
      ? { start: PB.wrapPosition(promoter.end + 1, n),
          end: PB.wrapPosition(terminator.start - 1, n) }
      : { start: PB.wrapPosition(terminator.end + 1, n),
          end: PB.wrapPosition(promoter.start - 1, n) };
    return { found: true, promoter: promoter, terminator: terminator,
             from: from, to: to, between: between, strand: promoter.strand,
             flank: flank, wholePlasmid: wholePlasmid };
  }

  /* How far past the terminator a stop codon is still looked for.
   *
   * A coding sequence that runs into the terminator region is unusual but
   * real, and until 2026-08-29 such an ORF was invisible here: only the span
   * strictly between the two motifs was scanned, so an ORF whose stop fell
   * outside it was silently skipped and the protein view stayed empty. Start
   * codons are still only accepted between the promoter and the terminator;
   * this only widens where the stop may be. */
  var STOP_SEARCH_NT = 300;

  /* Every ORF between the promoter and the terminator, read on the promoter's
   * strand, longest first.
   *
   * Several ATGs in one frame usually share the same stop codon. Taking the
   * longest means taking the most upstream ATG, and that is a guess: if the
   * real start codon is one of the later ones, every residue number shifts.
   * The mutation-notation mode addresses residues by number, so a wrong guess
   * there does not look wrong, it just designs the wrong primers. All the
   * candidates are returned so the caller can put the choice on screen.
   */
  function orfCandidates(record, between, strand, downstream) {
    if (!between) { return []; }
    var n = record.seq.length;
    downstream = (downstream === undefined || downstream === null)
      ? STOP_SEARCH_NT : downstream;
    var betweenLength = windowPositions(record, between.start, between.end).length;
    // Never let the scan bite its own tail.
    downstream = Math.min(downstream, n - betweenLength);
    if (downstream < 0) { downstream = 0; }

    // The scan runs in the direction of transcription. On the minus strand
    // that is descending plasmid coordinates, so the extension goes below.
    var scanStart = strand > 0
      ? between.start : PB.wrapPosition(between.start - downstream, n);
    var scanEnd = strand > 0
      ? PB.wrapPosition(between.end + downstream, n) : between.end;
    var span = record.subseq(scanStart, scanEnd);
    var region = strand > 0 ? span : PB.reverseComplement(span);
    // In region coordinates the between-part always comes first.
    var startLimit = betweenLength;

    var found = [];
    for (var frame = 0; frame < 3; frame++) {
      for (var i = frame; i + 3 <= region.length && i < startLimit; i += 3) {
        if (region.substr(i, 3) !== 'ATG') { continue; }
        for (var j = i + 3; j + 3 <= region.length; j += 3) {
          var codon = region.substr(j, 3);
          if (codon === 'TAA' || codon === 'TAG' || codon === 'TGA') {
            found.push({ offset: i, length: j + 3 - i,
                         beyondTerminator: j + 3 > startLimit });
            break;
          }
        }
      }
    }
    found.sort(function (a, b) { return b.length - a.length; });
    return found.map(function (entry) {
      return featureFor(record, entry, scanStart, region.length, strand);
    });
  }

  /* One candidate, as a CDS feature in plasmid coordinates. */
  function featureFor(record, entry, scanStart, regionLength, strand) {
    var n = record.seq.length, segments;
    if (strand > 0) {
      var start = PB.wrapPosition(scanStart + entry.offset, n);
      var end = PB.wrapPosition(start + entry.length - 1, n);
      segments = start <= end ? [[start, end]] : [[start, n], [1, end]];
    } else {
      var tail = regionLength - entry.offset - entry.length;
      var lo = PB.wrapPosition(scanStart + tail, n);
      var hi = PB.wrapPosition(lo + entry.length - 1, n);
      segments = lo <= hi ? [[lo, hi]] : [[lo, n], [1, hi]];
    }
    var feature = new PB.Feature('CDS', segments, strand,
                                 { label: 'T7 発現領域の ORF' });
    feature.residues = entry.length / 3 - 1;
    feature.beyondTerminator = !!entry.beyondTerminator;
    return feature;
  }

  /* A CDS the file itself declares, inside the expression window.
   *
   * An annotation the user wrote is better evidence than the longest-ORF
   * guess, and it is what ApE shows them. Preferred whenever it is there.
   */
  function annotatedCds(record, between, strand) {
    if (!between) { return null; }
    var inWindow = {};
    windowPositions(record, between.start, between.end)
      .forEach(function (p) { inWindow[p] = true; });
    var hits = record.cdsFeatures().filter(function (feature) {
      if (strand && feature.strand !== strand) { return false; }
      return feature.segments.every(function (segment) {
        return inWindow[segment[0]] && inWindow[segment[1]];
      });
    });
    if (!hits.length) { return null; }
    hits.sort(function (a, b) {
      return PB.translate(record.featureSequence(b)).protein.length
           - PB.translate(record.featureSequence(a)).protein.length;
    });
    return hits[0];
  }

  /* The ORF to show: the file's own CDS if it declares one, else the longest
   * candidate. Kept for callers that do not offer a choice. */
  function findOrf(record, between, strand) {
    var declared = annotatedCds(record, between, strand);
    if (declared) { return declared; }
    var candidates = orfCandidates(record, between, strand);
    return candidates.length ? candidates[0] : null;
  }

  /* Positions in a window, in display order, wrapping the origin. */
  function windowPositions(record, from, to) {
    var n = record.seq.length;
    var count = ((to - from) % n + n) % n + 1;
    var positions = [];
    for (var i = 0; i < count; i++) {
      positions.push(PB.wrapPosition(from + i, n));
    }
    return positions;
  }

  /* Render a stretch of DNA as one span per base, so a click can name the
   * exact position. Marked regions get a class for colouring. */
  function renderDna(record, positions, marks, perRow) {
    perRow = perRow || 60;
    var classOf = {};
    (marks || []).forEach(function (mark) {
      windowPositions(record, mark.start, mark.end).forEach(function (p) {
        classOf[p] = mark.cls;
      });
    });
    var rows = [];
    for (var offset = 0; offset < positions.length; offset += perRow) {
      var slice = positions.slice(offset, offset + perRow);
      var label = String(slice[0]);
      while (label.length < 6) { label = ' ' + label; }
      var cells = '';
      slice.forEach(function (p, index) {
        if (index && index % 10 === 0) { cells += '<span class="gap"> </span>'; }
        cells += '<span class="b' + (classOf[p] ? ' ' + classOf[p] : '')
          + '" data-pos="' + p + '">' + record.seq[p - 1] + '</span>';
      });
      rows.push('<div class="srow"><span class="pos">' + label
                + '</span> ' + cells + '</div>');
    }
    return rows.join('');
  }

  /* Render a protein as one span per residue, each carrying the plasmid
   * positions of its codon so a selection maps straight back to DNA. */
  function renderProtein(record, feature, perRow) {
    perRow = perRow || 30;
    var map = PB.cdsPositionMap(feature);
    var protein = PB.translate(record.featureSequence(feature)).protein;
    var body = protein.replace(/\*+$/, '');
    var rows = [];
    for (var offset = 0; offset < body.length; offset += perRow) {
      var label = String(offset + 1);
      while (label.length < 6) { label = ' ' + label; }
      var cells = '';
      for (var i = 0; i < perRow && offset + i < body.length; i++) {
        var residue = offset + i + 1;
        if (i && i % 10 === 0) { cells += '<span class="gap"> </span>'; }
        var codon = PB.residueToPlasmid(map, residue);
        cells += '<span class="r" data-res="' + residue
          + '" data-codon="' + codon.join(',') + '">' + body[residue - 1]
          + '</span>';
      }
      rows.push('<div class="srow"><span class="pos">' + label
                + '</span> ' + cells + '</div>');
    }
    return { html: rows.join(''), protein: protein, map: map };
  }

  return {
    T7_PROMOTER: T7_PROMOTER, T7_TERMINATOR: T7_TERMINATOR,
    FLANK_BASES: FLANK_BASES,
    findMotif: findMotif, expressionWindow: expressionWindow,
    findOrf: findOrf, orfCandidates: orfCandidates,
    annotatedCds: annotatedCds, STOP_SEARCH_NT: STOP_SEARCH_NT,
    windowPositions: windowPositions,
    renderDna: renderDna, renderProtein: renderProtein
  };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = GEVIEW; }
