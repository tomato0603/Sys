/* Primer Builder core -- inverse PCR primer design.
 *
 * A direct port of the Python package in this project. Every algorithm and
 * every constant is the same, so the two implementations can be compared
 * against each other; that comparison is part of the test suite.
 *
 * Nothing here touches the network. See the 出典 section of the page.
 */
var PB = (function () {
  'use strict';

  // =====================================================================
  // Constants
  // =====================================================================

  // Takara Bio PrimeSTAR HS recommended conditions.
  // Source: catalog.takara-bio.co.jp (unitid=U100004616), 2026-08-25.
  var DENATURE_TEMP_C = 98, DENATURE_TIME_SEC = 10;
  var ANNEAL_TEMP_C = 55, ANNEAL_TIME_SHORT_SEC = 5, ANNEAL_TIME_LONG_SEC = 15;
  var ANNEAL_TM_THRESHOLD_C = 55;
  var EXTENSION_TEMP_C_3STEP = 72, EXTENSION_TEMP_C_2STEP = 68;
  var EXTENSION_SEC_PER_KB = 60, PCR_CYCLES = 30;
  var PRIMER_AMOUNT_PMOL = '10-15 pmol (each)';
  var TM_FORMULA_MAX_LEN = 25;
  var EXTENSION_ROUNDING_SEC = 30;

  // GenBank / ApE text layout, measured from a file ApE wrote.
  var LOCUS_LENGTH_COL = 40, LOCUS_MOLTYPE_COL = 47;
  var LOCUS_TOPOLOGY_COL = 58, LOCUS_DATE_COL = 71;
  var FEATURE_KEY_COL = 5, FEATURE_LOCATION_COL = 21;
  var GENBANK_LINE_WIDTH = 79;
  var ORIGIN_BASES_PER_LINE = 60, ORIGIN_GROUP_SIZE = 10;
  var GENBANK_MONTHS = ['JAN','FEB','MAR','APR','MAY','JUN',
                        'JUL','AUG','SEP','OCT','NOV','DEC'];

  var UNAMBIGUOUS = 'ACGT';
  var IUPAC = 'ACGTUNRYKMSWBDHV';
  var COMPLEMENT = {A:'T', C:'G', G:'C', T:'A', U:'A', N:'N', R:'Y', Y:'R',
                    K:'M', M:'K', S:'S', W:'W', B:'V', V:'B', D:'H', H:'D'};

  // NCBI genetic code, translation table 1, built from the canonical
  // TCAG-ordered amino acid string so it can be checked by eye.
  var CODE_BASES = 'TCAG';
  var CODE_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  var GENETIC_CODE = {};
  (function () {
    for (var i = 0; i < 4; i++) {
      for (var j = 0; j < 4; j++) {
        for (var k = 0; k < 4; k++) {
          GENETIC_CODE[CODE_BASES[i] + CODE_BASES[j] + CODE_BASES[k]] =
            CODE_AA[i * 16 + j * 4 + k];
        }
      }
    }
  })();
  var STOP_CODON = '*';

  /* The edited plasmid carries a marker over the sequence the edit put in.
   * Without it the file opens in ApE looking like the parent, and the one
   * place that matters is the one place with no label on it. Asked for
   * 2026-09-25. The colour is what ApE shows the feature in. */
  var EDIT_MARKER_LABEL = '改変部位';
  var EDIT_MARKER_COLOR = '#ff8040';

  var MAX_HOMOPOLYMER_RUN = 4, MAX_GC_RUN = 6;
  var DEFAULT_MAX_PRIMER_LEN = 60, DEFAULT_MIN_ANNEAL_LEN = 18;
  var DEFAULT_MAX_ANNEAL_LEN = TM_FORMULA_MAX_LEN;
  var THREE_PRIME_CHECK_LEN = 15;
  var CODON_CONTEXT_NT = 12;
  var RARE_CODON_FRACTION = 0.15;

  function PbError(message) {
    var error = new Error(message);
    error.name = 'PbError';
    return error;
  }

  // =====================================================================
  // Sequence helpers
  // =====================================================================

  function reverseComplement(seq) {
    var out = '';
    var upper = seq.toUpperCase();
    for (var i = upper.length - 1; i >= 0; i--) {
      var base = COMPLEMENT[upper[i]];
      if (!base) { throw PbError('相補塩基が分からない文字です: ' + upper[i]); }
      out += base;
    }
    return out;
  }

  function translate(dna) {
    var s = dna.toUpperCase().replace(/U/g, 'T');
    var trailing = s.length % 3;
    var body = trailing ? s.slice(0, s.length - trailing) : s;
    var protein = '';
    for (var i = 0; i < body.length; i += 3) {
      var aa = GENETIC_CODE[body.substr(i, 3)];
      protein += (aa === undefined ? 'X' : aa);
    }
    return { protein: protein, trailing: trailing };
  }

  /* Occurrences on both strands, overlaps included, palindrome counted once.
   * Every type II restriction site is palindromic, so counting both strands
   * naively doubles every hit and makes "exactly one site" unreachable. */
  function countOccurrences(sequence, query, circular) {
    query = (query || '').toUpperCase();
    if (!query) { return 0; }
    var haystack = sequence.toUpperCase();
    if (circular) { haystack += haystack.slice(0, query.length - 1); }
    var targets = [query];
    var complement = reverseComplement(query);
    if (complement !== query) { targets.push(complement); }
    var total = 0;
    for (var t = 0; t < targets.length; t++) {
      var index = haystack.indexOf(targets[t]);
      while (index !== -1) {
        total++;
        index = haystack.indexOf(targets[t], index + 1);
      }
    }
    return total;
  }

  function longestHomopolymerRun(seq) {
    var best = 0, run = 0, previous = '';
    for (var i = 0; i < seq.length; i++) {
      run = (seq[i] === previous) ? run + 1 : 1;
      previous = seq[i];
      if (run > best) { best = run; }
    }
    return best;
  }

  function longestGcRun(seq) {
    var best = 0, run = 0;
    for (var i = 0; i < seq.length; i++) {
      run = (seq[i] === 'G' || seq[i] === 'C') ? run + 1 : 0;
      if (run > best) { best = run; }
    }
    return best;
  }

  function gcPercent(seq) {
    if (!seq.length) { return 0; }
    var gc = 0;
    for (var i = 0; i < seq.length; i++) {
      if (seq[i] === 'G' || seq[i] === 'C') { gc++; }
    }
    return 100 * gc / seq.length;
  }

  function wrapPosition(position, length) {
    return ((position - 1) % length + length) % length + 1;
  }

  // =====================================================================
  // Record
  // =====================================================================

  /* Coordinates are 1-based and inclusive, matching GenBank and ApE. */
  function Record(fields) {
    this.name = fields.name || 'sequence';
    this.seq = (fields.seq || '').toUpperCase();
    this.isCircular = !!fields.isCircular;
    this.features = fields.features || [];
    this.sourceFormat = fields.sourceFormat || 'unknown';
    this.headerLines = fields.headerLines || [];
    this.declaredLength = fields.declaredLength === undefined
      ? null : fields.declaredLength;
    this.warnings = fields.warnings || [];
  }
  Record.prototype.length = function () { return this.seq.length; };
  Record.prototype.topology = function () {
    return this.isCircular ? 'circular' : 'linear';
  };
  Record.prototype.baseCounts = function () {
    var counts = {};
    for (var i = 0; i < this.seq.length; i++) {
      counts[this.seq[i]] = (counts[this.seq[i]] || 0) + 1;
    }
    return counts;
  };
  Record.prototype.ambiguousCount = function () {
    var counts = this.baseCounts(), total = 0;
    for (var base in counts) {
      if (UNAMBIGUOUS.indexOf(base) === -1) { total += counts[base]; }
    }
    return total;
  };
  /* GC over unambiguous bases only, so a handful of N's does not drag the
   * value toward zero. Returns null when there is nothing to divide by. */
  Record.prototype.gcContent = function () {
    var counts = this.baseCounts();
    var acgt = (counts.A || 0) + (counts.C || 0) + (counts.G || 0) + (counts.T || 0);
    if (!acgt) { return null; }
    return ((counts.G || 0) + (counts.C || 0)) / acgt;
  };
  /* start > end means the span crosses the origin, legal only on a circle. */
  Record.prototype.subseq = function (start, end) {
    var n = this.seq.length;
    if (start < 1 || start > n || end < 1 || end > n) {
      throw PbError('位置 ' + start + '..' + end + ' は 1..' + n + ' の外です');
    }
    if (start <= end) { return this.seq.slice(start - 1, end); }
    if (!this.isCircular) {
      throw PbError('範囲 ' + start + '..' + end
                    + ' は原点をまたぎますが、配列は直鎖です');
    }
    return this.seq.slice(start - 1) + this.seq.slice(0, end);
  };
  Record.prototype.findAll = function (query) {
    var q = (query || '').toUpperCase();
    var n = this.seq.length;
    if (!q || q.length > n) { return []; }
    var haystack = this.seq + (this.isCircular ? this.seq.slice(0, q.length - 1) : '');
    var hits = [], index = haystack.indexOf(q);
    while (index !== -1) {
      hits.push(index + 1);
      index = haystack.indexOf(q, index + 1);
    }
    return hits;
  };
  Record.prototype.cdsFeatures = function () {
    return this.features.filter(function (f) {
      return f.key.toUpperCase() === 'CDS';
    });
  };
  Record.prototype.featureSequence = function (feature) {
    var joined = '';
    for (var i = 0; i < feature.segments.length; i++) {
      joined += this.subseq(feature.segments[i][0], feature.segments[i][1]);
    }
    return feature.strand < 0 ? reverseComplement(joined) : joined;
  };

  function Feature(key, segments, strand, qualifiers) {
    this.key = key;
    this.segments = segments;
    this.strand = strand;
    this.qualifiers = qualifiers || {};
  }
  Feature.prototype.label = function () {
    var order = ['gene', 'locus_tag', 'product', 'label', 'note'];
    for (var i = 0; i < order.length; i++) {
      if (this.qualifiers[order[i]]) { return this.qualifiers[order[i]]; }
    }
    return this.key;
  };
  Feature.prototype.spanLength = function () {
    var total = 0;
    for (var i = 0; i < this.segments.length; i++) {
      total += this.segments[i][1] - this.segments[i][0] + 1;
    }
    return total;
  };
  Feature.prototype.locationText = function () {
    var parts = this.segments.map(function (s) {
      return s[0] === s[1] ? String(s[0]) : s[0] + '..' + s[1];
    });
    var text = parts.length === 1 ? parts[0] : 'join(' + parts.join(',') + ')';
    return this.strand < 0 ? 'complement(' + text + ')' : text;
  };

  // =====================================================================
  // Reading the template: ApE (.ape) or plain DNA text. Nothing else.
  //
  // ApE's native .ape format is GenBank text with extra /ApEinfo_*
  // qualifiers, so the GenBank reader handles it.
  //
  // FASTA is refused, and so is anything that looks like a protein sequence.
  // This tool edits DNA in a plasmid: it needs the plasmid, and a file of
  // amino acids cannot supply one. Refusing loudly beats guessing.
  // =====================================================================

  function detectFormat(text, filename) {
    var isApe = /\.ape$/i.test(filename || '');
    var lines = text.split(/\r?\n/);
    for (var i = 0; i < lines.length; i++) {
      var stripped = lines[i].trim();
      if (!stripped) { continue; }
      if (stripped[0] === '>') {
        throw PbError('FASTA 形式は対象外です。'
          + 'ApE の .ape ファイル、またはプレーンテキストの塩基配列を'
          + '渡してください。');
      }
      if (/^LOCUS/i.test(stripped)) { return 'genbank'; }
      break;
    }
    if (isApe) {
      throw PbError('.ape ですが LOCUS 行で始まっていません。'
        + 'ApE の native 形式は GenBank 形式のはずなので、別の形式の可能性が'
        + 'あります。内容を推測して読むと壊れた結果を出すため、中止します。');
    }
    return 'plain';
  }

  /* Plain text must be DNA. A protein sequence shares several letters with
   * the IUPAC nucleotide codes, so the character check alone would let some
   * through; require the letters to be overwhelmingly A/C/G/T instead.
   *
   * This runs BEFORE the character check, so a protein sequence is reported
   * as a protein sequence rather than as "the letter I is not a base", which
   * says nothing about what to do next. */
  var DNA_FRACTION_REQUIRED = 0.9;

  function requireDna(text, where) {
    var letters = String(text).toUpperCase().replace(/[^A-Z]/g, '');
    if (!letters.length) { return; }
    var acgt = 0;
    for (var i = 0; i < letters.length; i++) {
      if (UNAMBIGUOUS.indexOf(letters[i]) !== -1 || letters[i] === 'U') {
        acgt++;
      }
    }
    if (acgt / letters.length < DNA_FRACTION_REQUIRED) {
      throw PbError(where + ' が塩基配列に見えません（A/C/G/T が '
        + Math.round(100 * acgt / letters.length) + '% しかありません）。'
        + 'タンパク質配列のファイルはこのツールの対象外です。'
        + '改変対象のプラスミドを、ApE の .ape ファイルか'
        + 'プレーンテキストの塩基配列で渡してください。'
        + 'アミノ酸配列を使いたい場合は、「3. どう改変するか」の'
        + '「目的タンパク質の配列」欄に貼ってください。');
    }
  }

  function cleanSequence(raw, where) {
    var seq = raw.replace(/[\s\d]/g, '').toUpperCase();
    for (var i = 0; i < seq.length; i++) {
      if (IUPAC.indexOf(seq[i]) === -1) {
        throw PbError(where + ' に塩基として解釈できない文字があります: '
                      + JSON.stringify(seq[i]));
      }
    }
    return seq;
  }

  /* Supported: ranges, single positions, fuzzy markers, complement(),
   * join(), order() and nesting. Anything else throws: a wrong CDS boundary
   * would shift the reading frame of every downstream calculation. */
  function parseLocation(text) {
    var source = String(text).trim().replace(/\s+/g, '');
    if (!source) { throw PbError('location が空です'); }

    function parse(node, sign) {
      var prefixes = ['complement(', 'join(', 'order('];
      for (var p = 0; p < prefixes.length; p++) {
        var prefix = prefixes[p];
        if (node.toLowerCase().indexOf(prefix) === 0) {
          var inner = node.slice(prefix.length);
          if (inner.slice(-1) !== ')') {
            throw PbError('location の括弧が閉じていません: ' + text);
          }
          inner = inner.slice(0, -1);
          if (prefix === 'complement(') { return parse(inner, -sign); }
          var segments = [], signs = [], parts = splitTopLevel(inner);
          for (var i = 0; i < parts.length; i++) {
            var sub = parse(parts[i], sign);
            segments = segments.concat(sub.segments);
            signs.push(sub.sign);
          }
          var resultSign = signs.some(function (x) { return x < 0; }) ? -1 : sign;
          return { segments: segments, sign: resultSign };
        }
      }
      var match = /^[<>]?(\d+)(?:\.\.[<>]?(\d+))?$/.exec(node);
      if (!match) { throw PbError('解釈できない location です: ' + text); }
      var start = parseInt(match[1], 10);
      var end = match[2] ? parseInt(match[2], 10) : start;
      return { segments: [[start, end]], sign: sign };
    }

    var parsed = parse(source, 1);
    if (!parsed.segments.length) {
      throw PbError('location から座標を取り出せません: ' + text);
    }
    parsed.segments.forEach(function (s) {
      if (s[0] > s[1]) {
        throw PbError('location の始点が終点より後ろです'
          + '（原点をまたぐ CDS は join(...) で表記してください）: ' + text);
      }
    });
    return { segments: parsed.segments, strand: parsed.sign };
  }

  function splitTopLevel(text) {
    var parts = [], depth = 0, current = '';
    for (var i = 0; i < text.length; i++) {
      var ch = text[i];
      if (ch === '(') { depth++; }
      else if (ch === ')') { depth--; }
      if (ch === ',' && depth === 0) { parts.push(current); current = ''; }
      else { current += ch; }
    }
    parts.push(current);
    return parts.filter(function (p) { return p.length; });
  }

  function parseGenbank(text) {
    var lines = text.split(/\r?\n/);
    var name = 'unnamed', declaredLength = null;
    var isCircular = false, topologyStated = false;
    var features = [], headerLines = [], seqChunks = [], warnings = [];
    var section = 'header';
    var pendingKey = null, pendingLocation = [], pendingQuals = {};
    var currentQual = null, currentQualParts = [];

    function flushQualifier() {
      if (currentQual !== null) {
        var value = currentQualParts.join('').trim();
        if (value.length >= 2 && value[0] === '"'
            && value.slice(-1) === '"') {
          value = value.slice(1, -1);
        }
        pendingQuals[currentQual] = value;
      }
      currentQual = null; currentQualParts = [];
    }
    function flushFeature() {
      flushQualifier();
      if (pendingKey !== null) {
        var parsed = parseLocation(pendingLocation.join(''));
        features.push(new Feature(pendingKey, parsed.segments, parsed.strand,
                                  pendingQuals));
      }
      pendingKey = null; pendingLocation = []; pendingQuals = {};
    }

    for (var i = 0; i < lines.length; i++) {
      var line = lines[i];
      if (line.indexOf('//') === 0) { break; }
      var upper = line.toUpperCase();

      if (section === 'header') {
        var locus = /^LOCUS\s+(\S+)\s+(\d+)\s+bp/i.exec(line);
        if (locus) {
          name = locus[1];
          declaredLength = parseInt(locus[2], 10);
          if (/\bcircular\b/i.test(line)) { isCircular = true; topologyStated = true; }
          else if (/\blinear\b/i.test(line)) { isCircular = false; topologyStated = true; }
          continue;
        }
        if (upper.indexOf('FEATURES') === 0) { section = 'features'; continue; }
        if (upper.indexOf('ORIGIN') === 0) { section = 'origin'; continue; }
        // Keep DEFINITION / COMMENT and friends verbatim, trailing spaces
        // included, so writing the record back out is byte-identical.
        if (line.trim()) { headerLines.push(line); }
        continue;
      }

      if (section === 'features') {
        if (upper.indexOf('ORIGIN') === 0) { flushFeature(); section = 'origin'; continue; }
        if (line.length && line[0] !== ' ' && line[0] !== '\t' && line.trim()) {
          flushFeature(); section = 'header'; continue;
        }
        var stripped = line.trim();
        if (!stripped) { continue; }
        // A feature key starts at column 6 (index 5); qualifiers are indented
        // further, to column 22 (index 21).
        var isKey = line.length > 5 && line.slice(0, 5).trim() === ''
          && line[5] !== ' ';
        if (isKey) {
          flushFeature();
          var space = stripped.search(/\s/);
          if (space === -1) { pendingKey = stripped; pendingLocation = []; }
          else {
            pendingKey = stripped.slice(0, space);
            pendingLocation = [stripped.slice(space).trim()];
          }
          continue;
        }
        if (pendingKey === null) { continue; }
        var qual = /^\/([^=]+)(?:=([\s\S]*))?$/.exec(stripped);
        if (qual) {
          flushQualifier();
          currentQual = qual[1];
          currentQualParts = [qual[2] === undefined ? '' : qual[2]];
        } else if (currentQual !== null) {
          currentQualParts.push(' ' + stripped);
        } else {
          pendingLocation.push(stripped);
        }
        continue;
      }

      if (section === 'origin') { seqChunks.push(line); }
    }
    flushFeature();

    var seq = cleanSequence(seqChunks.join(''), 'GenBank の ORIGIN');
    if (!seq) { throw PbError('GenBank の ORIGIN に配列がありません'); }
    if (declaredLength !== null && declaredLength !== seq.length) {
      warnings.push('LOCUS 行の長さ (' + declaredLength + ' bp) と ORIGIN の実長 ('
        + seq.length + ' bp) が一致しません。実長を使用します。');
    }
    if (!topologyStated) {
      warnings.push('LOCUS 行に circular/linear の記載がありません。'
        + 'プラスミドを前提に circular として扱います。');
      isCircular = true;
    }
    return new Record({
      name: name, seq: seq, isCircular: isCircular, features: features,
      sourceFormat: 'genbank', headerLines: headerLines,
      declaredLength: declaredLength, warnings: warnings
    });
  }

  function parsePlain(text, name) {
    requireDna(text, '貼り付けた配列');
    var seq = cleanSequence(text, '配列テキスト');
    return new Record({
      name: name || '配列',
      seq: seq, isCircular: true, sourceFormat: 'plain',
      warnings: ['注釈のない塩基配列として読み込みました。'
                 + '改変対象の CDS の範囲を入力してください。']
    });
  }

  /* .ape files always say "linear" in the LOCUS line, even for plasmids, so
   * that field carries no information for them. Inverse PCR needs a circular
   * template, so circular is the useful default -- announced, not silent. */
  function loadFromText(text, filename, forceTopology) {
    var kind = detectFormat(text, filename);
    var record = kind === 'genbank'
      ? parseGenbank(text) : parsePlain(text, filename || '配列');
    if (kind === 'genbank') { requireDna(record.seq, 'ORIGIN の配列'); }
    if (forceTopology === 'circular') { record.isCircular = true; }
    else if (forceTopology === 'linear') { record.isCircular = false; }
    else if (/\.ape$/i.test(filename || '') && !record.isCircular) {
      record.isCircular = true;
      record.warnings.push('.ape の LOCUS 行は linear ですが、ApE は'
        + 'プラスミドでも常に linear と書き出すため、環状として扱います。');
    }
    return record;
  }

  // =====================================================================
  // Writing GenBank / ApE
  // =====================================================================

  function pad(text, width) {
    while (text.length < width) { text += ' '; }
    return text;
  }

  function genbankDate(date) {
    var d = date || new Date();
    var day = d.getDate();
    return (day < 10 ? '0' + day : String(day)) + '-'
      + GENBANK_MONTHS[d.getMonth()] + '-' + d.getFullYear();
  }

  function formatLocusLine(name, length, isCircular, dateText) {
    var prefix = 'LOCUS       ' + name;
    var number = String(length);
    var gap = LOCUS_LENGTH_COL - prefix.length - number.length;
    var line = prefix + new Array(Math.max(gap, 1) + 1).join(' ') + number + ' bp';
    line = pad(line, LOCUS_MOLTYPE_COL) + 'DNA';
    line = pad(line, LOCUS_TOPOLOGY_COL) + (isCircular ? 'circular' : 'linear');
    line = pad(line, LOCUS_DATE_COL) + dateText;
    return line;
  }

  /* Wrapping only ever happens at a space, and the reader rejoins
   * continuation lines with a single space, so values round-trip. */
  function formatQualifier(key, value) {
    var body;
    if (value === '') { body = '/' + key; }
    else if (/^\d+$/.test(value)) { body = '/' + key + '=' + value; }
    else { body = '/' + key + '="' + value + '"'; }
    var indent = new Array(FEATURE_LOCATION_COL + 1).join(' ');
    var width = GENBANK_LINE_WIDTH - FEATURE_LOCATION_COL;
    var lines = [], current = '';
    body.split(' ').forEach(function (token) {
      if (!current.length) { current = token; }
      else if (current.length + 1 + token.length <= width) {
        current = current + ' ' + token;
      } else { lines.push(indent + current); current = token; }
    });
    lines.push(indent + current);
    return lines;
  }

  function formatOrigin(seq) {
    var lines = [];
    for (var start = 0; start < seq.length; start += ORIGIN_BASES_PER_LINE) {
      var block = seq.substr(start, ORIGIN_BASES_PER_LINE).toLowerCase();
      var groups = [];
      for (var i = 0; i < block.length; i += ORIGIN_GROUP_SIZE) {
        groups.push(block.substr(i, ORIGIN_GROUP_SIZE));
      }
      var number = String(start + 1);
      while (number.length < 9) { number = ' ' + number; }
      lines.push(number + ' ' + groups.join(' '));
    }
    return lines;
  }

  function formatGenbank(record, dateText) {
    var out = [formatLocusLine(record.name, record.seq.length,
                               record.isCircular,
                               dateText || genbankDate())];
    if (record.headerLines && record.headerLines.length) {
      out = out.concat(record.headerLines);
    } else { out.push('DEFINITION  .'); }
    out.push('FEATURES             Location/Qualifiers');
    record.features.forEach(function (feature) {
      out.push(new Array(FEATURE_KEY_COL + 1).join(' ')
        + pad(feature.key, FEATURE_LOCATION_COL - FEATURE_KEY_COL)
        + feature.locationText());
      Object.keys(feature.qualifiers).forEach(function (key) {
        out = out.concat(formatQualifier(key, feature.qualifiers[key]));
      });
    });
    out.push('ORIGIN');
    out = out.concat(formatOrigin(record.seq));
    out.push('//');
    return out.join('\n') + '\n';
  }

  // =====================================================================
  // Tm and PCR conditions  (Takara Bio simple formula)
  // =====================================================================

  /* Tm = 2*(A+T) + 4*(G+C) - 5, stated to apply to 25 mer or shorter.
   * Ambiguous bases throw: the formula defines no weight for them, and a
   * guessed weight would corrupt every downstream decision. */
  function tmTakara(seq) {
    var s = (seq || '').toUpperCase().replace(/U/g, 'T');
    if (!s.length) { throw PbError('空の配列の Tm は計算できません'); }
    var at = 0, gc = 0;
    for (var i = 0; i < s.length; i++) {
      if (s[i] === 'A' || s[i] === 'T') { at++; }
      else if (s[i] === 'G' || s[i] === 'C') { gc++; }
      else { throw PbError('Tm は A/C/G/T のみで計算します。使えない塩基: ' + s[i]); }
    }
    return 2 * at + 4 * gc - 5;
  }

  function tmFormulaApplies(length) { return length <= TM_FORMULA_MAX_LEN; }

  function annealingTimeSec(primerLength, tmAnneal) {
    if (primerLength > TM_FORMULA_MAX_LEN) { return ANNEAL_TIME_SHORT_SEC; }
    return tmAnneal >= ANNEAL_TM_THRESHOLD_C
      ? ANNEAL_TIME_SHORT_SEC : ANNEAL_TIME_LONG_SEC;
  }

  /* The '> 25 mer -> 5 s' rule assumes a long primer is a strong one. An
   * inverse PCR primer is long because of its 5' tail, so warn when the two
   * rules disagree. */
  function annealingTimeNote(primerLength, tmAnneal) {
    if (primerLength > TM_FORMULA_MAX_LEN && tmAnneal < ANNEAL_TM_THRESHOLD_C) {
      return 'プライマー長が ' + primerLength + ' mer なので規則上は '
        + ANNEAL_TIME_SHORT_SEC + ' 秒ですが、3\' 領域の Tm が ' + tmAnneal
        + '℃ (< ' + ANNEAL_TM_THRESHOLD_C + '℃) です。増幅しない場合は '
        + ANNEAL_TIME_LONG_SEC + ' 秒を試してください。';
    }
    return null;
  }

  function extensionTimeSec(productBp) {
    if (productBp <= 0) { throw PbError('産物長は正の整数で指定してください'); }
    return productBp / 1000 * EXTENSION_SEC_PER_KB;
  }

  /* Rounded up to a 30 s grid: a cycler cannot be set to 161.16 s, and
   * too little extension time fails outright while too much rarely hurts. */
  function practicalExtensionTimeSec(productBp) {
    return Math.ceil(extensionTimeSec(productBp) / EXTENSION_ROUNDING_SEC)
      * EXTENSION_ROUNDING_SEC;
  }

  function formatSeconds(seconds) {
    var total = Math.round(seconds);
    var minutes = Math.floor(total / 60), rest = total % 60;
    if (minutes && rest) { return minutes + ' 分 ' + rest + ' 秒'; }
    if (minutes) { return minutes + ' 分'; }
    return rest + ' 秒';
  }

  // =====================================================================
  // Codon usage and selection
  //
  // Measured usage for E. coli: each codon's share within its own amino acid
  // (the shares for one amino acid add up to 1).
  //
  // PROVENANCE: read once out of python_codon_tables (short name
  // "e_coli"); the one-shot script that did it was removed 2026-09-28. Table
  // python_codon_tables ships as e_coli_316407 -- Escherichia coli K-12
  // W3110, taxID 316407 -- which derives from the Kazusa codon usage
  // database. The same numbers are in data/20260907_ecoli-codon-usage.tsv,
  // and tools/parity.py checks that this engine and the Python one still make
  // the same choices. Nothing here was written from memory.
  //
  // Until 2026-09-07 this was a rank ORDER. Rank is ordinal, so it could not
  // say that CGC and CGT are all but interchangeable while AGG is twenty
  // times rarer, and that is exactly what codon choice depends on.
  // =====================================================================

  var ECOLI_CODON_USAGE = {
    // *: TAA 100%, TGA 45%, TAG 11%
    '*': { TAA: 0.6400, TGA: 0.2900, TAG: 0.0700 },
    // A: GCG 100%, GCC 75%, GCA 58%, GCT 44%
    'A': { GCG: 0.3600, GCC: 0.2700, GCA: 0.2100, GCT: 0.1600 },
    // C: TGC 100%, TGT 79%
    'C': { TGC: 0.5600, TGT: 0.4400 },
    // D: GAT 100%, GAC 59%
    'D': { GAT: 0.6300, GAC: 0.3700 },
    // E: GAA 100%, GAG 45%
    'E': { GAA: 0.6900, GAG: 0.3100 },
    // F: TTT 100%, TTC 75%
    'F': { TTT: 0.5700, TTC: 0.4300 },
    // G: GGC 100%, GGT 83%, GGG 37%, GGA 27%
    'G': { GGC: 0.4100, GGT: 0.3400, GGG: 0.1500, GGA: 0.1100 },
    // H: CAT 100%, CAC 75%
    'H': { CAT: 0.5700, CAC: 0.4300 },
    // I: ATT 100%, ATC 82%, ATA 14%
    'I': { ATT: 0.5100, ATC: 0.4200, ATA: 0.0700 },
    // K: AAA 100%, AAG 32%
    'K': { AAA: 0.7600, AAG: 0.2400 },
    // L: CTG 100%, TTA 26%, TTG 26%, CTC 20%, CTT 20%, CTA 8%
    'L': { CTG: 0.5000, TTA: 0.1300, TTG: 0.1300, CTC: 0.1000, CTT: 0.1000, CTA: 0.0400 },
    // M: ATG 100%
    'M': { ATG: 1.0000 },
    // N: AAC 100%, AAT 82%
    'N': { AAC: 0.5500, AAT: 0.4500 },
    // P: CCG 100%, CCA 36%, CCT 30%, CCC 23%
    'P': { CCG: 0.5300, CCA: 0.1900, CCT: 0.1600, CCC: 0.1200 },
    // Q: CAG 100%, CAA 54%
    'Q': { CAG: 0.6500, CAA: 0.3500 },
    // R: CGC 100%, CGT 95%, CGG 25%, CGA 15%, AGA 10%, AGG 5%
    'R': { CGC: 0.4000, CGT: 0.3800, CGG: 0.1000, CGA: 0.0600, AGA: 0.0400, AGG: 0.0200 },
    // S: AGC 100%, AGT 54%, TCC 54%, TCG 54%, TCT 54%, TCA 43%
    'S': { AGC: 0.2800, AGT: 0.1500, TCC: 0.1500, TCG: 0.1500, TCT: 0.1500, TCA: 0.1200 },
    // T: ACC 100%, ACG 61%, ACT 36%, ACA 30%
    'T': { ACC: 0.4400, ACG: 0.2700, ACT: 0.1600, ACA: 0.1300 },
    // V: GTG 100%, GTT 70%, GTC 59%, GTA 41%
    'V': { GTG: 0.3700, GTT: 0.2600, GTC: 0.2200, GTA: 0.1500 },
    // W: TGG 100%
    'W': { TGG: 1.0000 },
    // Y: TAT 100%, TAC 75%
    'Y': { TAT: 0.5700, TAC: 0.4300 },
  };

  // Comparison order for display. Derived, never edited by hand.
  var ECOLI_CODON_PREFERENCE = (function () {
    var out = {};
    Object.keys(ECOLI_CODON_USAGE).forEach(function (aa) {
      var codons = ECOLI_CODON_USAGE[aa];
      out[aa] = Object.keys(codons).sort(function (a, b) {
        return codons[b] - codons[a] || (a < b ? -1 : 1);
      });
    });
    return out;
  })();

  /* How far below the most-used codon a synonym may fall and still be used.
   * At 0.80 only Gly, Ile, Asn and Arg have a second usable codon, so the
   * room this leaves the primer search is real but small. */
  var DEFAULT_MIN_USAGE = 0.80;

  /* Above this many combinations, stop enumerating and fall back to the
   * dynamic program. Never reached for a peptide a primer can carry; it is
   * here so a pathological input cannot hang the page. */
  var MAX_CODON_COMBINATIONS = 2048;

  var ECOLI_RARE_CODONS = ['ATA','AGG','AGA','CTA','CCC','GGA','CGG','CGA'];

  /* 2026-09-28: supplying a measured table was removed. The user settled on
   * E. coli only, so the way in went with it. The numbers themselves are
   * unchanged and still come from data/20260907_ecoli-codon-usage.tsv. */
  function getCodonPreference(host) {
    if ((host || 'ecoli').toLowerCase() !== 'ecoli') {
      throw PbError('組み込みのコドン表は ecoli のみです。');
    }
    return { preference: ECOLI_CODON_PREFERENCE, rare: ECOLI_RARE_CODONS };
  }

  function violatesConstraints(context, addition, avoid, limit) {
    var combined = context + addition;
    var window = limit > 1
      ? combined.slice(Math.max(0, combined.length - (addition.length + limit - 1)))
      : addition;
    if (longestHomopolymerRun(window) > MAX_HOMOPOLYMER_RUN) { return true; }
    if (longestGcRun(window) > MAX_GC_RUN) { return true; }
    for (var i = 0; i < avoid.length; i++) {
      if (window.indexOf(avoid[i]) !== -1) { return true; }
    }
    return false;
  }

  /* Codon shares for a host. */
  function usageShares(host) {
    if ((host || 'ecoli').toLowerCase() !== 'ecoli') {
      throw PbError('組み込みのコドン表は ecoli のみです。');
    }
    return ECOLI_CODON_USAGE;
  }

  /* Codon -> its share of the most-used codon for the same amino acid.
   * 1.0 for the most-used codon. This is what selection maximises. */
  function relativeUsage(shares) {
    var relative = {};
    Object.keys(shares).forEach(function (aa) {
      var codons = shares[aa], best = 0;
      Object.keys(codons).forEach(function (c) {
        if (codons[c] > best) { best = codons[c]; }
      });
      Object.keys(codons).forEach(function (c) {
        relative[c] = best ? codons[c] / best : 0;
      });
    });
    return relative;
  }

  /* How bad a codon is, in thousandths, as an integer.
   *
   * Both engines have to pick the same codons. Summing floats and comparing
   * totals invites a disagreement in the last bit; summing integers cannot.
   * The Python side uses int((1 - w) * 1000 + 0.5), which is this. */
  function usageCost(relativeValue) {
    return Math.floor((1 - relativeValue) * 1000 + 0.5);
  }

  /* Violations created by inserting these codons, counted per junction. */
  function countViolations(upstream, codons, downstream, avoid, limit, keep) {
    var context = upstream.slice(Math.max(0, upstream.length - keep));
    var total = 0;
    codons.forEach(function (codon) {
      if (violatesConstraints(context, codon, avoid, limit)) { total++; }
      var combined = context + codon;
      context = combined.slice(Math.max(0, combined.length - keep));
    });
    if (downstream && violatesConstraints(context, downstream.slice(0, limit),
                                          avoid, limit)) {
      total++;
    }
    return total;
  }

  /* Pick codons for a peptide, favouring the codons the host actually uses.
   *
   * Expression comes first. Only codons used at least `minUsage` as often as
   * the best codon for the same amino acid are considered at all; within that
   * band the choice minimises, in order:
   *
   *   1. constraint violations (long single-base runs, long GC runs, newly
   *      created copies of an avoided site)
   *   2. total distance from the most-used codon
   *
   * Measured usage replaced a rank order on 2026-09-07. There is no separate
   * rare-codon rule any more: a rare codon simply falls outside the band.
   *
   * `candidates` asks for the whole ranked list rather than one answer, so
   * the primer search can trade among strings that are all inside the band.
   * `allowRare` is kept and now means "no band at all".
   */
  function chooseCodons(peptide, options) {
    options = options || {};
    var shares = usageShares(options.host);
    var relative = relativeUsage(shares);
    var minUsage = options.minUsage === undefined || options.minUsage === null
      ? DEFAULT_MIN_USAGE : Number(options.minUsage);
    if (options.allowRare) { minUsage = 0; }
    var maxCost = usageCost(minUsage);

    var avoid = (options.avoid || []).map(function (x) {
      return String(x).toUpperCase();
    }).filter(function (x) { return x.length; });
    var upstream = options.upstream || '', downstream = options.downstream || '';
    var residues = (peptide || '').toUpperCase().split('');

    residues.forEach(function (aa, index) {
      if (!shares[aa]) {
        throw PbError('アミノ酸として解釈できません: ' + JSON.stringify(aa)
                      + '（位置 ' + (index + 1) + '）');
      }
    });
    if (!residues.length) {
      return { peptide: '', dna: '', picks: [], warnings: [],
               level: '新規コドンなし', violations: 0,
               alternatives: [''], usageMean: 0, minUsage: minUsage };
    }

    // Inside the band, cheapest first. The most-used codon always qualifies,
    // so no amino acid can end up with nothing to choose from.
    var allowed = {};
    residues.forEach(function (aa) {
      if (allowed[aa]) { return; }
      allowed[aa] = Object.keys(shares[aa]).filter(function (c) {
        return usageCost(relative[c]) <= maxCost;
      }).sort(function (a, b) {
        return usageCost(relative[a]) - usageCost(relative[b])
               || (a < b ? -1 : 1);
      });
    });

    var longestAvoid = 0;
    avoid.forEach(function (x) { longestAvoid = Math.max(longestAvoid, x.length); });
    var limit = Math.max(MAX_HOMOPOLYMER_RUN + 1, MAX_GC_RUN + 1, longestAvoid);
    var keep = limit + 3;

    var combinations = 1;
    residues.forEach(function (aa) { combinations *= allowed[aa].length; });

    var scored = [];
    if (combinations <= MAX_CODON_COMBINATIONS) {
      // Enumerating is exact and produces a ranked list, which the dynamic
      // program cannot. The band keeps the count small.
      var index = new Array(residues.length);
      for (var i = 0; i < residues.length; i++) { index[i] = 0; }
      for (;;) {
        var combo = residues.map(function (aa, at) { return allowed[aa][index[at]]; });
        var cost = 0;
        combo.forEach(function (c) { cost += usageCost(relative[c]); });
        scored.push({
          violations: countViolations(upstream, combo, downstream, avoid,
                                      limit, keep),
          cost: cost, dna: combo.join('')
        });
        var at = residues.length - 1;
        while (at >= 0 && index[at] + 1 >= allowed[residues[at]].length) {
          index[at] = 0; at--;
        }
        if (at < 0) { break; }
        index[at]++;
      }
      scored.sort(function (a, b) {
        return a.violations - b.violations || a.cost - b.cost
               || (a.dna < b.dna ? -1 : 1);
      });
    } else {
      // Too many to enumerate: fall back to the layered dynamic program over
      // (position, recent context). One answer, no ranked list.
      var states = {};
      states[upstream.slice(Math.max(0, upstream.length - keep))] =
        { violations: 0, cost: 0, path: [] };
      residues.forEach(function (aa) {
        var next = {};
        Object.keys(states).forEach(function (context) {
          var state = states[context];
          allowed[aa].forEach(function (codon) {
            var entry = {
              violations: state.violations
                + (violatesConstraints(context, codon, avoid, limit) ? 1 : 0),
              cost: state.cost + usageCost(relative[codon]),
              path: state.path.concat([codon])
            };
            var combined = context + codon;
            var key = combined.slice(Math.max(0, combined.length - keep));
            var current = next[key];
            if (!current || entry.violations < current.violations
                || (entry.violations === current.violations
                    && entry.cost < current.cost)) {
              next[key] = entry;
            }
          });
        });
        states = next;
      });
      var best = null;
      Object.keys(states).forEach(function (context) {
        var state = states[context];
        var total = state.violations;
        if (downstream && violatesConstraints(context, downstream.slice(0, limit),
                                              avoid, limit)) {
          total++;
        }
        var entry = { violations: total, cost: state.cost,
                      dna: state.path.join('') };
        if (!best || entry.violations < best.violations
            || (entry.violations === best.violations
                && (entry.cost < best.cost
                    || (entry.cost === best.cost && entry.dna < best.dna)))) {
          best = entry;
        }
      });
      scored = [best];
    }

    var chosen = scored[0];
    var codons = [];
    for (var at2 = 0; at2 < chosen.dna.length; at2 += 3) {
      codons.push(chosen.dna.substr(at2, 3));
    }

    var order = {};
    residues.forEach(function (aa) {
      if (order[aa]) { return; }
      order[aa] = Object.keys(shares[aa]).sort(function (a, b) {
        return shares[aa][b] - shares[aa][a] || (a < b ? -1 : 1);
      });
    });

    var picks = [], warnings = [], totalRelative = 0;
    codons.forEach(function (codon, at3) {
      var aa = residues[at3];
      var rank = order[aa].indexOf(codon) + 1;
      var share = relative[codon];
      totalRelative += share;
      var reason = (rank === 1 ? '最頻コドン' : '第 ' + rank + ' 候補')
        + '（最頻の ' + Math.round(100 * share) + '%）';
      if (share < RARE_CODON_FRACTION) {
        reason += ' / レアコドン';
        warnings.push('位置 ' + (at3 + 1) + ' (' + aa + '): ' + codon
          + ' は最頻コドンの ' + Math.round(100 * share) + '% しか使われていません。');
      }
      picks.push({ aa: aa, codon: codon, rank: rank, reason: reason });
    });

    var usageMean = codons.length ? totalRelative / codons.length : 0;
    var level;
    if (chosen.violations) {
      level = '制約違反 ' + chosen.violations + ' 箇所 (前後の配列のため回避できません)';
      warnings.push('同一塩基 ' + (MAX_HOMOPOLYMER_RUN + 1) + ' 以上の連続、GC '
        + (MAX_GC_RUN + 1) + ' 以上の連続、または回避配列を '
        + chosen.violations + ' 箇所で避けられませんでした。'
        + '前後の元配列がその原因です。');
    } else {
      level = '制約違反なし';
    }
    level += '（使用頻度は最頻コドン比で平均 ' + Math.round(100 * usageMean)
      + '%、下限 ' + Math.round(100 * minUsage) + '%）';

    var wanted = Math.max(1, parseInt(options.candidates, 10) || 1);
    return { peptide: residues.join(''), dna: chosen.dna, picks: picks,
             warnings: warnings, level: level, violations: chosen.violations,
             alternatives: scored.slice(0, wanted).map(function (e) { return e.dna; }),
             usageMean: usageMean, minUsage: minUsage };
  }

  // =====================================================================
  // Protein edits
  //
  // Residue numbers are 1-based in the translated CDS, so residue 1 is the
  // initiator methionine. An insertion occupies no original residue and is
  // stored as start = N + 1, end = N, meaning "before residue start".
  // =====================================================================

  function Edit(start, end, original, replacement, source) {
    this.start = start; this.end = end;
    this.original = original || ''; this.replacement = replacement || '';
    this.source = source || '';
  }
  /* Decided from the span, not from `original`: original is only filled in by
   * validateEdits, so a deletion parsed straight from notation would
   * otherwise look like an insertion. */
  Edit.prototype.kind = function () {
    if (this.end < this.start) { return 'insertion'; }
    if (!this.replacement) { return 'deletion'; }
    if (this.start === this.end && this.replacement.length === 1) {
      return 'substitution';
    }
    return 'replacement';
  };
  Edit.prototype.describe = function () {
    var residues = this.original || (this.end - this.start + 1) + ' aa';
    switch (this.kind()) {
      case 'insertion':
        return '残基 ' + (this.start - 1) + ' と ' + this.start + ' の間に '
          + this.replacement + ' を挿入 (' + this.replacement.length + ' aa)';
      case 'deletion':
        return '残基 ' + this.start + '-' + this.end + ' (' + residues
          + ') を欠失 (' + (this.end - this.start + 1) + ' aa)';
      case 'substitution':
        return '残基 ' + this.start + ' を ' + this.original + ' から '
          + this.replacement + ' へ置換';
      default:
        return '残基 ' + this.start + '-' + this.end + ' (' + residues + ') を '
          + this.replacement + ' へ置換 (' + (this.end - this.start + 1)
          + ' aa -> ' + this.replacement.length + ' aa)';
    }
  };

  function parseEdits(text) {
    var edits = [];
    String(text || '').trim().split(/[,\s]+/).forEach(function (token) {
      if (!token) { return; }
      var upper = token.toUpperCase();
      var m = /^([A-Z])(\d+)([A-Z*])$/.exec(upper);
      if (m) {
        var position = parseInt(m[2], 10);
        edits.push(new Edit(position, position, m[1], m[3], token));
        return;
      }
      m = /^ins@(\d+):([A-Za-z*]+)$/i.exec(token);
      if (m) {
        var after = parseInt(m[1], 10);
        edits.push(new Edit(after + 1, after, '', m[2].toUpperCase(), token));
        return;
      }
      m = /^del(\d+)(?:-(\d+))?$/i.exec(token);
      if (m) {
        var start = parseInt(m[1], 10);
        var end = m[2] ? parseInt(m[2], 10) : start;
        if (end < start) { throw PbError('欠失の終点が始点より前です: ' + token); }
        edits.push(new Edit(start, end, '', '', token));
        return;
      }
      m = /^rep(\d+)-(\d+):([A-Za-z*]+)$/i.exec(token);
      if (m) {
        var s = parseInt(m[1], 10), e = parseInt(m[2], 10);
        if (e < s) { throw PbError('置換の終点が始点より前です: ' + token); }
        edits.push(new Edit(s, e, '', m[3].toUpperCase(), token));
        return;
      }
      throw PbError('変異表記として解釈できません: ' + JSON.stringify(token)
        + '\n使える書き方: K42A / ins@50:GGSGGS / del45-50 / rep45-50:GS');
    });
    if (!edits.length) { throw PbError('編集が 1 つも指定されていません'); }
    return sortEdits(edits);
  }

  function sortEdits(edits) {
    return edits.slice().sort(function (a, b) {
      return a.start - b.start || a.end - b.end;
    });
  }

  /* Catches an off-by-one in residue numbering, or a CDS range that is not
   * what the user thought. Raises instead of correcting: either the notation
   * or the CDS is wrong and only the user knows which. */
  function validateEdits(edits, protein) {
    var length = protein.length, cursor = 0;
    sortEdits(edits).forEach(function (edit) {
      if (edit.kind() === 'insertion') {
        if (edit.start < 1 || edit.start > length + 1) {
          throw PbError((edit.source || '挿入') + ': 挿入位置（残基 '
            + (edit.start - 1) + ' の後ろ）がタンパク質の範囲 1-' + length
            + ' の外です');
        }
      } else {
        if (edit.start < 1 || edit.start > length
            || edit.end < 1 || edit.end > length) {
          throw PbError((edit.source || '編集') + ': 残基 ' + edit.start + '-'
            + edit.end + ' はタンパク質の範囲 1-' + length + ' の外です');
        }
        var actual = protein.slice(edit.start - 1, edit.end);
        if (edit.original && edit.original !== actual) {
          throw PbError((edit.source || '編集') + ': 残基 ' + edit.start
            + ' は ' + edit.original + ' ではなく ' + actual
            + ' です。残基番号か CDS の範囲を確認してください。');
        }
        edit.original = actual;
      }
      if (edit.start < cursor) {
        throw PbError('編集の範囲が重なっています: ' + (edit.source || '編集'));
      }
      cursor = edit.end + 1;
    });
  }

  function applyEdits(protein, edits) {
    var out = '', cursor = 1;
    sortEdits(edits).forEach(function (edit) {
      if (edit.start < cursor) { throw PbError('編集の範囲が重なっています'); }
      out += protein.slice(cursor - 1, edit.start - 1) + edit.replacement;
      cursor = edit.end + 1;
    });
    return out + protein.slice(cursor - 1);
  }

  // --- alignment ------------------------------------------------------
  // Identity based, not BLOSUM62: the two sequences are a protein and its own
  // mutant, so what is wanted is the smallest edit script, not evolutionary
  // similarity. It also needs no substitution matrix written from memory.
  var ALIGN_MATCH = 2, ALIGN_MISMATCH = -3;
  var ALIGN_GAP_OPEN = -8, ALIGN_GAP_EXTEND = -1;

  function alignProteins(a, b) {
    a = a.toUpperCase(); b = b.toUpperCase();
    var n = a.length, m = b.length, NEG = -Infinity;
    function grid(fill) {
      var rows = [];
      for (var i = 0; i <= n; i++) {
        rows.push(new Array(m + 1).fill(fill));
      }
      return rows;
    }
    var M = grid(NEG), X = grid(NEG), Y = grid(NEG);
    var pM = grid(null), pX = grid(null), pY = grid(null);
    M[0][0] = 0;
    for (var i = 1; i <= n; i++) {
      X[i][0] = ALIGN_GAP_OPEN + ALIGN_GAP_EXTEND * i;
      pX[i][0] = i > 1 ? 'X' : 'M';
    }
    for (var j = 1; j <= m; j++) {
      Y[0][j] = ALIGN_GAP_OPEN + ALIGN_GAP_EXTEND * j;
      pY[0][j] = j > 1 ? 'Y' : 'M';
    }
    function best(options) {
      var top = options[0];
      for (var k = 1; k < options.length; k++) {
        if (options[k][0] > top[0]) { top = options[k]; }
      }
      return top;
    }
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
        var score = a[i - 1] === b[j - 1] ? ALIGN_MATCH : ALIGN_MISMATCH;
        var pick = best([[M[i-1][j-1], 'M'], [X[i-1][j-1], 'X'], [Y[i-1][j-1], 'Y']]);
        M[i][j] = pick[0] + score; pM[i][j] = pick[1];
        pick = best([[M[i-1][j] + ALIGN_GAP_OPEN + ALIGN_GAP_EXTEND, 'M'],
                     [X[i-1][j] + ALIGN_GAP_EXTEND, 'X'],
                     [Y[i-1][j] + ALIGN_GAP_OPEN + ALIGN_GAP_EXTEND, 'Y']]);
        X[i][j] = pick[0]; pX[i][j] = pick[1];
        pick = best([[M[i][j-1] + ALIGN_GAP_OPEN + ALIGN_GAP_EXTEND, 'M'],
                     [X[i][j-1] + ALIGN_GAP_OPEN + ALIGN_GAP_EXTEND, 'X'],
                     [Y[i][j-1] + ALIGN_GAP_EXTEND, 'Y']]);
        Y[i][j] = pick[0]; pY[i][j] = pick[1];
      }
    }
    var state = best([[M[n][m], 'M'], [X[n][m], 'X'], [Y[n][m], 'Y']])[1];
    var top = [], bottom = [];
    i = n; j = m;
    while (i > 0 || j > 0) {
      if (state === 'M') {
        top.push(a[i - 1]); bottom.push(b[j - 1]);
        state = pM[i][j]; i--; j--;
      } else if (state === 'X') {
        top.push(a[i - 1]); bottom.push('-');
        state = pX[i][j]; i--;
      } else {
        top.push('-'); bottom.push(b[j - 1]);
        state = pY[i][j]; j--;
      }
    }
    return { top: top.reverse().join(''), bottom: bottom.reverse().join('') };
  }

  /* Contiguous differences are merged: in inverse PCR they ride on the same
   * primer tail anyway. */
  function diffProteins(original, target) {
    var aligned = alignProteins(original, target);
    var edits = [], position = 0;
    var runStart = null, runOriginal = [], runReplacement = [];
    function flush() {
      if (runStart === null) { return; }
      edits.push(new Edit(runStart + 1, runStart + runOriginal.length,
                          runOriginal.join(''), runReplacement.join(''), '差分'));
    }
    for (var i = 0; i < aligned.top.length; i++) {
      var t = aligned.top[i], b = aligned.bottom[i];
      if (t === b && t !== '-') {
        flush(); runStart = null; runOriginal = []; runReplacement = [];
        position++;
        continue;
      }
      if (runStart === null) { runStart = position; runOriginal = []; runReplacement = []; }
      if (t !== '-') { runOriginal.push(t); position++; }
      if (b !== '-') { runReplacement.push(b); }
    }
    flush();
    return edits;
  }

  // =====================================================================
  // Grouping edits into rounds
  //
  // One round of inverse PCR has exactly one junction, so the whole new
  // sequence for the edited span rides on the two 5' tails. The unchanged
  // residues *between* two edits are part of that span and are carried too.
  // =====================================================================

  function tailBudgetNt(maxPrimerLen, minAnnealLen) {
    maxPrimerLen = maxPrimerLen || DEFAULT_MAX_PRIMER_LEN;
    minAnnealLen = minAnnealLen || DEFAULT_MIN_ANNEAL_LEN;
    if (minAnnealLen >= maxPrimerLen) {
      throw PbError('アニーリング領域の最小長がプライマー長の上限以上です');
    }
    return 2 * (maxPrimerLen - minAnnealLen);
  }

  function spanNewPeptide(protein, edits, spanStart, spanEnd) {
    var out = '', cursor = spanStart;
    sortEdits(edits).forEach(function (edit) {
      out += protein.slice(cursor - 1, edit.start - 1) + edit.replacement;
      cursor = edit.end + 1;
    });
    return out + protein.slice(cursor - 1, spanEnd);
  }

  /* Residues inside the span that are not being changed keep their ORIGINAL
   * codons. Re-encoding them would change the DNA without changing the
   * protein -- silent synonymous substitutions nobody asked for, which can
   * destroy an existing restriction site or annotation. */
  function buildEditGroup(protein, cdsDna, members, options) {
    options = options || {};
    members = sortEdits(members);
    var spanStart = Math.min.apply(null, members.map(function (e) { return e.start; }));
    var spanEnd = Math.max.apply(null, members.map(function (e) { return e.end; }));
    var originalPeptide = spanEnd >= spanStart
      ? protein.slice(spanStart - 1, spanEnd) : '';
    var newPeptide = spanNewPeptide(protein, members, spanStart, spanEnd);

    var prefix = cdsDna
      ? cdsDna.slice(Math.max(0, (spanStart - 1) * 3 - CODON_CONTEXT_NT),
                     (spanStart - 1) * 3)
      : '';
    var built = [], picks = [], warnings = [], chosenPeptide = [];
    var totalViolations = 0, totalRares = 0, keptNt = 0, chosenNt = 0;
    var cursor = spanStart;

    members.forEach(function (edit) {
      if (edit.start > cursor && cdsDna) {
        var kept = cdsDna.slice((cursor - 1) * 3, (edit.start - 1) * 3);
        built.push(kept); keptNt += kept.length;
      }
      if (edit.replacement) {
        var downstream = cdsDna
          ? cdsDna.slice(edit.end * 3, edit.end * 3 + CODON_CONTEXT_NT) : '';
        var choice = chooseCodons(edit.replacement, {
          host: options.host,
          upstream: prefix + built.join(''), downstream: downstream,
          allowRare: options.allowRare, avoid: options.avoid,
          minUsage: options.minUsage
        });
        built.push(choice.dna); chosenNt += choice.dna.length;
        chosenPeptide.push(edit.replacement);
        picks = picks.concat(choice.picks);
        warnings = warnings.concat(choice.warnings);
        totalViolations += choice.violations;
        totalRares += choice.picks.filter(function (p) {
          return p.reason.indexOf('レアコドン') !== -1;
        }).length;
      }
      cursor = edit.end + 1;
    });
    if (cdsDna && spanEnd >= cursor) {
      var tail = cdsDna.slice((cursor - 1) * 3, spanEnd * 3);
      built.push(tail); keptNt += tail.length;
    }
    var newDna = built.join('');

    /* Compare against the ORIGINAL span: a long GC stretch that was already
     * in the plasmid is not something this design introduced. */
    if (cdsDna && newDna) {
      var suffix = cdsDna.slice(spanEnd * 3, spanEnd * 3 + CODON_CONTEXT_NT);
      var after = prefix + newDna + suffix;
      var before = prefix + cdsDna.slice((spanStart - 1) * 3, spanEnd * 3) + suffix;
      var runAfter = longestHomopolymerRun(after);
      if (runAfter > MAX_HOMOPOLYMER_RUN && runAfter > longestHomopolymerRun(before)) {
        warnings.push('接合部を含めると同一塩基が ' + runAfter + ' 連続します（元の配列では '
          + longestHomopolymerRun(before) + '）。');
      }
      var gcAfter = longestGcRun(after);
      if (gcAfter > MAX_GC_RUN && gcAfter > longestGcRun(before)) {
        warnings.push('接合部を含めると GC が ' + gcAfter + ' 連続します（元の配列では '
          + longestGcRun(before) + '）。');
      }
      (options.avoid || []).forEach(function (site) {
        var upper = String(site).toUpperCase();
        if (upper && after.indexOf(upper) !== -1 && before.indexOf(upper) === -1) {
          warnings.push('回避対象の配列 ' + upper + ' ができています。');
        }
      });
    }

    var level;
    if (!picks.length) { level = '新規コドンなし'; }
    else if (totalViolations) {
      level = '制約違反 ' + totalViolations + ' 箇所（前後の元配列のため回避できません）';
    } else if (totalRares) {
      level = '制約違反なし（レアコドン ' + totalRares + ' 個を使用）';
    } else { level = '制約違反なし'; }

    return {
      edits: members, spanStart: spanStart, spanEnd: spanEnd,
      originalPeptide: originalPeptide, newPeptide: newPeptide, newDna: newDna,
      keptNt: keptNt, chosenNt: chosenNt,
      codonChoice: { picks: picks, warnings: warnings, level: level,
                     violations: totalViolations },
      requiredTailNt: newDna.length,
      removedNt: originalPeptide.length * 3,
      feasible: function (budget) { return this.requiredTailNt <= budget; }
    };
  }

  function groupEdits(protein, edits, options) {
    options = options || {};
    var budget = options.budget || tailBudgetNt();
    var ordered = sortEdits(edits), groups = [], current = [];
    function build(members) {
      return buildEditGroup(protein, options.cdsDna || '', members, options);
    }
    ordered.forEach(function (edit) {
      if (!current.length) { current = [edit]; return; }
      var trial = build(current.concat([edit]));
      if (trial.requiredTailNt <= budget) { current.push(edit); }
      else { groups.push(build(current)); current = [edit]; }
    });
    if (current.length) { groups.push(build(current)); }
    return groups;
  }

  /* Plasmid positions for each base of the CDS, in translation order. For a
   * minus-strand feature the list runs backwards along the plasmid, which is
   * what makes residue numbering work the same way on both strands. */
  function cdsPositionMap(feature) {
    var positions = [];
    feature.segments.forEach(function (segment) {
      for (var p = segment[0]; p <= segment[1]; p++) { positions.push(p); }
    });
    if (feature.strand < 0) { positions.reverse(); }
    return positions;
  }

  function residueToPlasmid(positions, residue) {
    var offset = (residue - 1) * 3;
    if (offset < 0 || offset + 3 > positions.length) {
      throw PbError('残基 ' + residue + ' は CDS の範囲外です（CDS は '
        + Math.floor(positions.length / 3) + ' 残基）');
    }
    return positions.slice(offset, offset + 3);
  }

  // =====================================================================
  // Primer design
  //
  // The primers sit back to back. Reading the ligated circle across the
  // junction gives
  //     ... left_pos | reverseComplement(tailR) + tailF | right_pos ...
  // so the insert is split at one point: tailF carries the right part as it
  // stands, tailR the reverse complement of the left part.
  // =====================================================================

  var WEIGHT_TM_DIFFERENCE = 3.0, WEIGHT_TM_BELOW_TARGET = 4.0;
  var WEIGHT_TM_ABOVE_MAX = 1.0, WEIGHT_GC_OUTSIDE = 0.5;
  var WEIGHT_TERMINAL_NOT_GC = 3.0, WEIGHT_HOMOPOLYMER = 5.0;
  var WEIGHT_LENGTH = 0.2, PENALTY_NON_UNIQUE = 1000.0;
  var TARGET_TM_MIN = ANNEAL_TM_THRESHOLD_C, TARGET_TM_MAX = 68;
  var TARGET_GC_MIN = 40, TARGET_GC_MAX = 60;
  var MAX_HOMOPOLYMER_IN_PRIMER = 4;

  /* A CDS on the minus strand runs backwards along the plasmid, so both the
   * interval and the replacement sequence are flipped for it. */
  function editSite(record, feature, group) {
    var positions = cdsPositionMap(feature);
    var length = record.seq.length;
    var plus = feature.strand > 0;
    var leftPos, rightPos, removed;

    if (group.spanEnd < group.spanStart) {
      var index = (group.spanStart - 1) * 3;
      if (index >= positions.length) {
        throw PbError('挿入位置が CDS の範囲を超えています（残基 '
                      + group.spanStart + '）');
      }
      var anchor = positions[index];
      if (plus) { leftPos = wrapPosition(anchor - 1, length); rightPos = anchor; }
      else { leftPos = anchor; rightPos = wrapPosition(anchor + 1, length); }
      removed = '';
    } else {
      var first = (group.spanStart - 1) * 3, last = group.spanEnd * 3 - 1;
      if (last >= positions.length) {
        throw PbError('残基 ' + group.spanEnd + ' は CDS の範囲を超えています');
      }
      var startPos = positions[first], endPos = positions[last];
      if (plus) {
        leftPos = wrapPosition(startPos - 1, length);
        rightPos = wrapPosition(endPos + 1, length);
        removed = record.subseq(startPos, endPos);
      } else {
        leftPos = wrapPosition(endPos - 1, length);
        rightPos = wrapPosition(startPos + 1, length);
        removed = record.subseq(endPos, startPos);
      }
    }
    return {
      leftPos: leftPos, rightPos: rightPos, removed: removed,
      insert: plus ? group.newDna : reverseComplement(group.newDna),
      removedLength: removed.length,
      insertLength: (plus ? group.newDna : group.newDna).length
    };
  }

  function countBindingSites(record, seq) {
    if (!seq) { return 0; }
    return countOccurrences(record.seq, seq, record.isCircular);
  }

  function measurePrimer(name, tail, anneal, record) {
    var sequence = tail + anneal;
    var tm = tmTakara(anneal);
    var gc = gcPercent(anneal);
    var terminal = anneal.length ? anneal[anneal.length - 1] : '';
    var run = longestHomopolymerRun(sequence);
    var sites = countBindingSites(record, anneal);
    var threePrime = anneal.length >= THREE_PRIME_CHECK_LEN
      ? countBindingSites(record, anneal.slice(-THREE_PRIME_CHECK_LEN))
      : sites;

    var penalties = {};
    if (tm < TARGET_TM_MIN) {
      penalties['Tm が ' + TARGET_TM_MIN + '℃ 未満'] =
        WEIGHT_TM_BELOW_TARGET * (TARGET_TM_MIN - tm);
    }
    if (tm > TARGET_TM_MAX) {
      penalties['Tm が ' + TARGET_TM_MAX + '℃ 超'] =
        WEIGHT_TM_ABOVE_MAX * (tm - TARGET_TM_MAX);
    }
    if (gc < TARGET_GC_MIN) {
      penalties['GC が ' + TARGET_GC_MIN + '% 未満'] =
        WEIGHT_GC_OUTSIDE * (TARGET_GC_MIN - gc);
    }
    if (gc > TARGET_GC_MAX) {
      penalties['GC が ' + TARGET_GC_MAX + '% 超'] =
        WEIGHT_GC_OUTSIDE * (gc - TARGET_GC_MAX);
    }
    if (terminal && 'GC'.indexOf(terminal) === -1) {
      penalties["3' 末端が " + terminal] = WEIGHT_TERMINAL_NOT_GC;
    }
    if (run > MAX_HOMOPOLYMER_IN_PRIMER) {
      penalties['同一塩基が ' + run + ' 連続'] =
        WEIGHT_HOMOPOLYMER * (run - MAX_HOMOPOLYMER_IN_PRIMER);
    }
    penalties['長さ ' + sequence.length + ' nt'] = WEIGHT_LENGTH * sequence.length;
    if (sites !== 1) {
      penalties['結合部位が ' + sites + ' 箇所'] = PENALTY_NON_UNIQUE;
    } else if (threePrime !== 1) {
      penalties["3' 側 " + THREE_PRIME_CHECK_LEN + ' nt が ' + threePrime
                + ' 箇所に一致'] = PENALTY_NON_UNIQUE;
    }
    var total = 0;
    Object.keys(penalties).forEach(function (k) { total += penalties[k]; });

    return { name: name, sequence: sequence, tail: tail, anneal: anneal,
             tm: tm, gcPercent: gc, terminalBase: terminal, longestRun: run,
             bindingSites: sites, threePrimeSites: threePrime,
             penalties: penalties, penaltyTotal: total,
             length: sequence.length,
             unique: sites === 1 && threePrime === 1 };
  }

  /* Where to cut the insert between the two tails: as close to half as the
   * room on each primer allows, or null when it cannot fit. */
  function splitPoint(insertLength, forwardRoom, reverseRoom) {
    var lowest = Math.max(0, insertLength - forwardRoom);
    var highest = Math.min(insertLength, reverseRoom);
    if (lowest > highest) { return null; }
    var midpoint = Math.floor(insertLength / 2);
    return Math.min(Math.max(midpoint, lowest), highest);
  }

  function designPrimerPairs(record, site, options) {
    options = options || {};
    var maxPrimerLen = options.maxPrimerLen || DEFAULT_MAX_PRIMER_LEN;
    var minAnnealLen = options.minAnnealLen || DEFAULT_MIN_ANNEAL_LEN;
    var maxAnnealLen = options.maxAnnealLen || DEFAULT_MAX_ANNEAL_LEN;
    var limit = options.limit || 5;
    if (!record.isCircular) {
      throw PbError('インバースPCR は環状の鋳型が前提です。'
        + '環状として読み込み直してください。');
    }
    var length = record.seq.length, insert = site.insert, results = [];

    for (var fl = minAnnealLen; fl <= maxAnnealLen; fl++) {
      var forwardRoom = maxPrimerLen - fl;
      if (forwardRoom < 0) { continue; }
      var annealF = record.subseq(site.rightPos,
                                  wrapPosition(site.rightPos + fl - 1, length));
      for (var rl = minAnnealLen; rl <= maxAnnealLen; rl++) {
        var reverseRoom = maxPrimerLen - rl;
        if (reverseRoom < 0) { continue; }
        var split = splitPoint(insert.length, forwardRoom, reverseRoom);
        if (split === null) { continue; }
        var annealR = reverseComplement(record.subseq(
          wrapPosition(site.leftPos - rl + 1, length), site.leftPos));

        var forward = measurePrimer('F', insert.slice(split), annealF, record);
        var reverse = measurePrimer('R',
          reverseComplement(insert.slice(0, split)), annealR, record);
        var difference = Math.abs(forward.tm - reverse.tm);
        var pairPenalty = WEIGHT_TM_DIFFERENCE * difference;
        results.push({
          site: site, forward: forward, reverse: reverse, split: split,
          productLength: length - site.removedLength + insert.length,
          tmDifference: difference,
          score: forward.penaltyTotal + reverse.penaltyTotal + pairPenalty,
          usable: forward.unique && reverse.unique
        });
      }
    }
    results.sort(function (a, b) {
      return a.score - b.score
        || (a.forward.length + a.reverse.length) - (b.forward.length + b.reverse.length)
        || (a.forward.sequence < b.forward.sequence ? -1 : 1);
    });
    return results.slice(0, limit);
  }

  function annealingTimeForPair(pair) {
    return Math.max(annealingTimeSec(pair.forward.length, pair.forward.tm),
                    annealingTimeSec(pair.reverse.length, pair.reverse.tm));
  }

  // =====================================================================
  // Verification
  //
  // The check starts from the primer sequences, not from the plan that
  // produced them. The CDS is located in the rebuilt circle by searching for
  // its sequence, so the arithmetic that produced the primers is not reused
  // to grade itself.
  // =====================================================================

  function pcrProduct(record, pair) {
    var retained = record.subseq(pair.site.rightPos, pair.site.leftPos);
    return pair.forward.tail + retained + reverseComplement(pair.reverse.tail);
  }

  function verifyRound(record, feature, protein, group, pair, name) {
    var problems = [], notes = [];
    var product = pcrProduct(record, pair);
    var length = record.seq.length;
    var tailF = pair.forward.tail.length;
    var retainedLength = ((pair.site.leftPos - pair.site.rightPos) % length
                          + length) % length + 1;
    var offsetOfOne = ((1 - pair.site.rightPos) % length + length) % length;
    var circle = product, rotation = 0, keptOrigin = true;
    if (offsetOfOne >= retainedLength) {
      keptOrigin = false;
      notes.push('元の 1 番目の塩基が編集で除去されたため、出力の座標は'
        + 'フォワードプライマーの 5\' 末端を 1 番としています。');
    } else {
      rotation = tailF + offsetOfOne;
      circle = product.slice(rotation) + product.slice(0, rotation);
    }
    var total = circle.length;

    var expectedTotal = length - pair.site.removedLength + pair.site.insertLength;
    if (total !== expectedTotal) {
      problems.push('環状化後の長さが計算と合いません: ' + total + ' bp（期待 '
                    + expectedTotal + ' bp）');
    }

    var newRecord = new Record({
      name: name || record.name, seq: circle, isCircular: true,
      features: [], sourceFormat: record.sourceFormat,
      headerLines: record.headerLines.slice()
    });

    var plus = feature.strand > 0;
    var cdsDna = record.featureSequence(feature);
    var spanDnaStart = (group.spanStart - 1) * 3;
    var spanDnaEnd = group.spanEnd * 3;
    var newCdsDna = cdsDna.slice(0, spanDnaStart) + group.newDna
      + cdsDna.slice(spanDnaEnd);
    var target = plus ? newCdsDna : reverseComplement(newCdsDna);

    var hits = newRecord.findAll(target);
    var newCds = null, observed = '';
    if (hits.length !== 1) {
      problems.push('改変後の CDS 配列が環状化後のプラスミドに ' + hits.length
        + ' 箇所見つかりました（1 箇所であるべきです）');
    } else {
      var start = hits[0], end = start + target.length - 1;
      var segments = end <= total ? [[start, end]]
        : [[start, total], [1, end - total]];
      newCds = new Feature('CDS', segments, feature.strand,
                           Object.assign({}, feature.qualifiers));
      var translated = translate(newRecord.featureSequence(newCds));
      observed = translated.protein;
      if (translated.trailing) {
        problems.push('改変後の CDS の長さが 3 の倍数ではありません（余り '
                      + translated.trailing + ' 塩基）');
      }
    }

    var expected = applyEdits(protein, group.edits);
    if (newCds && observed !== expected) {
      problems.push('翻訳結果が目的のタンパク質と一致しません（期待 '
        + expected.length + ' 残基 / 実際 ' + observed.length + ' 残基）');
    }
    if (observed) {
      var internal = observed.slice(0, -1).split(STOP_CODON).length - 1;
      if (internal) {
        problems.push('改変後の CDS に内部終止コドンが ' + internal + ' 個あります');
      }
      if (observed.slice(-1) !== STOP_CODON) {
        problems.push('改変後の CDS が終止コドンで終わっていません');
      }
    }

    // Carry over the features that survive untouched; drop the rest loudly.
    var mapping = {};
    for (var step = 0; step < retainedLength; step++) {
      var old = wrapPosition(pair.site.rightPos + step, length);
      mapping[old] = ((tailF + step - rotation) % total + total) % total + 1;
    }
    var dropped = [];
    record.features.forEach(function (existing) {
      if (existing === feature) { return; }
      if (existing.key.toLowerCase() === 'source'
          && existing.segments.length === 1
          && existing.segments[0][0] === 1
          && existing.segments[0][1] === length) {
        // `source` spans the whole record by definition; rebuild it.
        newRecord.features.push(new Feature('source', [[1, total]],
          existing.strand, Object.assign({}, existing.qualifiers)));
        return;
      }
      var moved = [], intact = true;
      for (var i = 0; i < existing.segments.length && intact; i++) {
        var segment = existing.segments[i];
        for (var p = segment[0]; p <= segment[1]; p++) {
          if (!(wrapPosition(p, length) in mapping)) { intact = false; break; }
        }
        if (intact) {
          moved.push([mapping[wrapPosition(segment[0], length)],
                      mapping[wrapPosition(segment[1], length)]]);
        }
      }
      var ordered = intact && moved.every(function (s) { return s[0] <= s[1]; });
      if (ordered) {
        newRecord.features.push(new Feature(existing.key, moved,
          existing.strand, Object.assign({}, existing.qualifiers)));
      } else { dropped.push(existing); }
    });
    /* Mark the bases the edit put in. They sit inside the new CDS, so the
     * CDS is walked in reading order to find where they landed -- backwards
     * along the plasmid when the CDS is on the minus strand. */
    var marker = null;
    if (newCds && group.newDna.length) {
      var reading = [];
      newCds.segments.forEach(function (seg) {
        for (var q = seg[0]; q <= seg[1]; q++) { reading.push(q); }
      });
      if (newCds.strand < 0) { reading.reverse(); }
      var span = reading.slice(spanDnaStart, spanDnaStart + group.newDna.length);
      if (span.length === group.newDna.length) {
        marker = editMarkerFeature(positionsToSegments(span, total),
                                   spanDnaEnd - spanDnaStart,
                                   group.newDna.length);
        if (marker) { newRecord.features.push(marker); }
      }
    }
    if (newCds) { newRecord.features.unshift(newCds); }
    if (dropped.length) {
      notes.push('編集部位にかかる、または原点をまたぐようになったフィーチャー '
        + dropped.length + ' 件を出力から外しました: '
        + dropped.map(function (f) { return f.label(); }).join(', '));
    }

    return { newRecord: newRecord, newCds: newCds, expectedProtein: expected,
             observedProtein: observed, productLength: product.length,
             editMarker: marker,
             droppedFeatures: dropped, problems: problems, notes: notes,
             ok: problems.length === 0 };
  }

  // =====================================================================
  // Diagnostic restriction sites
  //
  // Screening by digest needs a site the parent plasmid lacks and the edited
  // plasmid has exactly once. Only bases inside the replaced span may change,
  // and only synonymously.
  //
  // PROVENANCE: recognition sequences written from standard supplier
  // catalogues. Widely used and stable, but not checked against a database
  // from inside this program -- verify before ordering.
  // =====================================================================

  var BUILTIN_ENZYMES = {
    AflII:'CTTAAG', AgeI:'ACCGGT', ApaI:'GGGCCC', AscI:'GGCGCGCC',
    AvrII:'CCTAGG', BamHI:'GGATCC', BglII:'AGATCT', BsrGI:'TGTACA',
    ClaI:'ATCGAT', DraI:'TTTAAA', EcoRI:'GAATTC', EcoRV:'GATATC',
    FseI:'GGCCGGCC', HindIII:'AAGCTT', HpaI:'GTTAAC', KpnI:'GGTACC',
    MluI:'ACGCGT', NcoI:'CCATGG', NdeI:'CATATG', NheI:'GCTAGC',
    NotI:'GCGGCCGC', PacI:'TTAATTAA', PmeI:'GTTTAAAC', PstI:'CTGCAG',
    PvuII:'CAGCTG', SacI:'GAGCTC', SalI:'GTCGAC', SbfI:'CCTGCAGG',
    ScaI:'AGTACT', SmaI:'CCCGGG', SpeI:'ACTAGT', SphI:'GCATGC',
    SspI:'AATATT', StuI:'AGGCCT', SwaI:'ATTTAAAT', XbaI:'TCTAGA',
    XhoI:'CTCGAG'
  };

  function parseEnzymeTable(text) {
    var table = {};
    text.split(/\r?\n/).forEach(function (line, index) {
      var body = line.split('#')[0].trim();
      if (!body) { return; }
      var parts = body.replace(/,/g, ' ').split(/\s+/);
      if (parts.length < 2) {
        throw PbError((index + 1) + ' 行目: 「名前 認識配列」の形式ではありません: ' + line);
      }
      var site = parts[1].toUpperCase().replace(/U/g, 'T');
      for (var i = 0; i < site.length; i++) {
        if (UNAMBIGUOUS.indexOf(site[i]) === -1) {
          throw PbError((index + 1) + ' 行目: 認識配列に A/C/G/T 以外があります: '
                        + site[i]);
        }
      }
      table[parts[0]] = site;
    });
    if (!Object.keys(table).length) { throw PbError('酵素が 1 つもありません'); }
    return table;
  }

  function synonymousCodons(codon) {
    var aa = GENETIC_CODE[codon];
    return Object.keys(GENETIC_CODE).filter(function (c) {
      return GENETIC_CODE[c] === aa;
    });
  }

  /* Bases outside [mutableStart, mutableEnd] must already match, since the
   * primer annealing regions have to stay identical to the template. */
  function placeSite(newCds, recognition, offset, mutableStart, mutableEnd,
                     preference) {
    var end = offset + recognition.length - 1;
    if (end >= newCds.length) { return null; }
    if (end < mutableStart || offset > mutableEnd) { return null; }

    var required = {};
    for (var i = 0; i < recognition.length; i++) {
      var position = offset + i;
      if (position >= mutableStart && position <= mutableEnd) {
        required[position] = recognition[i];
      } else if (newCds[position] !== recognition[i]) { return null; }
    }

    var rewritten = newCds.split('');
    var changed = [];
    var codonIndexes = {};
    Object.keys(required).forEach(function (p) {
      codonIndexes[Math.floor(parseInt(p, 10) / 3)] = true;
    });
    var ordered = Object.keys(codonIndexes).map(Number).sort(function (a, b) {
      return a - b;
    });
    for (var k = 0; k < ordered.length; k++) {
      var codonIndex = ordered[k];
      var original = newCds.substr(codonIndex * 3, 3);
      if (original.length < 3 || !GENETIC_CODE[original]) { return null; }
      var constraints = {};
      Object.keys(required).forEach(function (p) {
        var position = parseInt(p, 10);
        if (Math.floor(position / 3) === codonIndex) {
          constraints[position % 3] = required[position];
        }
      });
      for (var slot = 0; slot < 3; slot++) {
        var absolute = codonIndex * 3 + slot;
        if (absolute < mutableStart || absolute > mutableEnd) {
          if (!(slot in constraints)) { constraints[slot] = newCds[absolute]; }
        }
      }
      var options = synonymousCodons(original).filter(function (candidate) {
        return Object.keys(constraints).every(function (slot) {
          return candidate[slot] === constraints[slot];
        });
      });
      if (!options.length) { return null; }
      var order = preference[GENETIC_CODE[original]] || options;
      options.sort(function (a, b) {
        var ia = order.indexOf(a), ib = order.indexOf(b);
        return (ia === -1 ? order.length : ia) - (ib === -1 ? order.length : ib);
      });
      var chosen = options[0];
      if (chosen !== original) { changed.push([codonIndex, original, chosen]); }
      for (var s = 0; s < 3; s++) { rewritten[codonIndex * 3 + s] = chosen[s]; }
    }
    var result = rewritten.join('');
    if (result.slice(offset, end + 1) !== recognition) { return null; }
    return { cds: result, changed: changed };
  }

  /* Only enzymes absent from the parent are considered, so a hit means the
   * edited plasmid is cut exactly once and the parent is not cut at all. */
  function findDiagnosticSites(record, feature, group, site, options) {
    options = options || {};
    var enzymes = options.enzymes || BUILTIN_ENZYMES;
    var preference = getCodonPreference(options.host).preference;
    var plus = feature.strand > 0;
    var cdsDna = record.featureSequence(feature);
    var spanStart = (group.spanStart - 1) * 3;
    var spanEnd = group.spanEnd * 3;
    var newCds = cdsDna.slice(0, spanStart) + group.newDna + cdsDna.slice(spanEnd);
    var mutableStart = spanStart;
    var mutableEnd = spanStart + group.newDna.length - 1;
    if (mutableEnd < mutableStart) { return []; }

    var results = [];
    Object.keys(enzymes).sort().forEach(function (name) {
      var recognition = enzymes[name];
      if (countOccurrences(record.seq, recognition, record.isCircular)) { return; }
      var from = Math.max(0, mutableStart - recognition.length + 1);
      for (var offset = from; offset <= mutableEnd; offset++) {
        var placed = placeSite(newCds, recognition, offset, mutableStart,
                               mutableEnd, preference);
        if (!placed || !placed.changed.length) { continue; }
        var variant = placed.cds.slice(mutableStart, mutableEnd + 1);
        if (variant.length !== group.newDna.length) { continue; }
        var insert = plus ? variant : reverseComplement(variant);
        var candidate = record.subseq(site.rightPos, site.leftPos) + insert;
        if (countOccurrences(candidate, recognition, true) !== 1) { continue; }
        var changedBases = 0;
        placed.changed.forEach(function (entry) {
          for (var i = 0; i < 3; i++) {
            if (entry[1][i] !== entry[2][i]) { changedBases++; }
          }
        });
        results.push({ enzyme: name, recognition: recognition, newDna: variant,
                       offsetInCds: offset, changedCodons: placed.changed,
                       changedBases: changedBases, siteCountAfter: 1 });
        break;
      }
    });
    results.sort(function (a, b) {
      return a.changedBases - b.changedBases
        || a.changedCodons.length - b.changedCodons.length
        || (a.enzyme < b.enzyme ? -1 : 1);
    });
    return results.slice(0, options.limit || 5);
  }

  // =====================================================================
  // Replacement mode: "this region becomes this sequence"
  //
  // The edit is stated at the DNA level -- a selected span of the plasmid and
  // what should be there instead -- so no CDS and no residue numbering are
  // needed. The primer search and the junction equation are the ones already
  // used by the mutation-notation path; only the entry point is different.
  //
  // The plasmid itself is still required: the 3' annealing regions come from
  // the template OUTSIDE the selection, so "old fragment + new fragment"
  // alone cannot produce primers.
  // =====================================================================

  /* Turn a selected span into an EditSite.
   *
   * end === start - 1 (or the wrapped equivalent) means an insertion before
   * `start` with nothing removed. A span may cross the origin, in which case
   * start > end.
   */
  /* Where a selected span sits, with no opinion about what replaces it.
   *
   * Split out of replacementSite because the flanking context is needed
   * BEFORE there is an insert: encoding a peptide looks at the bases on
   * either side, and asking replacementSite for that meant handing it an
   * empty insert, which it rightly refuses for an insertion. */
  function spanBounds(record, start, end) {
    var n = record.seq.length;
    if (!record.isCircular) {
      throw PbError('インバースPCR は環状の鋳型が前提です。');
    }
    if (!(start >= 1 && start <= n)) {
      throw PbError('選択の開始位置 ' + start + ' が 1..' + n + ' の外です');
    }
    var empty = (end === start - 1)
      || (wrapPosition(end + 1, n) === start && end === wrapPosition(start - 1, n));
    if (empty) {
      return { leftPos: wrapPosition(start - 1, n), rightPos: start,
               removed: '', removedLength: 0, empty: true };
    }
    if (!(end >= 1 && end <= n)) {
      throw PbError('選択の終了位置 ' + end + ' が 1..' + n + ' の外です');
    }
    var removed = record.subseq(start, end);
    if (removed.length > n - 2 * DEFAULT_MIN_ANNEAL_LEN) {
      throw PbError('選択範囲が長すぎます（' + removed.length + ' bp）。'
        + 'プライマーが結合する領域が残りません。');
    }
    return { leftPos: wrapPosition(start - 1, n),
             rightPos: wrapPosition(end + 1, n),
             removed: removed, removedLength: removed.length, empty: false };
  }

  function replacementSite(record, start, end, newDna) {
    var bounds = spanBounds(record, start, end);
    var insert = (newDna || '').toUpperCase().replace(/\s/g, '');
    for (var i = 0; i < insert.length; i++) {
      if (UNAMBIGUOUS.indexOf(insert[i]) === -1) {
        throw PbError('新しい配列に A/C/G/T 以外があります: ' + insert[i]);
      }
    }
    if (!insert.length && !bounds.removedLength) {
      throw PbError('何も除去せず、何も挿入しない指定です。'
        + '置き換える領域を選ぶか、新しい配列を入力してください。');
    }
    return { leftPos: bounds.leftPos, rightPos: bounds.rightPos,
             removed: bounds.removed, insert: insert,
             removedLength: bounds.removedLength, insertLength: insert.length };
  }

  /* Find a pasted selection in the plasmid and return its coordinates.
   *
   * Refuses anything but a single hit: a sequence that occurs twice gives no
   * way to know which copy was selected, and guessing would edit the wrong
   * one. A hit on the reverse strand is reported with its forward coordinates
   * so the rest of the pipeline stays on one strand.
   */
  function locateSelection(record, selected) {
    var query = (selected || '').toUpperCase().replace(/[\s\d]/g, '');
    if (!query) { throw PbError('選択した配列が空です'); }
    for (var i = 0; i < query.length; i++) {
      if (IUPAC.indexOf(query[i]) === -1) {
        throw PbError('選択した配列に塩基として解釈できない文字があります: '
                      + query[i]);
      }
    }
    if (query.length > record.seq.length) {
      throw PbError('選択した配列が鋳型より長いです');
    }
    var n = record.seq.length;
    var forward = record.findAll(query);
    var complement = reverseComplement(query);
    var reverse = complement === query ? [] : record.findAll(complement);
    var total = forward.length + reverse.length;
    if (total === 0) {
      throw PbError('選択した配列が鋳型の中に見つかりません（' + query.length
        + ' bp）。同じプラスミドの配列か、コピー漏れがないか確認してください。');
    }
    if (total > 1) {
      var places = forward.concat(reverse).slice(0, 6).join(', ');
      throw PbError('選択した配列が鋳型の中に ' + total + ' 箇所あります（位置 '
        + places + (total > 6 ? ' ほか' : '') + '）。どこを指しているか'
        + '決められないので、座標で指定するか、前後を含めて長めに選び直して'
        + 'ください。');
    }
    var start = forward.length ? forward[0] : reverse[0];
    var end = wrapPosition(start + query.length - 1, n);
    return { start: start, end: end, length: query.length,
             onReverseStrand: forward.length === 0,
             crossesOrigin: end < start };
  }

  /* Encode a peptide into DNA for a site, using the flanking template as
   * context so a run created at the junction is counted. */
  /* The cds argument as a Feature. Accepts a Feature or a location string. */
  function asCdsFeature(cds) {
    if (!cds) { return null; }
    if (typeof cds === 'string') {
      var parsed = parseLocation(cds);
      return new Feature('CDS', parsed.segments, parsed.strand, {});
    }
    return cds;
  }

  function encodePeptideForSite(record, site, peptide, options) {
    options = options || {};
    var n = record.seq.length;
    var before = record.subseq(
      wrapPosition(site.leftPos - CODON_CONTEXT_NT + 1, n), site.leftPos);
    var after = record.subseq(
      site.rightPos, wrapPosition(site.rightPos + CODON_CONTEXT_NT - 1, n));
    /* On the minus strand the CDS reads from high coordinates to low, so what
     * comes before the peptide in reading order is the plus-strand sequence
     * AFTER the site, reverse complemented -- and the codons chosen have to be
     * turned back to the plus strand before they can be inserted.
     *
     * Until 2026-09-11 this was ignored: the peptide was encoded as if the CDS
     * always read forwards, so asking for WW on a minus-strand CDS inserted
     * TGGTGG, which that CDS reads as PP. Nothing caught it because
     * verifyReplacement compares DNA, never translating. */
    var strand = options.strand === undefined ? 1 : options.strand;
    var upstream = strand < 0 ? reverseComplement(after) : before;
    var downstream = strand < 0 ? reverseComplement(before) : after;
    var choice = chooseCodons(peptide, {
      host: options.host,
      upstream: upstream, downstream: downstream,
      allowRare: options.allowRare, avoid: options.avoid,
      minUsage: options.minUsage,
      candidates: options.codonCandidates === undefined
        ? 8 : options.codonCandidates
    });
    if (strand < 0) {
      choice.dna = reverseComplement(choice.dna);
      choice.alternatives = choice.alternatives.map(function (x) {
        return reverseComplement(x);
      });
    }
    return choice;
  }

  /* Rebuild the circle FROM THE PRIMERS and compare it with the circle the
   * request describes, computed independently from the template.
   *
   * This is the whole check: if the two agree base for base, the primers
   * produce exactly the intended plasmid, and nothing outside the selection
   * has moved. It needs no CDS and no translation.
   */
  /* A misc_feature over the newly inserted bases.
   *
   * `segments` is already in GenBank order, so an origin-crossing span comes
   * in as [[start, total], [1, end]]. A pure deletion gets no marker: there
   * is no span left to mark, and a zero-length feature cannot be written. */
  function editMarkerFeature(segments, removedLength, insertLength) {
    if (!insertLength || !segments.length) { return null; }
    var note = removedLength
      ? '元の ' + removedLength + ' bp を ' + insertLength + ' bp に置き換え'
      : insertLength + ' bp を挿入';
    return new Feature('misc_feature', segments, 1, {
      label: EDIT_MARKER_LABEL,
      note: 'Primer Builder（インバースPCR）: ' + note,
      ApEinfo_fwdcolor: EDIT_MARKER_COLOR,
      ApEinfo_revcolor: EDIT_MARKER_COLOR
    });
  }

  /* Group plasmid positions into GenBank segments, origin crossing included.
   * The positions come in reading order, which for a minus-strand CDS runs
   * backwards, so they are sorted first. */
  function positionsToSegments(positions, total) {
    var sorted = positions.slice().sort(function (a, b) { return a - b; });
    var runs = [], start = sorted[0], prev = sorted[0];
    for (var i = 1; i < sorted.length; i++) {
      if (sorted[i] !== prev + 1) { runs.push([start, prev]); start = sorted[i]; }
      prev = sorted[i];
    }
    runs.push([start, prev]);
    if (runs.length === 2 && runs[0][0] === 1 && runs[1][1] === total) {
      runs = [runs[1], runs[0]];
    }
    return runs;
  }

  function verifyReplacement(record, site, pair) {
    var problems = [], notes = [];
    var retained = record.subseq(site.rightPos, site.leftPos);
    var expected = retained + site.insert;

    var product = pcrProduct(record, pair);
    var tailF = pair.forward.tail.length;
    var actual = product.slice(tailF) + product.slice(0, tailF);

    if (actual.length !== expected.length) {
      problems.push('環状化後の長さが合いません: ' + actual.length
                    + ' bp（期待 ' + expected.length + ' bp）');
    } else if (actual !== expected) {
      var at = 0;
      while (at < actual.length && actual[at] === expected[at]) { at++; }
      problems.push('環状化後の配列が期待と違います（' + (at + 1)
                    + ' 番目の塩基から食い違います）');
    }

    var expectedLength = record.seq.length - site.removedLength
      + site.insertLength;
    if (expected.length !== expectedLength) {
      problems.push('長さの計算が合いません: ' + expected.length + ' bp（期待 '
                    + expectedLength + ' bp）');
    }
    if (retained !== record.subseq(site.rightPos, site.leftPos)) {
      problems.push('選択範囲の外の配列が変化しています');
    }

    // Present the result with the original origin first, so the coordinates
    // still line up with the plasmid map the user already has.
    var n = record.seq.length;
    var offsetOfOne = ((1 - site.rightPos) % n + n) % n;
    var circle = actual, keptOrigin = true;
    if (offsetOfOne >= retained.length) {
      keptOrigin = false;
      notes.push('元の 1 番目の塩基が選択範囲に含まれていたため、'
        + '出力の座標は選択範囲の直後を 1 番としています。');
    } else {
      circle = actual.slice(offsetOfOne) + actual.slice(0, offsetOfOne);
    }

    /* Carry over the features that survive untouched; drop the rest loudly.
     *
     * The mutation-notation path does exactly this. Until 2026-08-29 this one
     * did not, so the same plasmid came out annotated from one mode and bare
     * from the other -- a difference nobody sees until the file is open in
     * ApE and every label is gone. */
    var total = circle.length;
    var shift = keptOrigin ? offsetOfOne : 0;
    var mapping = {};
    for (var step = 0; step < retained.length; step++) {
      mapping[wrapPosition(site.rightPos + step, n)] =
        ((step - shift) % total + total) % total + 1;
    }
    var carried = [], dropped = [];
    record.features.forEach(function (existing) {
      if (existing.key.toLowerCase() === 'source'
          && existing.segments.length === 1
          && existing.segments[0][0] === 1
          && existing.segments[0][1] === n) {
        // `source` spans the whole record by definition; rebuild it.
        carried.push(new Feature('source', [[1, total]], existing.strand,
          Object.assign({}, existing.qualifiers)));
        return;
      }
      var moved = [], intact = true;
      for (var i = 0; i < existing.segments.length && intact; i++) {
        var segment = existing.segments[i];
        for (var q = segment[0]; q <= segment[1]; q++) {
          if (!(wrapPosition(q, n) in mapping)) { intact = false; break; }
        }
        if (intact) {
          moved.push([mapping[wrapPosition(segment[0], n)],
                      mapping[wrapPosition(segment[1], n)]]);
        }
      }
      if (intact && moved.every(function (seg) { return seg[0] <= seg[1]; })) {
        carried.push(new Feature(existing.key, moved, existing.strand,
          Object.assign({}, existing.qualifiers)));
      } else { dropped.push(existing); }
    });
    if (dropped.length) {
      notes.push('編集部位にかかる、または原点をまたぐようになったフィーチャー '
        + dropped.length + ' 件を出力から外しました: '
        + dropped.map(function (f) { return f.label(); }).join(', '));
    }

    var newRecord = new Record({
      name: record.name, seq: circle, isCircular: true, features: carried,
      sourceFormat: record.sourceFormat,
      headerLines: record.headerLines.slice()
    });

    // Where the new sequence ended up, in the new record's coordinates.
    var insertStart = null, insertEnd = null;
    if (site.insertLength) {
      var offset = keptOrigin
        ? ((retained.length - offsetOfOne) % circle.length) : retained.length;
      insertStart = offset + 1;
      insertEnd = offset + site.insertLength;
      if (insertEnd > circle.length) { insertEnd -= circle.length; }
      if (insertStart <= circle.length && insertEnd >= insertStart) {
        if (newRecord.subseq(insertStart, insertEnd) !== site.insert) {
          problems.push('新しい配列が想定した位置にありません');
        }
      }
    }

    var marker = null;
    if (insertStart !== null) {
      marker = editMarkerFeature(
        insertStart <= insertEnd ? [[insertStart, insertEnd]]
          : [[insertStart, circle.length], [1, insertEnd]],
        site.removedLength, site.insertLength);
      if (marker) { newRecord.features.push(marker); }
    }

    return { newRecord: newRecord, expectedSequence: expected,
             insertStart: insertStart, insertEnd: insertEnd,
             editMarker: marker,
             keptOrigin: keptOrigin, productLength: product.length,
             problems: problems, notes: notes, ok: problems.length === 0 };
  }

  /* Is an amino-acid replacement in frame?
   *
   * A peptide is only the peptide you asked for if its codons start on a
   * codon boundary. Replacing a span whose length is a multiple of three is
   * not enough: the span has to BEGIN on a boundary too, and a pure insertion
   * removes nothing at all, so the length test says nothing about it. That
   * gap is reachable whenever coordinates are typed by hand.
   *
   * The frame can only be checked against a CDS. Without one, say so rather
   * than let silence pass for approval.
   */
  function checkReadingFrame(record, cds, start, removedLength) {
    if (!cds) {
      return { note: '読み枠の基準となる CDS が指定されていないため、'
        + '挿入位置がコドンの境界にあるかを確認できませんでした。'
        + 'ずれていると、指定したアミノ酸配列とは違うものが翻訳されます。'
        + 'CDS の範囲を指定すると確認します。' };
    }
    var feature = asCdsFeature(cds);
    var positions = cdsPositionMap(feature);
    var offset = positions.indexOf(start);
    if (offset === -1) {
      return { error: '選択の開始位置 ' + start + ' は指定された CDS ('
        + feature.locationText() + ') の中にありません。'
        + 'アミノ酸で指定するには、改変する場所が CDS の中にある必要があります。' };
    }
    // On the minus strand the CDS runs backwards, so the base that begins a
    // codon is the one at the HIGH end of the span, not `start`.
    if (feature.strand < 0) {
      var last = removedLength ? start + removedLength - 1 : start - 1;
      offset = positions.indexOf(wrapPosition(last, record.seq.length));
      if (offset === -1) {
        return { error: '選択範囲が指定された CDS ('
          + feature.locationText() + ') からはみ出しています。' };
      }
    }
    if (offset % 3 !== 0) {
      return { error: '選択の開始位置がコドンの境界にありません (CDS の '
        + (offset + 1) + ' 塩基目、残基 ' + (Math.floor(offset / 3) + 1)
        + ' の ' + (offset % 3 + 1) + ' 番目の塩基)。'
        + 'このまま挿入するとアミノ酸配列が読み枠からずれ、'
        + '指定したペプチドにはなりません。'
        + '残基 ' + (Math.floor(offset / 3) + 1) + ' の先頭か、残基 '
        + (Math.floor(offset / 3) + 2) + ' の先頭を選び直してください。' };
    }
    if (removedLength % 3 !== 0) {
      return { error: '置き換える範囲が ' + removedLength
        + ' bp で 3 の倍数ではありません。コドンの境界で始まっていても、'
        + '途中で終わると以降の読み枠がずれます。' };
    }
    return { note: '読み枠を確認しました (CDS の残基 '
      + (Math.floor(offset / 3) + 1) + ' の位置、コドンの境界です)。' };
  }

  /* The whole replacement design: coordinates in, primers and a verified
   * plasmid out. */
  function designReplacement(record, options) {
    options = options || {};
    var start = options.start, end = options.end;
    if (options.selected && (start === undefined || start === null)) {
      var located = locateSelection(record, options.selected);
      start = located.start;
      end = located.end;
      options.located = located;
    }
    if (start === undefined || start === null) {
      throw PbError('置き換える領域が指定されていません');
    }
    if (end === undefined || end === null) { end = start - 1; }

    var notes = [];
    var newSequence = (options.newSequence || '').replace(/\s/g, '');
    var codonChoice = null;
    if (options.asProtein) {
      var bounds = spanBounds(record, start, end);
      var removedLength = bounds.removedLength;
      var frame = checkReadingFrame(record, options.cds, start, removedLength);
      if (frame.error) { throw PbError(frame.error); }
      if (frame.note) { notes.push(frame.note); }

      var peptide = newSequence.toUpperCase();
      if (peptide) {
        var cdsFeature = asCdsFeature(options.cds);
        var cdsStrand = cdsFeature ? cdsFeature.strand : 1;
        codonChoice = encodePeptideForSite(record, bounds, peptide, {
          host: options.host,
          allowRare: options.allowRare, avoid: options.avoid,
          minUsage: options.minUsage,
          codonCandidates: options.codonCandidates,
          strand: cdsStrand
        });
        newSequence = codonChoice.dna;
        // The DNA about to be inserted must read back as the peptide asked
        // for. This cannot fail while the code is right; it is here because it
        // silently did fail for two weeks.
        var reading = cdsStrand > 0 ? newSequence : reverseComplement(newSequence);
        var readsAs = translate(reading).protein;
        if (readsAs !== peptide) {
          throw PbError('内部エラー: 変換した塩基配列が指定のアミノ酸配列に'
            + 'なりません (指定 ' + peptide + ' / 変換後の翻訳 ' + readsAs
            + ')。報告してください。');
        }
        notes.push('アミノ酸 ' + peptide.length + ' 残基を ' + newSequence.length
                   + ' bp に変換しました (' + codonChoice.level + ')。');
      }
      if (removedLength % 3 !== 0) {
        notes.push('選択範囲が ' + removedLength + ' bp で 3 の倍数ではありません。'
          + 'アミノ酸で指定する場合、読み枠がずれる可能性があります。');
      }
    }

    var site = replacementSite(record, start, end, newSequence);
    var maxPrimerLen = options.maxPrimerLen || DEFAULT_MAX_PRIMER_LEN;
    var minAnnealLen = options.minAnnealLen || DEFAULT_MIN_ANNEAL_LEN;
    var maxAnnealLen = options.maxAnnealLen || DEFAULT_MAX_ANNEAL_LEN;
    var budget = tailBudgetNt(maxPrimerLen, minAnnealLen);
    var problems = [];
    if (site.insertLength > budget) {
      problems.push('新しい配列が ' + site.insertLength + ' bp で、5\' テールの上限 '
        + budget + ' nt を超えています。1 ラウンドの逆PCRでは入りません。');
      return { record: record, site: site, pair: null, alternatives: [],
               verification: null, diagnostics: [], adoptedDiagnostic: null,
               codonChoice: codonChoice, budget: budget, notes: notes,
               problems: problems, ok: false,
               constructName: options.constructName || record.name };
    }

    function search(forSite) {
      return designPrimerPairs(record, forSite, {
        maxPrimerLen: maxPrimerLen, minAnnealLen: minAnnealLen,
        maxAnnealLen: maxAnnealLen, limit: options.candidates || 5
      });
    }

    var pairs = search(site);
    /* Every string in codonChoice.alternatives is inside the usage band, so
     * swapping between them costs nothing in expression terms while the
     * primer may be measurably better in one of them. Until 2026-09-07 the
     * codons were fixed before the primer search ever ran, so that freedom
     * went unused. Expression still comes first: only strings the band
     * already accepted are tried here. */
    if (codonChoice && codonChoice.alternatives
        && codonChoice.alternatives.length > 1) {
      var bestScore = pairs.length ? pairs[0].score : null;
      codonChoice.alternatives.slice(1).forEach(function (candidateDna) {
        if (candidateDna === site.insert) { return; }
        var trialSite = replacementSite(record, start, end, candidateDna);
        if (trialSite.insertLength > budget) { return; }
        var trial = search(trialSite);
        if (!trial.length) { return; }
        if (bestScore === null || trial[0].score < bestScore) {
          bestScore = trial[0].score;
          pairs = trial;
          site = trialSite;
        }
      });
      if (site.insert !== codonChoice.dna) {
        notes.push('同じ使用頻度の帯の中で、プライマーがより良くなるコドンに'
                   + '差し替えました: ' + codonChoice.dna + ' → ' + site.insert);
      }
    }
    if (!pairs.length) {
      problems.push('指定した長さの範囲でプライマーを構成できません。');
      return { record: record, site: site, pair: null, alternatives: [],
               verification: null, diagnostics: [], adoptedDiagnostic: null,
               codonChoice: codonChoice, budget: budget, notes: notes,
               problems: problems, ok: false,
               constructName: options.constructName || record.name };
    }
    var best = pairs[0];
    if (!best.usable) {
      problems.push('最良候補でも鋳型上に結合部位が 1 箇所になりません'
                    + '（非特異結合のおそれ）。');
    }
    var verification = verifyReplacement(record, site, best);

    return {
      record: record, site: site, pair: best, alternatives: pairs.slice(1),
      verification: verification, codonChoice: codonChoice, budget: budget,
      notes: notes, problems: problems,
      annealingTimeSec: annealingTimeForPair(best),
      start: start, end: end,
      located: options.located || null,
      constructName: options.constructName || record.name,
      forwardName: 'F', reverseName: 'R',
      ok: problems.length === 0 && best.usable && verification.ok
    };
  }

  // =====================================================================
  // Pipeline: edits -> primers -> verified construct
  // Nothing here calculates; it calls the phases in order.
  // =====================================================================

  function resolveEdits(protein, mutations, targetProtein) {
    if (!!mutations === !!targetProtein) {
      throw PbError('変異表記と目的タンパク質配列は、どちらか一方を指定してください');
    }
    var edits = mutations
      ? parseEdits(mutations)
      : diffProteins(protein.replace(/\*+$/, ''),
                     targetProtein.replace(/\*+$/, ''));
    validateEdits(edits, protein);
    return edits;
  }

  function designRun(record, feature, options) {
    options = options || {};
    var cdsDna = record.featureSequence(feature);
    var translated = translate(cdsDna);
    var protein = translated.protein;
    var notes = [];
    if (translated.trailing) {
      notes.push('CDS の長さが 3 の倍数ではありません（余り '
        + translated.trailing + ' 塩基）。CDS の範囲を確認してください。');
    }

    var edits = resolveEdits(protein, options.mutations, options.targetProtein);
    var maxPrimerLen = options.maxPrimerLen || DEFAULT_MAX_PRIMER_LEN;
    var minAnnealLen = options.minAnnealLen || DEFAULT_MIN_ANNEAL_LEN;
    var maxAnnealLen = options.maxAnnealLen || DEFAULT_MAX_ANNEAL_LEN;
    var budget = tailBudgetNt(maxPrimerLen, minAnnealLen);
    var groups = groupEdits(protein, edits, {
      budget: budget, host: options.host,
      cdsDna: cdsDna, allowRare: options.allowRare, avoid: options.avoid,
      minUsage: options.minUsage
    });

    var rounds = [];
    groups.forEach(function (group, index) {
      var problems = [];
      if (!group.feasible(budget)) {
        problems.push('必要テール長 ' + group.requiredTailNt + ' nt が上限 '
          + budget + ' nt を超えています');
        rounds.push({ index: index + 1, group: group, site: null, pair: null,
                      alternatives: [], verification: null, diagnostics: [],
                      adoptedDiagnostic: null, problems: problems, ok: false });
        return;
      }
      var site = editSite(record, feature, group);
      var diagnostics = findDiagnosticSites(record, feature, group, site, {
        enzymes: options.enzymes, host: options.host,
      });
      var adopted = null;
      if (options.diagnosticEnzyme) {
        diagnostics.forEach(function (candidate) {
          if (!adopted && candidate.enzyme.toLowerCase()
              === options.diagnosticEnzyme.toLowerCase()) { adopted = candidate; }
        });
        if (!adopted) {
          problems.push('診断用サイト ' + options.diagnosticEnzyme
                        + ' はこのラウンドでは導入できません');
        } else {
          group.newDna = adopted.newDna;
          site = editSite(record, feature, group);
        }
      }
      var pairs = designPrimerPairs(record, site, {
        maxPrimerLen: maxPrimerLen, minAnnealLen: minAnnealLen,
        maxAnnealLen: maxAnnealLen, limit: options.candidates || 5
      });
      if (!pairs.length) {
        problems.push('指定した長さの範囲でプライマーを構成できません');
        rounds.push({ index: index + 1, group: group, site: site, pair: null,
                      alternatives: [], verification: null,
                      diagnostics: diagnostics, adoptedDiagnostic: adopted,
                      problems: problems, ok: false });
        return;
      }
      var best = pairs[0];
      if (!best.usable) {
        problems.push('最良候補でも鋳型上に結合部位が 1 箇所になりません'
                      + '（非特異結合のおそれ）');
      }
      var verification = verifyRound(record, feature, protein, group, best,
                                     record.name);
      var label = 'R' + (index + 1);
      rounds.push({
        index: index + 1, label: label,
        forwardName: label + '_F', reverseName: label + '_R',
        group: group, site: site, pair: best, alternatives: pairs.slice(1),
        verification: verification, diagnostics: diagnostics,
        adoptedDiagnostic: adopted, problems: problems,
        annealingTimeSec: annealingTimeForPair(best),
        ok: problems.length === 0 && best.usable && verification.ok
      });
    });

    return {
      record: record, feature: feature, protein: protein, edits: edits,
      groups: groups, rounds: rounds, budget: budget,
      maxPrimerLen: maxPrimerLen, minAnnealLen: minAnnealLen,
      maxAnnealLen: maxAnnealLen,
      constructName: options.constructName || record.name,
      notes: notes,
      newProtein: applyEdits(protein, edits),
      ok: rounds.length > 0 && rounds.every(function (r) { return r.ok; })
    };
  }

  function primerRows(run) {
    var rows = [];
    run.rounds.forEach(function (round) {
      if (!round.pair) { return; }
      [[round.forwardName, round.pair.forward],
       [round.reverseName, round.pair.reverse]].forEach(function (entry) {
        rows.push({ name: run.constructName + '_' + entry[0],
                    sequence: entry[1].sequence, length: entry[1].length,
                    tm: entry[1].tm, gc: entry[1].gcPercent,
                    tail: entry[1].tail.length, anneal: entry[1].anneal.length });
      });
    });
    return rows;
  }

  return {
    // constants
    GENETIC_CODE: GENETIC_CODE, STOP_CODON: STOP_CODON,
    TM_FORMULA_MAX_LEN: TM_FORMULA_MAX_LEN,
    MAX_GC_RUN: MAX_GC_RUN, MAX_HOMOPOLYMER_RUN: MAX_HOMOPOLYMER_RUN,
    DEFAULT_MAX_PRIMER_LEN: DEFAULT_MAX_PRIMER_LEN,
    DEFAULT_MIN_ANNEAL_LEN: DEFAULT_MIN_ANNEAL_LEN,
    DEFAULT_MAX_ANNEAL_LEN: DEFAULT_MAX_ANNEAL_LEN,
    THREE_PRIME_CHECK_LEN: THREE_PRIME_CHECK_LEN,
    PENALTY_NON_UNIQUE: PENALTY_NON_UNIQUE,
    ECOLI_RARE_CODONS: ECOLI_RARE_CODONS,
    ECOLI_CODON_PREFERENCE: ECOLI_CODON_PREFERENCE,
    BUILTIN_ENZYMES: BUILTIN_ENZYMES,
    LOCUS_TOPOLOGY_COL: LOCUS_TOPOLOGY_COL, LOCUS_DATE_COL: LOCUS_DATE_COL,
    GENBANK_LINE_WIDTH: GENBANK_LINE_WIDTH,
    ANNEAL_TEMP_C: ANNEAL_TEMP_C, DENATURE_TEMP_C: DENATURE_TEMP_C,
    DENATURE_TIME_SEC: DENATURE_TIME_SEC, PCR_CYCLES: PCR_CYCLES,
    EXTENSION_TEMP_C_3STEP: EXTENSION_TEMP_C_3STEP,
    EXTENSION_TEMP_C_2STEP: EXTENSION_TEMP_C_2STEP,
    PRIMER_AMOUNT_PMOL: PRIMER_AMOUNT_PMOL,
    // sequence
    Record: Record, Feature: Feature,
    reverseComplement: reverseComplement, translate: translate,
    countOccurrences: countOccurrences, gcPercent: gcPercent,
    longestHomopolymerRun: longestHomopolymerRun, longestGcRun: longestGcRun,
    wrapPosition: wrapPosition,
    // io
    detectFormat: detectFormat, parseLocation: parseLocation,
    parseGenbank: parseGenbank, parsePlain: parsePlain,
    requireDna: requireDna,
    EDIT_MARKER_LABEL: EDIT_MARKER_LABEL,
    positionsToSegments: positionsToSegments,
    loadFromText: loadFromText, formatGenbank: formatGenbank,
    formatLocusLine: formatLocusLine, formatQualifier: formatQualifier,
    formatOrigin: formatOrigin, genbankDate: genbankDate,
    // thermo
    tmTakara: tmTakara, tmFormulaApplies: tmFormulaApplies,
    annealingTimeSec: annealingTimeSec, annealingTimeNote: annealingTimeNote,
    extensionTimeSec: extensionTimeSec,
    practicalExtensionTimeSec: practicalExtensionTimeSec,
    formatSeconds: formatSeconds,
    // codons and edits
    getCodonPreference: getCodonPreference,
    chooseCodons: chooseCodons, Edit: Edit, parseEdits: parseEdits,
    validateEdits: validateEdits, applyEdits: applyEdits,
    alignProteins: alignProteins, diffProteins: diffProteins,
    tailBudgetNt: tailBudgetNt, buildEditGroup: buildEditGroup,
    groupEdits: groupEdits, cdsPositionMap: cdsPositionMap,
    residueToPlasmid: residueToPlasmid,
    // design and verify
    editSite: editSite, designPrimerPairs: designPrimerPairs,
    splitPoint: splitPoint, countBindingSites: countBindingSites,
    annealingTimeForPair: annealingTimeForPair,
    pcrProduct: pcrProduct, verifyRound: verifyRound,
    spanBounds: spanBounds,
    replacementSite: replacementSite, locateSelection: locateSelection,
    checkReadingFrame: checkReadingFrame,
    encodePeptideForSite: encodePeptideForSite,
    verifyReplacement: verifyReplacement,
    designReplacement: designReplacement,
    parseEnzymeTable: parseEnzymeTable, placeSite: placeSite,
    findDiagnosticSites: findDiagnosticSites,
    // pipeline
    designRun: designRun, primerRows: primerRows
  };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = PB; }

/* The Markdown design report, downloaded alongside the .ape. */
PB.designMarkdown = function (run, dateText) {
  var lines = [];
  var stamp = dateText || (function () {
    var d = new Date();
    var m = d.getMonth() + 1, day = d.getDate();
    return d.getFullYear() + '-' + (m < 10 ? '0' + m : m) + '-'
      + (day < 10 ? '0' + day : day);
  })();
  var strip = function (p) { return p.replace(/\*+$/, ''); };

  lines.push('# ' + run.constructName + ' 設計レポート', '');
  lines.push('- 作成日: ' + stamp);
  lines.push("- 手法: インバースPCR（T4 DNA ligase で環状化）");
  lines.push('- 鋳型: ' + run.record.name + '  ' + run.record.length() + ' bp  '
             + run.record.topology());
  lines.push('- CDS: ' + run.feature.locationText() + '  '
             + strip(run.protein).length + ' 残基 + 終止', '');

  if (run.notes.length) {
    lines.push('## 注意', '');
    run.notes.forEach(function (n) { lines.push('- ' + n); });
    lines.push('');
  }

  lines.push('## 編集内容', '', '| # | 指定 | 内容 |', '| --- | --- | --- |');
  run.edits.forEach(function (edit, index) {
    lines.push('| ' + (index + 1) + ' | `' + (edit.source || '差分') + '` | '
               + edit.describe() + ' |');
  });
  lines.push('', '改変後のタンパク質: ' + strip(run.protein).length + ' 残基 → '
             + strip(run.newProtein).length + ' 残基', '');

  lines.push('## プライマー', '');
  lines.push('| 名前 | 配列 (5′→3′) | 長さ | Tm | GC% |',
             '| --- | --- | --- | --- | --- |');
  PB.primerRows(run).forEach(function (row) {
    lines.push('| ' + row.name + ' | `' + row.sequence + '` | ' + row.length
      + ' nt | ' + row.tm + '℃ | ' + row.gc.toFixed(0) + '% |');
  });
  lines.push('', 'Tm は 3′ アニーリング領域について、タカラバイオ簡易式 '
             + '`Tm = 2×(A+T) + 4×(G+C) − 5` で算出。', '');

  run.rounds.forEach(function (round) {
    lines.push('## ラウンド ' + round.index, '');
    if (!round.pair) {
      round.problems.forEach(function (p) { lines.push('- **' + p + '**'); });
      lines.push('');
      return;
    }
    var pair = round.pair, site = pair.site;
    lines.push('編集: ' + round.group.edits.map(function (e) {
      return '`' + (e.source || '差分') + '`';
    }).join(', '), '');
    lines.push('| 項目 | 値 |', '| --- | --- |');
    lines.push('| 除去する領域 | プラスミド ' + site.leftPos + ' の直後から '
      + site.rightPos + ' の直前 (' + site.removedLength + ' bp) |');
    lines.push('| 挿入する配列 | ' + site.insertLength + ' bp `' + site.insert
      + '` |');
    lines.push('| テールの分配 | R に ' + pair.split + ' nt / F に '
      + (site.insertLength - pair.split) + ' nt |');
    lines.push('| Tm 差 | ' + pair.tmDifference + '℃ |');
    lines.push('| 予想産物長 | ' + pair.productLength + ' bp |', '');

    /* A blank form, not our numbers.
     *
     * Until 2026-09-25 this printed the temperatures, times and cycle count,
     * and a step-by-step protocol underneath. The user asked for both to go:
     * they know the protocol, and they want to work the conditions out
     * themselves. What the tool owes them is the primers. */
    lines.push('### PCR 条件（記入欄）', '');
    lines.push('| 段階 | 温度 | 時間 |', '| --- | --- | --- |');
    lines.push('| 熱変性 |  |  |');
    lines.push('| アニーリング |  |  |');
    lines.push('| 伸長 |  |  |');
    lines.push('| サイクル数 |  |  |', '');

    lines.push('### 検算', '');
    var v = round.verification;
    if (!v) { lines.push('- 実施できませんでした。'); }
    else if (v.ok) {
      lines.push('プライマー配列から環状化後のプラスミドを組み立て直し、'
        + 'CDS を配列検索で特定して翻訳したところ、'
        + '**目的のタンパク質と完全に一致しました**。', '');
      lines.push('| 項目 | 値 |', '| --- | --- |');
      lines.push('| 環状化後の長さ | ' + v.newRecord.length() + ' bp |');
      lines.push('| 改変後の CDS | ' + v.newCds.locationText() + ' |');
      lines.push('| 翻訳結果 | ' + strip(v.observedProtein).length
                 + ' 残基（内部終止コドンなし） |');
    } else {
      lines.push('**検算で問題が見つかりました。**');
      v.problems.forEach(function (p) { lines.push('- ' + p); });
    }
    if (v) { v.notes.forEach(function (n) { lines.push('- 注意: ' + n); }); }
    lines.push('');

    if (round.diagnostics.length) {
      lines.push('### 診断用制限酵素サイトの候補', '');
      lines.push('親プラスミドに存在せず、改変後に 1 箇所だけ現れる酵素です。'
                 + 'アミノ酸配列は変えません。', '');
      lines.push('| 酵素 | 認識配列 | 変更する塩基数 | 採用 |',
                 '| --- | --- | --- | --- |');
      round.diagnostics.forEach(function (c) {
        var adopted = round.adoptedDiagnostic
          && round.adoptedDiagnostic.enzyme === c.enzyme;
        lines.push('| ' + c.enzyme + ' | `' + c.recognition + '` | '
                   + c.changedBases + ' | ' + (adopted ? '採用' : '-') + ' |');
      });
      lines.push('');
    }

    if (round.alternatives.length) {
      lines.push('### 代替候補', '');
      lines.push('| 順位 | F 配列 | R 配列 | Tm 差 | スコア |',
                 '| --- | --- | --- | --- | --- |');
      round.alternatives.forEach(function (a, index) {
        lines.push('| ' + (index + 2) + ' | `' + a.forward.sequence + '` | `'
          + a.reverse.sequence + '` | ' + a.tmDifference + '℃ | '
          + a.score.toFixed(1) + ' |');
      });
      lines.push('', 'スコアは減点方式で、小さいほど良い候補です。', '');
    }
  });

  lines.push('## 出典', '');
  lines.push('- Tm: タカラバイオ PrimeSTAR HS 製品基礎情報'
    + '（catalog.takara-bio.co.jp, unitid=U100004616, 2026-08-25 確認）');
  lines.push('- 遺伝暗号: NCBI 遺伝暗号表 1');
  lines.push('- 制限酵素の認識配列: 一般的な供給元カタログの記載（発注前に確認してください）',
             '');
  lines.push('このレポートは設計支援の結果です。最終確認は個人で'
             + '行ってください。');
  return lines.join('\n') + '\n';
};

/* COMMENT lines recording what a round did, in ApE's style, appended after
 * the parent's own header so the construction history is not lost. */
PB.editCommentLines = function (run, round) {
  var lines = ['COMMENT     Parent: ' + run.record.name,
    "COMMENT     Method: inverse_PCR (blunt-end ligation with T4 ligase)"];
  round.group.edits.forEach(function (edit) {
    lines.push('COMMENT     Edit: ' + (edit.source || edit.describe()));
  });
  lines.push('COMMENT     Primer_F: ' + round.pair.forward.sequence);
  lines.push('COMMENT     Primer_R: ' + round.pair.reverse.sequence);
  if (round.adoptedDiagnostic) {
    lines.push('COMMENT     Diagnostic_site: ' + round.adoptedDiagnostic.enzyme);
  }
  return lines;
};

/* Report for the replacement mode. The check it reports is at the DNA level:
 * the circle rebuilt from the primers is compared with the circle the request
 * describes, computed independently from the template. */
PB.replacementMarkdown = function (result, dateText) {
  var lines = [];
  var stamp = dateText || (function () {
    var d = new Date(), m = d.getMonth() + 1, day = d.getDate();
    return d.getFullYear() + '-' + (m < 10 ? '0' + m : m) + '-'
      + (day < 10 ? '0' + day : day);
  })();
  var site = result.site, pair = result.pair, v = result.verification;

  lines.push('# ' + result.constructName + ' 設計レポート（領域の置き換え）', '');
  lines.push('- 作成日: ' + stamp);
  lines.push("- 手法: インバースPCR（T4 DNA ligase で環状化）");
  lines.push('- 鋳型: ' + result.record.name + '  ' + result.record.length()
             + ' bp  ' + result.record.topology(), '');

  lines.push('## 置き換える内容', '', '| 項目 | 値 |', '| --- | --- |');
  if (site.removedLength) {
    lines.push('| 選択した領域 | ' + result.start + '..' + result.end + '（'
               + site.removedLength + ' bp） |');
    lines.push('| 元の配列 | `' + site.removed + '` |');
  } else {
    lines.push('| 挿入位置 | ' + site.leftPos + ' と ' + site.rightPos
               + ' の間（除去なし） |');
  }
  lines.push('| 新しい配列 | ' + (site.insertLength
    ? site.insertLength + ' bp `' + site.insert + '`' : '（なし＝欠失）') + ' |');
  if (result.codonChoice) {
    lines.push('| アミノ酸から変換 | ' + result.codonChoice.peptide + '（'
               + result.codonChoice.level + '） |');
  }
  lines.push('| 予想産物長 | ' + pair.productLength + ' bp |', '');
  result.notes.forEach(function (note) { lines.push('- ' + note); });
  if (result.notes.length) { lines.push(''); }

  lines.push('## プライマー', '');
  lines.push('| 名前 | 配列 (5′→3′) | 長さ | Tm | GC% |',
             '| --- | --- | --- | --- | --- |');
  [[result.constructName + '_F', pair.forward],
   [result.constructName + '_R', pair.reverse]].forEach(function (entry) {
    lines.push('| ' + entry[0] + ' | `' + entry[1].sequence + '` | '
      + entry[1].length + ' nt | ' + entry[1].tm + '℃ | '
      + entry[1].gcPercent.toFixed(0) + '% | 必要 |');
  });
  lines.push('', 'Tm は 3′ アニーリング領域について、タカラバイオ簡易式 '
             + '`Tm = 2×(A+T) + 4×(G+C) − 5` で算出。'
             + 'テールの分配は R に ' + pair.split + ' nt / F に '
             + (site.insertLength - pair.split) + ' nt。', '');

  // A blank form, not our numbers -- see the note in designMarkdown.
  lines.push('## PCR 条件（記入欄）', '');
  lines.push('| 段階 | 温度 | 時間 |', '| --- | --- | --- |');
  lines.push('| 熱変性 |  |  |');
  lines.push('| アニーリング |  |  |');
  lines.push('| 伸長 |  |  |');
  lines.push('| サイクル数 |  |  |', '');

  lines.push('## 検算', '');
  if (v && v.ok) {
    lines.push('プライマーの配列から環状化後のプラスミドを組み立て直し、'
      + '**目的の配列と塩基 1 個まで一致すること**を確認しました。'
      + '選択範囲の外は 1 塩基も変化していません。', '');
    lines.push('| 項目 | 値 |', '| --- | --- |');
    lines.push('| 環状化後の長さ | ' + v.newRecord.length() + ' bp |');
    if (v.insertStart) {
      lines.push('| 新しい配列の位置 | ' + v.insertStart + '..' + v.insertEnd
                 + ' |');
    }
  } else {
    lines.push('**検算で問題が見つかりました。**');
    (result.problems.concat(v ? v.problems : [])).forEach(function (p) {
      lines.push('- ' + p);
    });
  }
  if (v) { v.notes.forEach(function (n) { lines.push('- 注意: ' + n); }); }
  lines.push('');

  if (result.alternatives.length) {
    lines.push('## 代替候補', '');
    lines.push('| 順位 | F 配列 | R 配列 | Tm 差 | スコア |',
               '| --- | --- | --- | --- | --- |');
    result.alternatives.forEach(function (a, index) {
      lines.push('| ' + (index + 2) + ' | `' + a.forward.sequence + '` | `'
        + a.reverse.sequence + '` | ' + a.tmDifference + '℃ | '
        + a.score.toFixed(1) + ' |');
    });
    lines.push('');
  }

  lines.push('## 出典', '');
  lines.push('- Tm: タカラバイオ PrimeSTAR HS 製品基礎情報'
    + '（catalog.takara-bio.co.jp, unitid=U100004616, 2026-08-25 確認）');
  lines.push('- 遺伝暗号: NCBI 遺伝暗号表 1', '');
  lines.push('このレポートは設計支援の結果です。最終確認は個人で'
             + '行ってください。');
  return lines.join('\n') + '\n';
};

PB.replacementCommentLines = function (result) {
  var lines = ['COMMENT     Parent: ' + result.record.name,
    "COMMENT     Method: inverse_PCR (blunt-end ligation with T4 ligase)"];
  if (result.site.removedLength) {
    lines.push('COMMENT     Replaced: ' + result.start + '..' + result.end
               + ' (' + result.site.removedLength + ' bp)');
  } else {
    lines.push('COMMENT     Inserted_between: ' + result.site.leftPos + ' and '
               + result.site.rightPos);
  }
  lines.push('COMMENT     New_sequence: '
             + (result.site.insert || '(deleted)'));
  if (result.codonChoice) {
    lines.push('COMMENT     From_peptide: ' + result.codonChoice.peptide);
  }
  lines.push('COMMENT     Primer_F: ' + result.pair.forward.sequence);
  lines.push('COMMENT     Primer_R: ' + result.pair.reverse.sequence);
  return lines;
};
