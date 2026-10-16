
/* ===================================================================
 * Primer Builder: the interface.
 *
 * The flow is: read the plasmid, choose how the edit is stated, then design.
 * There are two ways to state it. This unit owns the first -- look at the
 * sequence and select the range, with coordinates never typed in -- and hands
 * the second, mutation notation, to PBMUT in mutation-ui.js. Loading, the
 * conditions, the request form and the self-test are shared by both.
 *
 * Primer design itself is the engine's job (designReplacement / designRun).
 * =================================================================== */
(function () {
  'use strict';

  var state = {
    record: null, text: '', filename: '',
    window: null, orf: null, orfChoice: 0, orfCandidates: [],
    protein: '', map: null,
    // plasmid position -> residue number, for painting the protein track
    // when the selection was made on the DNA track
    residueAt: null,
    selection: null,        // {mode:'dna'|'protein', start, end, resStart, resEnd}
    dragging: false, moved: false, anchor: null, pending: null,
    result: null,
    // Protein Builder sequences kept in this browser, and the ones currently
    // coloured in the two sequence tracks. Several may be coloured at once;
    // they share one colour and are told apart by name on hover.
    pbuild: [], pbHighlights: []
  };

  function esc(t) {
    return String(t === undefined || t === null ? '' : t)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }
  function $(id) { return document.getElementById(id); }
  function show(id, on) { $(id).style.display = on ? '' : 'none'; }
  function banner(kind, title, body) {
    return '<div class="banner ' + kind + '"><strong>' + esc(title) + '</strong>'
      + (body || '') + '</div>';
  }
  var SAY_TITLE = { ok: '完了', warn: '注意' };
  function say(text, kind) {
    $('messages').innerHTML = text
      ? banner(kind || 'bad',
               SAY_TITLE[kind] || '入力を確認してください', esc(text))
      : '';
  }
  function copyBtn(text, label) {
    return '<button type="button" class="tiny" data-copy="' + esc(text) + '">'
      + esc(label || 'コピー') + '</button>';
  }
  function num(id, fallback) {
    var raw = ($(id).value || '').trim();
    if (!raw) { return fallback; }
    if (!/^\d+$/.test(raw)) {
      throw new Error($(id).getAttribute('data-label') + ' は整数で入力してください');
    }
    return parseInt(raw, 10);
  }

  // ---- 1. load ---------------------------------------------------------

  /* Once a sequence is in, the paste box has done its job.
   *
   * It used to stay open under the drop zone after a file was loaded, which
   * reads as if the file had not been accepted (2026-09-25). The button
   * brings it back for someone who wants to paste a different sequence. */
  function showPasteArea(on) {
    show('pasteArea', on);
    show('loadPasted', on);
    show('showPaste', !on);
  }

  function readFile(file) {
    var reader = new FileReader();
    reader.onload = function (e) {
      state.text = String(e.target.result);
      state.filename = file.name;
      $('filename').textContent = file.name;
      $('pasted').value = '';
      loadSequence();
    };
    reader.onerror = function () { say('ファイルを読めませんでした。'); };
    reader.readAsText(file);
  }

  function loadSequence() {
    var text = state.text || $('pasted').value;
    if (!text.trim()) {
      say('配列を指定してください。ファイルをドラッグするか、貼り付けてください。');
      return;
    }
    try {
      state.record = PB.loadFromText(text, state.filename || 'pasted');
    } catch (error) {
      state.record = null;
      say(error.message);
      hideAllSteps();
      return;
    }
    state.selection = null; state.result = null;
    state.orfChoice = 0;
    $('result').innerHTML = '';
    $('messages').innerHTML = state.record.warnings.map(function (w) {
      return banner('warn', '注意', esc(w));
    }).join('');

    var record = state.record;
    $('summary').innerHTML = '<table>'
      + '<tr><th>名称</th><td>' + esc(record.name) + '</td></tr>'
      + '<tr><th>長さ</th><td>' + record.length() + ' bp</td></tr>'
      + '<tr><th>トポロジー</th><td>' + esc(record.topology()) + '</td></tr>'
      + '<tr><th>フィーチャー</th><td>' + record.features.length + ' 件</td></tr>'
      + '</table>';
    $('load').classList.add('done');
    showPasteArea(false);
    renderViews();
    PBMUT.setRecord(record);
    applyMode();
    renderSelection();
  }

  // ---- which way the edit is stated ------------------------------------
  //
  // Both modes end in the same primer search; they differ only in how the
  // change is described. Switching clears the result, because a result from
  // the other mode would look as if it belonged to this one.

  // The folds inside viewSection are not listed here: showing or hiding the
  // view already covers them.
  var MODE_STEPS = ['viewSection', 'assembleStep', 'editSection',
                    'mnCds', 'mnEdit', 'designStep'];

  function hideAllSteps() {
    MODE_STEPS.forEach(function (id) { show(id, false); });
  }

  function isMutationMode() {
    var checked = document.querySelector('input[name=mode]:checked');
    return !!checked && checked.value === 'mutation';
  }

  function applyMode() {
    var mutation = isMutationMode(), loaded = !!state.record;
    ['viewSection', 'editSection'].forEach(function (id) {
      show(id, loaded && !mutation);
    });
    // Only the combined build has something to put in it (assemble.js).
    show('assembleStep', loaded && !mutation && !!window.PBASM);
    show('mnCds', loaded && mutation);
    show('mnEdit', loaded && mutation);
    show('designStep', loaded);
    show('diagnosticWrap', mutation);
    $('designHint').textContent = mutation
      ? '残基番号は上の翻訳の番号です。離れた変異は自動でラウンドが分かれます。'
      : '置き換える範囲は、上の配列ビューで選んでください。';
    state.result = null;
    $('result').innerHTML = '';
    $('textPanel').innerHTML = '';
    PBMUT.clearResult();
  }

  // ---- 2. views --------------------------------------------------------

  function renderViews() {
    var record = state.record;
    var promoter = ($('t7promoter').value || '').toUpperCase().replace(/\s/g, '')
      || GEVIEW.T7_PROMOTER;
    var terminator = ($('t7terminator').value || '').toUpperCase()
      .replace(/\s/g, '') || GEVIEW.T7_TERMINATOR;
    var w;
    try { w = GEVIEW.expressionWindow(record, promoter, terminator); }
    catch (error) { w = { found: false }; }
    var previousOrf = state.orf ? state.orf.locationText() : null;
    state.window = w;
    state.orf = null; state.protein = ''; state.map = null;
    state.residueAt = null;

    var html = '', chooserHtml = '';
    if (w.found) {
      var marks = [
        { start: w.promoter.start, end: w.promoter.end, cls: 'prom' },
        { start: w.terminator.start, end: w.terminator.end, cls: 'term' }
      ];
      var positions = GEVIEW.windowPositions(record, w.from, w.to);
      /* State the range actually on screen, not the rule used to derive it.
       * The old wording always claimed "flank bases before the promoter
       * through flank bases after the terminator", but expressionWindow falls
       * back to the whole plasmid when that span would wrap onto itself --
       * which is what the bundled example does -- so the sentence was simply
       * wrong there. */
      html += '<div class="banner ok"><strong>T7 promoter と terminator を'
        + '見つけました</strong>'
        + 'promoter ' + w.promoter.start + '..' + w.promoter.end
        + (w.promoter.strand > 0 ? '' : '（逆方向鎖）') + '／'
        + 'terminator ' + w.terminator.start + '..' + w.terminator.end
        + '<br>表示範囲 ' + w.from + '..' + w.to
        + '（' + positions.length + ' bp'
        + (w.wholePlasmid ? '・プラスミド全長' : '') + '）'
        + '／色: <span class="promkey">promoter</span> '
        + '<span class="termkey">terminator</span></div>';
      html += '<div class="seqview" id="dnaWindow">'
        + GEVIEW.renderDna(record, positions, marks) + '</div>';

      /* Which ORF is on screen decides every residue number, so it is not
       * decided silently. A CDS the file declares wins outright; otherwise
       * the candidates go on screen and the longest is merely the default. */
      var declared = GEVIEW.annotatedCds(record, w.between, w.strand);
      var candidates = declared
        ? [] : GEVIEW.orfCandidates(record, w.between, w.strand);
      state.orfCandidates = candidates;
      if (declared) {
        state.orf = declared;
      } else if (candidates.length) {
        if (state.orfChoice >= candidates.length) { state.orfChoice = 0; }
        state.orf = candidates[state.orfChoice];
      } else {
        state.orf = null;
      }
      if (state.orf) {
        var rendered = GEVIEW.renderProtein(record, state.orf);
        state.protein = rendered.protein;
        state.map = rendered.map;
        /* Built once here rather than in paint(), which runs on every mouse
         * move during a drag. map is the CDS positions in reading order, so
         * three consecutive entries are one codon. */
        state.residueAt = {};
        for (var at = 0; at < rendered.map.length; at++) {
          state.residueAt[rendered.map[at]] = Math.floor(at / 3) + 1;
        }
        html += '<h3>promoter と terminator に挟まれたタンパク質</h3>'
          + '<p class="hint">' + state.orf.segments.map(function (s) {
              return s[0] + '..' + s[1];
            }).join(', ')
          + '（' + (state.orf.strand > 0 ? '' : '逆方向鎖、')
          + state.protein.replace(/\*+$/, '').length + ' 残基）</p>'
          + '<div class="seqview" id="proteinWindow">' + rendered.html + '</div>';
        if (declared) {
          chooserHtml = '<p class="hint">読み枠は<strong>ファイルに書かれている '
            + 'CDS 注釈</strong>です（' + esc(declared.label()) + '）。'
            + '注釈があるときは推測より優先します。</p>';
        } else if (candidates.length > 1) {
          chooserHtml = '<label for="orfChoice">開始コドン'
            + '<small>同じ終止コドンを共有する ATG が ' + candidates.length
            + ' 個あります。選び方で残基番号が変わります。'
            + '既定はいちばん長いものです</small></label>'
            + '<select id="orfChoice">';
          candidates.forEach(function (item, index) {
            chooserHtml += '<option value="' + index + '"'
              + (index === state.orfChoice ? ' selected' : '') + '>'
              + item.segments[0][0] + ' から開始（' + item.residues + ' 残基）'
              + '</option>';
          });
          chooserHtml += '</select>';
        }
        if (state.orf.beyondTerminator) {
          html += banner('warn', '終止コドンが terminator の内側にあります',
            'この読み枠は terminator の位置を越えてから終わっています。'
            + '範囲が正しいか確認してください。');
        }
      } else {
        html += banner('warn', 'ORF が見つかりません',
          'promoter と terminator の間に ATG から始まる読み枠が見つかりませんでした。'
          + '塩基配列のほうから範囲を選んでください。');
      }
    } else {
      html += banner('warn', 'T7 promoter / terminator が見つかりません',
        (w.promoter ? '' : 'promoter が見つかりません。')
        + (w.terminator ? '' : 'terminator が見つかりません。')
        + '下の「読み取りの設定」で調整できます。'
        + '全長配列のほうから範囲を選ぶこともできます。');
    }

    $('views').innerHTML = html;
    // The start codon lives with the recognised sequences, not in the middle
    // of the tracks: both decide how the sequence is read (2026-09-27).
    $('orfChooser').innerHTML = chooserHtml;

    // Its own container, below the feature fold rather than inside #views.
    var all = GEVIEW.windowPositions(record, 1, record.length());
    $('fullView').innerHTML = '<details id="fullDetails"><summary>全長配列を表示（'
      + record.length() + ' bp）</summary>'
      + '<div class="seqview" id="dnaFull">'
      + GEVIEW.renderDna(record, all, w.found ? [
          { start: w.promoter.start, end: w.promoter.end, cls: 'prom' },
          { start: w.terminator.start, end: w.terminator.end, cls: 'term' }
        ] : []) + '</div></details>';
    // A protein-track selection is expressed in residue numbers of a
    // particular ORF. If the ORF moved or vanished, those numbers no longer
    // mean anything, so the selection has to go rather than silently point
    // somewhere else.
    var currentOrf = state.orf ? state.orf.locationText() : null;
    if (state.selection && state.selection.mode === 'protein'
        && currentOrf !== previousOrf) {
      state.selection = null;
      state.pending = null;
    }
    var chooser = $('orfChoice');
    if (chooser) {
      chooser.addEventListener('change', function () {
        state.orfChoice = parseInt(chooser.value, 10) || 0;
        renderViews();
        renderSelection();
      });
    }
    wireSelection();
    // The imported list says where each sequence sits in the CDS on screen,
    // so it has to be redrawn whenever the reading frame changes.
    renderLibraryList();
    /* This rebuilt every span, so nothing on screen is painted any more.
     * Repainting here rather than at each call site: changing the start
     * codon used to drop the imported sequence's colour and never bring it
     * back, because that path did not repaint (found 2026-09-27). */
    paint();
  }

  // ---- selection -------------------------------------------------------

  function wireSelection() {
    ['dnaWindow', 'dnaFull'].forEach(function (id) {
      var node = $(id);
      if (node) { attachDna(node); }
    });
    var protein = $('proteinWindow');
    if (protein) { attachProtein(protein); }
    // renderViews rebuilt these, so the hover handlers go on again with them
    ['dnaWindow', 'dnaFull', 'proteinWindow'].forEach(function (id) {
      var node = $(id);
      if (!node) { return; }
      node.addEventListener('mousemove', showTip);
      node.addEventListener('mouseleave', hideTip);
    });
  }

  function posFrom(event) {
    var target = event.target;
    if (!target || !target.getAttribute) { return null; }
    var pos = target.getAttribute('data-pos');
    return pos === null ? null : parseInt(pos, 10);
  }
  function resFrom(event) {
    var target = event.target;
    if (!target || !target.getAttribute) { return null; }
    var res = target.getAttribute('data-res');
    return res === null ? null : parseInt(res, 10);
  }

  /* The gestures a text editor uses, because the edit being described is the
   * same one: put the caret somewhere and type to insert, select a stretch
   * and type to replace it.
   *
   *   click              -> caret at the nearer edge of the base clicked
   *   drag               -> range, anchor..wherever it ended
   *   shift + click      -> range, anchor..clicked
   *   double click       -> range of exactly the base clicked (point mutation)
   *
   * This replaced a click-then-click range gesture on 2026-09-08. A plain
   * click had to become the caret for the model to hold, so extending a
   * selection without dragging moved to shift+click, and picking out a single
   * base moved to a double click.
   */
  function attachRange(node, read, apply, applyCaret) {
    node.addEventListener('mousedown', function (event) {
      var value = read(event);
      if (value === null) { return; }
      event.preventDefault();
      if (event.shiftKey && state.anchor !== null
          && state.anchor !== undefined) {
        state.dragging = false;
        state.moved = false;
        apply(state.anchor, value);
        return;
      }
      state.dragging = true;
      state.moved = false;
      state.anchor = value;
    });
    node.addEventListener('mousemove', function (event) {
      if (!state.dragging) { return; }
      var value = read(event);
      if (value === null) { return; }
      if (value === state.anchor && !state.moved) { return; }
      state.moved = true;
      apply(state.anchor, value);
    });
    node.addEventListener('mouseup', function (event) {
      if (!state.dragging) { return; }
      state.dragging = false;
      if (state.moved) { return; }
      applyCaret(event, read(event));
    });
    node.addEventListener('dblclick', function (event) {
      var value = read(event);
      if (value === null) { return; }
      event.preventDefault();
      state.anchor = value;
      apply(value, value);
    });
  }

  function attachDna(node) {
    attachRange(node, posFrom, setDnaSelection, setDnaCaret);
  }

  /* No caret on the protein track. "Before residue r" is unambiguous in
   * reading order but not in plasmid coordinates -- on a minus-strand CDS the
   * codon's plasmid positions run backwards -- and amino-acid insertion is
   * what the mutation-notation mode's ins@ is for. Selecting residues still
   * replaces them. */
  function attachProtein(node) {
    attachRange(node, resFrom, setProteinSelection, function () {
      state.selection = null;
      state.pending = null;
      paint();
      renderSelection();
    });
  }

  /* Which side of the base was clicked decides which boundary the caret
   * lands on, the way it does in any text field. */
  function setDnaCaret(event, position) {
    if (position === null) { return; }
    var node = event.target, rect = node.getBoundingClientRect();
    var after = rect.width > 0
      && (event.clientX - rect.left) >= rect.width / 2;
    var n = state.record.length();
    var start = after ? position + 1 : position;
    if (start > n) { start = 1; }
    state.selection = { mode: 'caret', start: start };
    setInputKind('dna');
    paint();
    renderSelection();
  }

  function setDnaSelection(a, b) {
    var start = Math.min(a, b), end = Math.max(a, b);
    state.selection = { mode: 'dna', start: start, end: end };
    setInputKind('dna');
    paint();
    renderSelection();
  }

  /* Back to just-loaded: the plasmid stays, everything about the edit goes.
   * Conditions are left alone -- they have their own 規定値に戻す, and losing
   * a tuned primer length because the selection was wrong would be its own
   * annoyance. */
  function resetDesign() {
    state.selection = null;
    state.anchor = null;
    state.pending = null;
    state.result = null;
    $('newSequence').value = '';
    $('result').innerHTML = '';
    $('textPanel').innerHTML = '';
    say('');
    paint();
    renderSelection();
    $('viewSection').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  /* What the field holds is now chosen, not inferred from which track the
   * selection came from (asked for 2026-09-27). Selecting on a track still
   * sets the likely answer, and it can be changed afterwards: a range picked
   * on the DNA can be replaced with amino acids, and the other way round. */
  function isProteinInput() {
    var picked = document.querySelector('input[name=seqkind]:checked');
    return !!picked && picked.value === 'protein';
  }
  function setInputKind(kind) {
    var target = document.querySelector(
      'input[name=seqkind][value="' + kind + '"]');
    if (target) { target.checked = true; }
  }

  var DNA_ONLY = /^[ACGT]+$/;
  var AMINO_ONLY = /^[ACDEFGHIKLMNPQRSTVWY]+$/;

  function strayLetters(text, allowed) {
    var seen = {}, out = [];
    text.split('').forEach(function (letter) {
      if (allowed.test(letter) || seen[letter]) { return; }
      seen[letter] = true;
      out.push(letter);
    });
    return out.join(' ');
  }

  /* Refuse the case that can actually be told apart: letters that the chosen
   * kind has no meaning for. The other direction cannot be caught here --
   * `ATGCAT` is a valid peptide as well as a valid stretch of DNA -- so that
   * one gets a note under the field instead of a refusal. */
  function inputError(text, asProtein) {
    var upper = text.toUpperCase();
    if (!upper) { return ''; }
    if (asProtein) {
      if (AMINO_ONLY.test(upper)) { return ''; }
      return 'アミノ酸として読めない文字が入っています: '
        + strayLetters(upper, /[ACDEFGHIKLMNPQRSTVWY]/)
        + '。1 文字表記で入力してください。'
        + '塩基配列を入れるつもりなら「塩基配列で入力」を選んでください。';
    }
    if (DNA_ONLY.test(upper)) { return ''; }
    return '塩基配列に使えない文字が入っています: '
      + strayLetters(upper, /[ACGT]/)
      + '。A/T/G/C で入力してください。'
      + 'アミノ酸配列を入れるつもりなら「アミノ酸配列で入力」を選んでください。';
  }

  /* The label under the field says whether the edit is an insertion or a
   * replacement, so typing into it changes what the label should read. */
  function wireNewSequenceEcho() {
    $('newSequence').addEventListener('input', function () {
      if (state.selection) { renderSelection(); }
    });
    Array.prototype.forEach.call(
      document.querySelectorAll('input[name=seqkind]'), function (input) {
        input.addEventListener('change', function () { renderSelection(); });
      });
  }

  /* Paste a stretch of the plasmid and have it found, rather than counting
   * to it. The engine has done this all along (`locateSelection`, which also
   * checks the reverse strand and refuses an ambiguous match); it was the
   * screen that had no way in, since the 2026-08-28 merge. */
  function findPasted() {
    if (!state.record) { say('先に配列を読み込んでください。'); return; }
    var text = $('findText').value || '';
    if (!text.trim()) { say('探す配列を貼り付けてください。'); return; }
    var found;
    try { found = PB.locateSelection(state.record, text); }
    catch (error) { say(error.message); return; }
    state.selection = { mode: 'dna', start: found.start, end: found.end };
    state.anchor = null;
    state.pending = null;
    setInputKind('dna');
    paint();
    renderSelection();
    say('見つけました: ' + found.start + '..' + found.end + '（'
        + found.length + ' bp'
        + (found.onReverseStrand ? '、逆方向鎖と一致' : '')
        + (found.crossesOrigin ? '、原点をまたぎます' : '') + '）', 'ok');
  }

  function setProteinSelection(a, b) {
    var first = Math.min(a, b), last = Math.max(a, b);
    var codonFirst = PB.residueToPlasmid(state.map, first);
    var codonLast = PB.residueToPlasmid(state.map, last);
    var all = codonFirst.concat(codonLast);
    var start = Math.min.apply(null, all), end = Math.max.apply(null, all);
    state.selection = { mode: 'protein', start: start, end: end,
                        resStart: first, resEnd: last };
    setInputKind('protein');
    paint();
    renderSelection();
  }

  function paint() {
    var selection = state.selection;
    Array.prototype.forEach.call(
      document.querySelectorAll('.b.sel, .r.sel, .b.caret, .b.pb, .r.pb'),
      function (node) {
        node.classList.remove('sel');
        node.classList.remove('caret');
        node.classList.remove('pb');
        for (var c = 0; c < PB_COLORS; c++) { node.classList.remove('pb' + c); }
        node.removeAttribute('data-pbname');
      });
    hideTip();
    paintImported();
    if (!selection) { return; }
    if (selection.mode === 'caret') {
      // Marks the base the insertion goes in front of.
      Array.prototype.forEach.call(
        document.querySelectorAll('.b[data-pos="' + selection.start + '"]'),
        function (node) { node.classList.add('caret'); });
      return;
    }
    /* Walked by length rather than start..end, because a selection found by
     * pasting a sequence can cross the origin and come back with end < start
     * (2026-09-28). A drag cannot produce that, which is why the plain loop
     * lasted this long. */
    var total = state.record.length();
    var span = ((selection.end - selection.start) % total + total) % total + 1;
    for (var step = 0; step < span; step++) {
      var p = (selection.start - 1 + step) % total + 1;
      Array.prototype.forEach.call(
        document.querySelectorAll('.b[data-pos="' + p + '"]'),
        function (node) { node.classList.add('sel'); });
    }
    /* Both tracks get painted whichever one the selection came from.
     * Painting only the protein track when the selection came from it left
     * the two views disagreeing: selecting residues highlighted the DNA, but
     * selecting DNA left the protein untouched. A residue is painted when the
     * selection covers any part of its codon. */
    var residues = {};
    if (selection.mode === 'protein') {
      for (var r = selection.resStart; r <= selection.resEnd; r++) {
        residues[r] = true;
      }
    } else if (state.residueAt) {
      for (var i = 0; i < span; i++) {
        var found = state.residueAt[(selection.start - 1 + i) % total + 1];
        if (found) { residues[found] = true; }
      }
    }
    Object.keys(residues).forEach(function (residue) {
      Array.prototype.forEach.call(
        document.querySelectorAll('.r[data-res="' + residue + '"]'),
        function (node) { node.classList.add('sel'); });
    });
  }

  /* The sequence picked from the Protein Builder list, coloured in both
   * tracks at once so the part to be changed is visible in either.
   *
   * Its place is worked out again on every paint instead of being kept as
   * coordinates: the plasmid, the reading frame and the start codon can all
   * change underneath it, and colour left over the old position would be a
   * lie. If it no longer sits in exactly one place, nothing is painted. */
  /* Each sequence keeps its own colour, and keeps it: the class comes from
   * where the sequence sits in the stored list, not from the order the
   * colours were switched on. Turning one off does not recolour the rest. */
  var PB_COLORS = 7;
  function colorClass(entry) {
    var at = state.pbuild.findIndex(function (kept) {
      return sameEntry(kept, entry);
    });
    return 'pb' + ((at < 0 ? 0 : at) % PB_COLORS);
  }

  function paintImported() {
    if (!state.pbHighlights.length || !state.protein || !state.map) { return; }
    // position -> {names, cls}. Where two sequences overlap the first one
    // keeps the colour and the tooltip names both.
    var atBase = {}, atResidue = {};
    function note(table, key, name, cls) {
      if (table[key]) { table[key].names += ' ／ ' + name; }
      else { table[key] = { names: name, cls: cls }; }
    }
    state.pbHighlights.forEach(function (entry) {
      var found = PBUILD.findIn(state.protein, entry.sequence);
      if (found.count !== 1) { return; }
      var cls = colorClass(entry);
      state.map.slice((found.resStart - 1) * 3, found.resEnd * 3)
        .forEach(function (position) {
          note(atBase, position, entry.name, cls);
        });
      for (var r = found.resStart; r <= found.resEnd; r++) {
        note(atResidue, r, entry.name, cls);
      }
    });
    function mark(selector, attribute, table) {
      Array.prototype.forEach.call(document.querySelectorAll(selector),
        function (node) {
          var hit = table[node.getAttribute(attribute)];
          if (!hit) { return; }
          node.classList.add('pb');
          node.classList.add(hit.cls);
          node.setAttribute('data-pbname', hit.names);
        });
    }
    mark('.seqview .b', 'data-pos', atBase);
    mark('.seqview .r', 'data-res', atResidue);
  }

  /* Where the three bases of a residue sit on the plasmid. They come out of
   * state.map in reading order, so a minus-strand CDS gives them descending
   * and a codon over the origin gives them out of one piece. */
  function codonText(residue) {
    if (!state.map) { return ''; }
    var positions = state.map.slice((residue - 1) * 3, residue * 3);
    if (positions.length !== 3) { return ''; }
    var sorted = positions.slice().sort(function (a, b) { return a - b; });
    return sorted[2] - sorted[0] === 2
      ? sorted[0] + '..' + sorted[2] + ' bp'
      : sorted.join(', ') + ' bp';
  }

  /* What the cursor is over: the imported sequence's name if there is one,
   * and always the numbers -- the residue and the bases it is made of, or
   * the plasmid position and the residue it belongs to. Asked for
   * 2026-09-27; counting letters by eye is what this replaces. */
  function tipText(target) {
    if (!target || !target.getAttribute) { return ''; }
    var parts = [];
    var name = target.getAttribute('data-pbname');
    if (name) { parts.push(name); }
    var residue = target.getAttribute('data-res');
    var position = target.getAttribute('data-pos');
    if (residue !== null) {
      var number = parseInt(residue, 10);
      var letter = (target.textContent || '').trim();
      parts.push('残基 ' + number + (letter ? ' ' + letter : ''));
      var bases = codonText(number);
      if (bases) { parts.push(bases); }
    } else if (position !== null) {
      parts.push(position + ' bp');
      var at = state.residueAt && state.residueAt[position];
      if (at) { parts.push('残基 ' + at); }
    }
    return parts.join(' ／ ');
  }

  /* A plain title= would have done, except that the sequence is one span per
   * base: moving along it leaves and enters elements constantly, and the
   * browser restarts its tooltip delay every time. This one appears at once
   * and follows the cursor. */
  function showTip(event) {
    var tip = $('pbTip');
    // Not while dragging out a selection: the range being chosen is the
    // thing to look at, and a box following the cursor sits on top of it.
    var text = state.dragging ? '' : tipText(event.target);
    if (!text) { tip.style.display = 'none'; return; }
    tip.textContent = text;
    tip.style.display = 'block';
    var room = document.documentElement.clientWidth - tip.offsetWidth - 8;
    tip.style.left = Math.max(4, Math.min(event.clientX + 12, room)) + 'px';
    tip.style.top = (event.clientY + 16) + 'px';
  }
  function hideTip() {
    var tip = $('pbTip');
    if (tip) { tip.style.display = 'none'; }
  }

  /* The edit is read off the selection and the field, so there is no
   * operation to pick by hand. Every case is distinguishable:
   *
   *   caret + text   -> insert, nothing removed
   *   caret + empty  -> nothing to do; runDesign says so
   *   range + text   -> remove the range, put the text in
   *   range + empty  -> remove the range
   */
  function inferredOperation() {
    var selection = state.selection;
    if (!selection) { return null; }
    var text = ($('newSequence').value || '').replace(/\s/g, '');
    if (selection.mode === 'caret') { return text ? 'insert' : null; }
    return text ? 'replace' : 'delete';
  }

  var OPERATION_LABEL = { insert: '挿入のみ', replace: '挿入+欠損',
                          'delete': '欠損のみ' };

  function renderSelection() {
    var selection = state.selection;
    /* The warning belongs where the work stops, not where the sequence is
     * shown: the field under it is the one that cannot be used yet
     * (2026-09-27). Set before anything returns. */
    /* No <strong> in a banner body: `.banner strong` is display:block for
     * the title, so any other one breaks the line. This notice read as six
     * lines because of it (2026-09-27). */
    $('editNotice').innerHTML = selection ? '' : banner('warn',
      '範囲が選ばれていません',
      '上の配列をなぞると範囲、クリックで挿入位置になります。');
    if (selection && selection.mode === 'caret') {
      $('selectionInfo').innerHTML = '<table><tr><th>挿入位置</th>'
        + '<td>プラスミド ' + selection.start + ' の直前</td></tr></table>'
        + '<div class="actions">'
        + '<button type="button" class="tiny" id="clearSelection">'
        + '選択を解除</button></div>';
      $('clearSelection').addEventListener('click', function () {
        state.selection = null; state.anchor = null;
        paint(); renderSelection();
      });
      describeField(selection.start + ' の直前に挿入します。何も除去しません。');
      return;
    }
    if (!selection) {
      $('selectionInfo').innerHTML = '';
      $('newSeqLabel').textContent = '新しい配列';
      $('newSeqHint').textContent = '';
      return;
    }
    var record = state.record;
    var dna = record.subseq(selection.start, selection.end);
    var rows = '<tr><th>選んだ範囲</th><td>プラスミド ' + selection.start + '..'
      + selection.end + '（' + dna.length + ' bp）</td></tr>';
    if (selection.mode === 'protein') {
      rows += '<tr><th>残基</th><td>' + selection.resStart + '..'
        + selection.resEnd + '（'
        + (selection.resEnd - selection.resStart + 1) + ' 残基）</td></tr>';
      var peptide = state.protein.slice(selection.resStart - 1, selection.resEnd);
      rows += '<tr><th>いまのアミノ酸</th><td class="seq">' + esc(peptide)
        + ' ' + copyBtn(peptide) + '</td></tr>';
    }
    rows += '<tr><th>いまの塩基配列</th><td class="seq">'
      + esc(dna.length > 120 ? dna.slice(0, 117) + '...' : dna)
      + ' ' + copyBtn(dna) + '</td></tr>';
    $('selectionInfo').innerHTML = '<table>' + rows + '</table>'
      + '<div class="actions">'
      + '<button type="button" class="tiny" id="clearSelection">選択を解除</button>'
      + '<span class="hint">選択の種類: '
      + (selection.mode === 'protein' ? 'タンパク質配列' : '塩基配列')
      + '</span></div>';
    $('clearSelection').addEventListener('click', function () {
      state.selection = null; state.pending = null;
      paint(); renderSelection();
    });
    Array.prototype.forEach.call(
      $('selectionInfo').querySelectorAll('[data-copy]'), function (button) {
        button.addEventListener('click', function () {
          copy(button.getAttribute('data-copy'), button);
        });
      });

    describeField('空欄にすると、この範囲の欠損になります。');
  }

  /* One place decides what the field is asking for, because the answer now
   * comes from the radio rather than from the track. The note about a
   * peptide written only in A/C/G/T goes here too: it cannot be refused
   * (those are real amino acids), but it is worth saying. */
  function describeField(tail) {
    var asProtein = isProteinInput(), caret = state.selection
      && state.selection.mode === 'caret';
    var text = ($('newSequence').value || '').replace(/\s/g, '');
    $('newSeqLabel').textContent = (caret ? '挿入する' : '新しい')
      + (asProtein ? 'アミノ酸配列' : '塩基配列');
    var hint = asProtein
      ? 'アミノ酸 1 文字表記で入力してください（例 GGSGGS）。'
        + 'コドンは宿主の使用頻度から選びます。'
      : 'A/T/G/C で入力してください。';
    var problem = inputError(text, asProtein);
    if (problem) { hint = problem; }
    else if (asProtein && text.length >= 6 && DNA_ONLY.test(text.toUpperCase())) {
      hint += ' ※ A/T/G/C だけの配列です。塩基配列のつもりなら'
        + '「塩基配列で入力」を選んでください。';
    }
    $('newSeqHint').textContent = hint + (problem ? '' : ' ' + tail);
  }

  // ---- feature editing -------------------------------------------------

  /* The sequences are kept in this browser, so the backup is read once
   * rather than every visit (asked for 2026-09-27). Same store as the saved
   * designs: localStorage, no network, gone when site data is cleared. */
  var PROTEIN_KEY = 'primerBuilder.proteins.v1';

  // null (not []) when the browser refuses storage, so "none kept" and
  // "cannot keep" stay tellable apart.
  function readProteins() {
    try {
      var raw = window.localStorage.getItem(PROTEIN_KEY);
      var list = raw ? JSON.parse(raw) : [];
      return Object.prototype.toString.call(list) === '[object Array]'
        ? list : [];
    } catch (error) { return null; }
  }

  function writeProteins(list) {
    try {
      window.localStorage.setItem(PROTEIN_KEY, JSON.stringify(list));
      return true;
    } catch (error) {
      say('保存できませんでした（ブラウザの保存領域が足りません）。');
      return false;
    }
  }

  function sameEntry(a, b) {
    return !!a && !!b && a.name === b.name && a.sequence === b.sequence;
  }

  function isHighlighted(entry) {
    return state.pbHighlights.some(function (kept) {
      return sameEntry(kept, entry);
    });
  }

  function toggleHighlight(entry) {
    state.pbHighlights = isHighlighted(entry)
      ? state.pbHighlights.filter(function (kept) {
          return !sameEntry(kept, entry);
        })
      : state.pbHighlights.concat([entry]);
    renderLibraryList();
    paint();
  }

  function dropProtein(index) {
    var stored = readProteins();
    var going = state.pbuild[index];
    state.pbHighlights = state.pbHighlights.filter(function (kept) {
      return !sameEntry(kept, going);
    });
    if (stored === null) {
      state.pbuild.splice(index, 1);
    } else {
      var entry = state.pbuild[index];
      writeProteins(stored.filter(function (kept) {
        return !(kept.name === entry.name && kept.sequence === entry.sequence);
      }));
      state.pbuild = readProteins() || [];
    }
    renderLibraryList();
  }

  /* The list says up front whether each sequence is in the CDS on screen,
   * so nobody has to press a button to find out. */
  /* The count goes in the fold's own summary, so it still says something
   * useful while the list is closed (asked for 2026-09-27). */
  function setBuildCount(kept, here) {
    var label = $('pbuildCount');
    if (!label) { return; }
    label.textContent = !kept ? ''
      : '（' + kept + ' 件' + (here ? '・うち ' + here + ' 件がこの CDS に' : '')
        + '）';
  }

  function renderLibraryList() {
    var box = $('pbuildList');
    if (!box) { return; }
    if (!state.pbuild.length) {
      setBuildCount(0, 0);
      box.innerHTML = '<p class="hint">まだありません。'
        + '下で .json を読み込むと、ここに残ります。</p>';
      return;
    }
    var rows = PBUILD.arrange(state.protein || '', state.pbuild);
    setBuildCount(state.pbuild.length, rows.filter(function (row) {
      return row.found.count === 1;
    }).length);
    /* Colour only. Writing these into the plasmid as features went with the
     * rest of the feature editor on 2026-09-27: 「ProteinBuilder から
     * タンパク質配列を選んでハイライトが表示できるならそれが一番です」. */
    var html = '<div class="scroll"><table><tr><th>名前</th><th>残基</th>'
      + '<th>いまの CDS の中</th><th>画面の色</th><th>一覧から</th></tr>';
    rows.forEach(function (row) {
      var found = row.found;
      var where = !state.protein ? '読み枠が出ていません'
        : found.count === 1 ? '残基 ' + found.resStart + '..' + found.resEnd
        : found.count === 0 ? 'なし' : found.count + ' 箇所';
      html += '<tr><td>' + esc(row.entry.name) + '</td><td>'
        + row.entry.sequence.length + '</td><td>' + esc(where) + '</td><td>'
        + (found.count === 1
            ? '<span class="chip ' + colorClass(row.entry) + '"></span>'
              + '<button type="button" class="tiny" data-pbcolor="' + row.index
              + '">' + (isHighlighted(row.entry) ? '消す' : 'つける')
              + '</button>' : '')
        + '</td><td><button type="button" class="tiny" data-pbdel="'
        + row.index + '">削除</button></td></tr>';
    });
    html += '</table></div><p class="hint">'
      + '色は配列ごとに違います。画面で見るためだけのもので、'
      + 'ページを閉じると消え、書き出すファイルには入りません。</p>';
    box.innerHTML = html;
    // The assembly step, when this build has one, picks from the same list.
    if (window.PBASM) { window.PBASM.refresh(); }
    Array.prototype.forEach.call(box.querySelectorAll('[data-pbcolor]'),
      function (button) {
        button.addEventListener('click', function () {
          toggleHighlight(state.pbuild[
            parseInt(button.getAttribute('data-pbcolor'), 10)]);
        });
      });
    Array.prototype.forEach.call(box.querySelectorAll('[data-pbdel]'),
      function (button) {
        button.addEventListener('click', function () {
          dropProtein(parseInt(button.getAttribute('data-pbdel'), 10));
        });
      });
  }

  function readBuildFile(file) {
    var reader = new FileReader();
    reader.onload = function (event) {
      var parsed;
      try { parsed = PBUILD.parseBackup(String(event.target.result)); }
      catch (error) { say(error.message); return; }
      var stored = readProteins();
      var merged = PBUILD.merge(stored === null ? state.pbuild : stored,
                                parsed.sequences);
      var kept = stored !== null && writeProteins(merged.list);
      state.pbuild = merged.list;
      renderLibraryList();
      var note = parsed.warnings.join(' ');
      say((merged.added ? merged.added + ' 件を読み込みました。'
                        : '新しい配列はありませんでした。')
          + (kept ? '' : ' このブラウザでは保存できないので、'
                         + 'このページを閉じるまでの間だけ使えます。')
          + (note ? ' ' + note : ''),
          kept && !note ? 'ok' : 'warn');
    };
    reader.onerror = function () { say('ファイルを読めませんでした。'); };
    reader.readAsText(file);
  }

  // ---- design ----------------------------------------------------------

  function runDesign() {
    // Drop the previous result first. Every early return below would
    // otherwise leave the last design on screen while the message says the
    // new one could not be made, which reads as if it had succeeded.
    state.result = null;
    $('result').innerHTML = '';

    var selection = state.selection;
    if (!selection) {
      say('先に上の配列で範囲を選んでください。');
      $('viewSection').scrollIntoView({ behavior: 'smooth', block: 'start' });
      return;
    }
    var newSequence = ($('newSequence').value || '').replace(/\s/g, '');
    var asProtein = isProteinInput();
    var problem = inputError(newSequence, asProtein);
    if (problem) { say(problem); return; }
    var operation = inferredOperation();
    if (!operation) {
      say('挿入する配列を入力してください。');
      return;
    }
    var start = selection.start, end = selection.end;
    if (operation === 'insert') {
      // end == start - 1 is how the engine spells "insert, remove nothing".
      // At the origin that would be position 0, so it wraps to the last base.
      end = start === 1 ? state.record.length() : start - 1;
    }

    var options;
    try {
      options = {
        start: start, end: end, newSequence: newSequence,
        // Chosen in the field above, not taken from the track the selection
        // came from, and it applies to insertions as well (2026-09-27).
        // Kept true for an empty field too: a deletion asked for in amino
        // acids still has to start and end on codon boundaries.
        asProtein: asProtein,
        /* The ORF found between the T7 promoter and terminator is the frame
         * of reference. Passing it lets the engine confirm that the edit
         * starts on a codon boundary instead of warning that it could not
         * check. Selections made in the protein track are codon-aligned by
         * construction, but the engine cannot see that from coordinates
         * alone. A Feature or a location string are both accepted. */
        cds: state.orf || null,
        maxPrimerLen: num('maxPrimerLen', PB.DEFAULT_MAX_PRIMER_LEN),
        minAnnealLen: num('minAnnealLen', PB.DEFAULT_MIN_ANNEAL_LEN),
        maxAnnealLen: num('maxAnnealLen', PB.DEFAULT_MAX_ANNEAL_LEN),
        // The field is a percentage; the engine wants a fraction.
        minUsage: num('minUsage', 80) / 100,
        candidates: 5
      };
    } catch (error) { say(error.message); return; }

    try {
      state.result = PB.designReplacement(state.record, options);
    } catch (error) {
      state.result = null;
      $('result').innerHTML = '';
      say(error.message);
      return;
    }
    say('');
    renderResult();
    $('result').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  /* No note steering the user to another method here. The tool used to say
   * that a pure insertion or deletion did not need inverse PCR and that
   * QuikChange would do; removed on 2026-09-08 on the user's instruction.
   * The method is their decision, not something the tool should reopen. */

  /* The PCR conditions as a blank form.
   *
   * Until 2026-09-25 the screen printed the temperatures, times and cycle
   * count, and a protocol underneath it. The user asked for both to go: they
   * know the protocol, and they want to work the conditions out themselves.
   * What the tool owes them is the primers. The rows stay so there is
   * somewhere to write the answer. */
  function pcrForm() {
    return '<h3>PCR 条件（記入欄）</h3><table>'
      + ['熱変性', 'アニーリング', '伸長', 'サイクル数'].map(function (stage) {
        return '<tr><th>' + stage
          + '</th><td><input type="text" class="fill"></td></tr>';
      }).join('')
      + '</table>';
  }

  function renderResult() {
    var result = state.result;
    var html = '';

    var verified = result.verification && !result.verification.problems.length;
    html += verified
      ? banner('ok', '設計と検算が通りました',
               'プライマーの配列から環状化後のプラスミドを組み立て直し、'
               + '狙った配列になることを確認しました。')
      : banner('bad', '問題があります',
               (result.problems || [])
                 .concat(result.verification ? result.verification.problems : [])
                 .map(esc).join('<br>'));
    (result.notes || [])
      .concat(result.verification ? (result.verification.notes || []) : [])
      .forEach(function (n) {
        html += banner('warn', '注意', esc(n));
      });

    /* No pair means the search found nothing that fits -- most often the
     * insert is longer than the two 5' tails can carry. Everything below
     * reads result.pair, so it has to stop here. It did not, and the page
     * threw instead of saying so (found 2026-09-28 while assembling a 38
     * residue region: 114 bp of insert against a 84 bp budget). */
    if (!result.pair) {
      html += '<table>'
        + '<tr><th>置き換える領域</th><td>プラスミド ' + result.start + '..'
        + result.end + '（' + result.site.removedLength + ' bp を除去）</td></tr>'
        + '<tr><th>入れる配列</th><td>' + result.site.insertLength
        + ' bp</td></tr>'
        + '<tr><th>テールに載せられる上限</th><td>' + result.budget
        + ' bp</td></tr></table>'
        + '<p class="hint">プライマー長の上限を上げるか、アニール域を短くするか、'
        + '入れる配列を短くしてください。'
        + '上限を上げると合成の収量と純度が落ちます。</p>';
      $('result').innerHTML = html;
      return;
    }

    html += '<table>'
      + '<tr><th>置き換える領域</th><td>プラスミド ' + result.start + '..'
      + result.end + '（' + result.site.removedLength + ' bp を除去）</td></tr>'
      + '<tr><th>入れる配列</th><td class="seq">' + result.site.insertLength
      + ' bp ' + esc(result.site.insert) + ' '
      + copyBtn(result.site.insert) + '</td></tr>'
      + '<tr><th>予想産物長</th><td>' + result.pair.productLength
      + ' bp</td></tr></table>';

    html += '<h3>プライマー</h3>'
      + '<div class="scroll"><table><tr><th></th><th>配列 (5′→3′)</th>'
      + '<th>長さ</th><th>Tm</th><th>GC</th><th>結合部位</th><th></th></tr>';
    [[result.forwardName || 'F', result.pair.forward],
     [result.reverseName || 'R', result.pair.reverse]].forEach(function (entry) {
      var primer = entry[1];
      html += '<tr><td class="mono">' + esc(entry[0]) + '</td>'
        + '<td class="seq">' + esc(primer.sequence) + '</td><td>'
        + primer.length + ' nt</td><td>' + primer.tm + '℃</td><td>'
        + primer.gcPercent.toFixed(0) + '%</td><td>' + primer.bindingSites
        + ' 箇所' + (primer.unique ? '' : ' ⚠') + '</td><td>'
        + copyBtn(primer.sequence) + '</td></tr>';
    });
    html += '</table></div>';
    html += '<p>' + copyBtn(
      (result.forwardName || 'F') + '\t' + result.pair.forward.sequence + '\n'
      + (result.reverseName || 'R') + '\t' + result.pair.reverse.sequence,
      '2 本まとめてコピー') + '</p>';

    html += pcrForm();

    if (verified) {
      /* Two different .ape files come out of this page, so they must not
       * share a name: they land in the same download folder. */
      // The marker this page just added is not one of the parent's, so it
      // is not counted among the features carried over.
      var carried = state.result.verification.newRecord.features.length
        - (state.result.verification.editMarker ? 1 : 0);
      html += '<p class="hint">書き出される <code>〜_designed.ape</code> は'
        + '<strong>塩基配列を改変したプラスミド</strong>です。'
        + 'フィーチャー（注釈）は座標を付け替えて引き継がれます（'
        + carried + ' 件）。編集部位にかかるものだけは外し、上に理由を出します。'
        + '改変した範囲には <code>' + esc(PB.EDIT_MARKER_LABEL)
        + '</code> の注釈が付きます。</p>'
        + '<div class="actions">'
        + '<button type="button" class="btn" id="downloadResult">'
        + '改変後のプラスミド (.ape) を書き出す</button>'
        + '<button type="button" class="btn" id="saveDesign">'
        + 'この設計を保存する</button></div>';
    }
    $('result').innerHTML = html;
    Array.prototype.forEach.call(
      $('result').querySelectorAll('[data-copy]'), function (button) {
        button.addEventListener('click', function () {
          copy(button.getAttribute('data-copy'), button);
        });
      });
    var download = $('downloadResult');
    if (download) {
      download.addEventListener('click', function () {
        saveText(fileBase() + '_designed.ape', designedText());
      });
      $('saveDesign').addEventListener('click', function () {
        var result = state.result;
        saveDesign({
          name: state.record.name || 'plasmid',
          base: fileBase(),
          summary: result.start + '..' + result.end + ' を '
            + result.site.insertLength + ' bp に置き換え',
          report: PB.replacementMarkdown(result),
          template: PB.formatGenbank(state.record),
          edited: designedText()
        });
      });
    }
  }

  function plasmidBase() {
    return (state.record.name || 'plasmid').replace(/[^\w.\-]/g, '_');
  }
  // The date goes in front of the file name only, never into the record's
  // own name: the LOCUS line is what ApE shows (2026-09-25).
  function fileBase() { return stamp() + '_' + plasmidBase(); }
  function designedText() {
    var record = state.result.verification.newRecord;
    record.name = plasmidBase() + '_designed';
    return PB.formatGenbank(record);
  }

  // ---- copy and save ---------------------------------------------------

  function copy(text, button) {
    function done() {
      var original = button.textContent;
      button.textContent = 'コピーしました';
      setTimeout(function () { button.textContent = original; }, 1200);
    }
    if (navigator.clipboard && window.isSecureContext) {
      navigator.clipboard.writeText(text).then(done, function () {
        legacyCopy(text, done);
      });
    } else { legacyCopy(text, done); }
  }
  function legacyCopy(text, done) {
    var area = document.createElement('textarea');
    area.value = text;
    area.style.position = 'fixed'; area.style.opacity = '0';
    document.body.appendChild(area);
    area.select();
    try { document.execCommand('copy'); done(); }
    catch (error) { say('コピーできませんでした。手動で選択してください。'); }
    document.body.removeChild(area);
  }
  /* Files land in one download folder and stay there for years. The date
   * in front is how the user tells this design from the next one; it is the
   * project's own file-naming rule (YYYYMMDD_name). Asked for 2026-09-25. */
  function stamp() {
    var d = new Date(), m = d.getMonth() + 1, day = d.getDate();
    return String(d.getFullYear()) + (m < 10 ? '0' + m : m)
      + (day < 10 ? '0' + day : day);
  }
  function saveText(name, text) {
    try {
      var blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
      var url = URL.createObjectURL(blob);
      var link = document.createElement('a');
      link.href = url; link.download = name;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
    } catch (error) {
      $('textPanel').innerHTML = '<h3>' + esc(name) + '</h3><p>'
        + copyBtn(text, '全文をコピー') + '</p><textarea rows="14" readonly>'
        + esc(text) + '</textarea>';
      Array.prototype.forEach.call(
        $('textPanel').querySelectorAll('[data-copy]'), function (button) {
          button.addEventListener('click', function () {
            copy(button.getAttribute('data-copy'), button);
          });
        });
    }
  }

  // ---- saved designs ---------------------------------------------------

  /* The primers, the template and the edited plasmid, kept in the browser.
   *
   * Asked for on 2026-09-25: a record of what was made, so a design can be
   * found again months later. It lives in this browser's localStorage and
   * goes nowhere else -- nothing on this page touches the network -- so
   * clearing the browser's site data clears it too.
   *
   * Each entry holds the three files as text rather than a reference to
   * anything on disk, because the point is to still have them when the
   * original files have been moved or renamed. That makes an entry roughly
   * the size of two plasmids, so HISTORY_LIMIT caps how many are kept. */
  var HISTORY_KEY = 'primerBuilder.history.v1';
  var HISTORY_LIMIT = 20;

  // null (not []) when the browser refuses storage, so the difference
  // between "no saved designs" and "cannot save" stays visible.
  function readHistory() {
    try {
      var raw = window.localStorage.getItem(HISTORY_KEY);
      var list = raw ? JSON.parse(raw) : [];
      return Object.prototype.toString.call(list) === '[object Array]'
        ? list : [];
    } catch (error) { return null; }
  }

  function timestamp() {
    var d = new Date();
    function two(n) { return n < 10 ? '0' + n : String(n); }
    return d.getFullYear() + '-' + two(d.getMonth() + 1) + '-'
      + two(d.getDate()) + ' ' + two(d.getHours()) + ':' + two(d.getMinutes());
  }

  function saveDesign(entry) {
    var list = readHistory();
    if (list === null) {
      say('このブラウザは保存領域を使えない設定になっています。'
          + '設計は保存できません。');
      return;
    }
    entry.id = String(Date.now()) + Math.random().toString(36).slice(2, 8);
    entry.savedAt = timestamp();
    list.unshift(entry);
    if (list.length > HISTORY_LIMIT) { list = list.slice(0, HISTORY_LIMIT); }
    // Out of room: drop the oldest and try again, so a save either happens
    // or says why -- never quietly does nothing.
    while (true) {
      try {
        window.localStorage.setItem(HISTORY_KEY, JSON.stringify(list));
        break;
      } catch (error) {
        if (list.length <= 1) {
          say('保存できませんでした（ブラウザの保存領域が足りません）。');
          return;
        }
        list = list.slice(0, list.length - 1);
      }
    }
    renderHistory();
    say('保存しました。下の「保存した設計」から取り出せます。', 'ok');
  }

  function dropDesign(id) {
    var list = readHistory();
    if (list === null) { return; }
    window.localStorage.setItem(HISTORY_KEY, JSON.stringify(
      list.filter(function (entry) { return entry.id !== id; })));
    renderHistory();
  }

  var HISTORY_FILES = [
    ['report', 'プライマー', '_design.md'],
    ['template', '鋳型', '_template.ape'],
    ['edited', '改変後', '_designed.ape']
  ];

  function renderHistory() {
    var box = $('historyList'), list = readHistory();
    if (list === null) {
      box.innerHTML = '<p class="hint">このブラウザでは保存できません。</p>';
      return;
    }
    if (!list.length) { box.innerHTML = '<p class="hint">まだありません。</p>'; return; }
    var html = '<div class="scroll"><table><tr><th>保存日時</th><th>名前</th>'
      + '<th>内容</th><th>取り出す</th><th></th></tr>';
    list.forEach(function (entry) {
      html += '<tr><td>' + esc(entry.savedAt) + '</td><td>' + esc(entry.name)
        + '</td><td>' + esc(entry.summary) + '</td><td>'
        + HISTORY_FILES.map(function (file) {
          return '<button type="button" class="tiny" data-take="' + file[0]
            + '" data-id="' + esc(entry.id) + '">' + file[1] + '</button>';
        }).join(' ')
        + '</td><td><button type="button" class="tiny" data-drop="'
        + esc(entry.id) + '">削除</button></td></tr>';
    });
    box.innerHTML = html + '</table></div>';

    Array.prototype.forEach.call(box.querySelectorAll('[data-take]'),
      function (button) {
        button.addEventListener('click', function () {
          var entry = (readHistory() || []).filter(function (e) {
            return e.id === button.getAttribute('data-id');
          })[0];
          if (!entry) { renderHistory(); return; }
          var kind = button.getAttribute('data-take');
          HISTORY_FILES.forEach(function (file) {
            if (file[0] === kind) { saveText(entry.base + file[2], entry[kind]); }
          });
        });
      });
    Array.prototype.forEach.call(box.querySelectorAll('[data-drop]'),
      function (button) {
        button.addEventListener('click', function () {
          if (window.confirm('この設計を消しますか。元に戻せません。')) {
            dropDesign(button.getAttribute('data-drop'));
          }
        });
      });
  }

  // ---- requests --------------------------------------------------------

  var requests = [];
  function today() {
    var d = new Date(), m = d.getMonth() + 1, day = d.getDate();
    return d.getFullYear() + '-' + (m < 10 ? '0' + m : m) + '-'
      + (day < 10 ? '0' + day : day);
  }
  function currentSettings() {
    var parts = [];
    function add(label, value) { if (value) { parts.push(label + ': ' + value); } }
    add('鋳型の長さ', state.record ? state.record.length() + ' bp' : '');
    add('モード', isMutationMode() ? '変異表記' : '範囲選択');
    if (isMutationMode()) {
      PBMUT.settings().forEach(function (line) { parts.push(line); });
    }
    /* A caret has a start and no end, so the old `start..end` printed
     * `65..undefined` in the request form (seen 2026-09-27). */
    add('選択範囲', !isMutationMode() && state.selection
        ? (state.selection.mode === 'caret'
             ? state.selection.start + ' の直前（挿入位置）'
             : state.selection.start + '..' + state.selection.end
               + '（'
               + (state.selection.mode === 'protein' ? 'タンパク質' : '塩基')
               + 'から選択）')
        : '');
    add('操作', !isMutationMode()
        ? (OPERATION_LABEL[inferredOperation()] || '') : '');
    add('プライマー長の上限', ($('maxPrimerLen').value || '') + ' nt');
    add('アニール域', ($('minAnnealLen').value || '') + '-'
        + ($('maxAnnealLen').value || '') + ' nt');
    add('使用頻度の下限', ($('minUsage').value || '') + ' %');
    if (state.result) {
      add('結果', state.result.verification
          && !state.result.verification.problems.length
          ? '検算まで通った' : '問題が残っている');
    }
    return parts;
  }
  function addRequest() {
    var text = ($('requestText').value || '').trim();
    if (!text) {
      $('requestNotice').innerHTML = banner('warn', 'リクエストが空です',
        '内容を入力してから「追加する」を押してください。');
      return;
    }
    requests.push({ date: today(), text: text,
                    settings: $('attachSettings').checked
                      ? currentSettings() : null });
    $('requestText').value = '';
    $('requestNotice').innerHTML = '';
    renderRequests();
  }
  function requestsMarkdown() {
    var lines = ['# Primer Builder へのリクエスト', '',
      'このファイルは Primer Builder.html のリクエスト欄から書き出したものです。',
      '配列データは含まれていません。', ''];
    requests.forEach(function (entry, index) {
      lines.push('## #' + (index + 1) + ' - ' + entry.date + ' - 未対応', '');
      lines.push(entry.text, '');
      if (entry.settings && entry.settings.length) {
        lines.push('そのときの設定:', '');
        entry.settings.forEach(function (line) { lines.push('- ' + line); });
        lines.push('');
      }
    });
    return lines.join('\n');
  }
  function renderRequests() {
    if (!requests.length) { $('requestList').innerHTML = ''; return; }
    var html = banner('warn', 'このページは何も保存しません',
      'タブを閉じると下の内容は消えます。'
      + '「まとめてコピー」か「REQUESTS.md として保存」で取り出してください。');
    html += '<div class="scroll"><table><tr><th>#</th><th>日付</th>'
      + '<th>内容</th><th>設定</th><th></th></tr>';
    requests.forEach(function (entry, index) {
      html += '<tr><td>' + (index + 1) + '</td><td>' + esc(entry.date)
        + '</td><td>' + esc(entry.text) + '</td><td>'
        + (entry.settings ? '添付あり' : '-') + '</td><td>'
        + '<button type="button" class="tiny" data-remove="' + index
        + '">削除</button></td></tr>';
    });
    html += '</table></div><div class="actions">'
      + copyBtn(requestsMarkdown(), 'まとめてコピー')
      + '<button type="button" class="btn" id="saveRequests">'
      + 'REQUESTS.md として保存</button></div>';
    $('requestList').innerHTML = html;
    Array.prototype.forEach.call(
      $('requestList').querySelectorAll('[data-copy]'), function (button) {
        button.addEventListener('click', function () {
          copy(button.getAttribute('data-copy'), button);
        });
      });
    Array.prototype.forEach.call(
      $('requestList').querySelectorAll('[data-remove]'), function (button) {
        button.addEventListener('click', function () {
          requests.splice(parseInt(button.getAttribute('data-remove'), 10), 1);
          renderRequests();
        });
      });
    $('saveRequests').addEventListener('click', function () {
      saveText('REQUESTS.md', requestsMarkdown());
    });
  }

  // ---- self test -------------------------------------------------------

  function runSelfTests() {
    var passed = 0, failures = [];
    PBTESTS.cases.forEach(function (item) {
      try { item.fn(); passed++; }
      catch (error) { failures.push([item.name, error.message]); }
    });
    var html = failures.length
      ? banner('bad', passed + ' / ' + PBTESTS.cases.length + ' 通過（'
               + failures.length + ' 件失敗）',
               'このツールの計算は信用できません。')
      : banner('ok', passed + ' / ' + PBTESTS.cases.length + ' すべて通過',
               'Tm・コドン選択・座標計算・プライマー構成・検算・T7 の認識が、'
               + 'すべて手計算した期待値と一致しました。');
    html += '<div class="scroll"><table><tr><th></th><th>項目</th>'
      + '<th>結果</th></tr>';
    PBTESTS.cases.forEach(function (item) {
      var failure = null;
      failures.forEach(function (f) { if (f[0] === item.name) { failure = f[1]; } });
      html += '<tr><td>' + (failure ? '❌' : '✅') + '</td><td>' + esc(item.name)
        + '</td><td>' + (failure ? esc(failure) : 'OK') + '</td></tr>';
    });
    $('testResult').innerHTML = html + '</table></div>';
  }

  /* Two examples, because the two modes need different things from the
   * plasmid: the mutation mode numbers residues against a CDS annotation,
   * and the range-selection mode looks for the T7 promoter and terminator. */
  function loadExample() {
    var mutation = isMutationMode();
    state.text = mutation ? PBTESTS.EXAMPLE_GENBANK : PBTESTS.T7_EXAMPLE;
    state.filename = mutation ? '20260826_example-plasmid.gb' : 'T7-example.ape';
    $('filename').textContent = mutation
      ? '同梱の合成プラスミド（CDS 注釈つき）'
      : 'T7 promoter/terminator つきの合成プラスミド';
    $('pasted').value = '';
    if (mutation) { $('mutations').value = 'D2A,ins@6:GGSGGS'; }
    loadSequence();
  }

  // Shared with the mutation mode, which lives in its own unit.
  window.PBUI = {
    copy: copy, saveText: saveText, pcrForm: pcrForm, stamp: stamp,
    saveDesign: saveDesign,
    // Read by assemble.js in the combined build; harmless without it.
    proteins: function () { return state.pbuild.slice(); },
    usePeptide: function (peptide) {
      $('newSequence').value = peptide;
      setInputKind('protein');
      renderSelection();
      $('editSection').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  };

  // ---- wiring ----------------------------------------------------------

  document.addEventListener('DOMContentLoaded', function () {
    var drop = $('drop'), picker = $('picker');
    showPasteArea(true);
    state.pbuild = readProteins() || [];
    renderLibraryList();
    renderHistory();
    picker.addEventListener('change', function () {
      if (picker.files.length) { readFile(picker.files[0]); }
    });
    ['dragenter', 'dragover'].forEach(function (name) {
      drop.addEventListener(name, function (event) {
        event.preventDefault(); drop.classList.add('over');
      });
    });
    ['dragleave', 'drop'].forEach(function (name) {
      drop.addEventListener(name, function (event) {
        event.preventDefault(); drop.classList.remove('over');
      });
    });
    drop.addEventListener('drop', function (event) {
      if (event.dataTransfer.files.length) { readFile(event.dataTransfer.files[0]); }
    });
    $('loadPasted').addEventListener('click', function () {
      state.text = ''; state.filename = '';
      $('filename').textContent = '貼り付けたテキスト';
      loadSequence();
    });
    var buildPicker = $('pbuildFile');
    buildPicker.addEventListener('change', function () {
      if (buildPicker.files.length) { readBuildFile(buildPicker.files[0]); }
    });
    $('findButton').addEventListener('click', findPasted);
    $('showPaste').addEventListener('click', function () {
      showPasteArea(true);
      $('pasted').focus();
    });
    $('exampleButton').addEventListener('click', loadExample);
    $('resetDesign').addEventListener('click', resetDesign);
    $('designButton').addEventListener('click', function () {
      if (isMutationMode()) { PBMUT.design(null); } else { runDesign(); }
    });
    Array.prototype.forEach.call(
      document.querySelectorAll('input[name=mode]'), function (input) {
        input.addEventListener('change', applyMode);
      });
    $('testButton').addEventListener('click', runSelfTests);
    $('addRequest').addEventListener('click', addRequest);
    $('rescan').addEventListener('click', function () {
      if (state.record) { renderViews(); renderSelection(); }
    });
    $('resetDefaults').addEventListener('click', function () {
      ['maxPrimerLen', 'minAnnealLen', 'maxAnnealLen', 'minUsage']
        .forEach(function (id) {
          $(id).value = $(id).getAttribute('data-default');
        });
      $('diagnostic').value = '';
      $('t7promoter').value = GEVIEW.T7_PROMOTER;
      $('t7terminator').value = GEVIEW.T7_TERMINATOR;
    });
    wireNewSequenceEcho();
    document.addEventListener('mouseup', function () { state.dragging = false; });
    $('t7promoter').value = GEVIEW.T7_PROMOTER;
    $('t7terminator').value = GEVIEW.T7_TERMINATOR;
    hideAllSteps();
  });
})();
